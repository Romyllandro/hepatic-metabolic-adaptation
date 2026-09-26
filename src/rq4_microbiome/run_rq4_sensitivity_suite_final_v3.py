#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

def tag(x): return str(x).replace(".","p")

def norm_species(x):
    import re
    return re.sub(r'[^a-z0-9]+','_',str(x).lower()).strip('_')


def run(cmd,log,resume=False,expected=None,dry=False):
    log=Path(log); log.parent.mkdir(parents=True,exist_ok=True)
    if resume and expected and Path(expected).exists():
        print(f"[RESUME] exists -> {expected}")
        return
    print("\n> "+" ".join(map(str,cmd)))
    if dry: return
    with log.open("w",encoding="utf-8") as f:
        p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,text=True)
    if p.returncode:
        tail=log.read_text(errors="replace").splitlines()[-120:]
        print("\n".join(tail))
        raise RuntimeError(f"Command failed: {' '.join(map(str,cmd))}")

def included_map(final_map,out,expected_n):
    d=pd.read_csv(final_map)
    inc=d[d["included"].astype(str).str.lower().isin(["true","1","yes"])].copy()
    inc=inc[inc["model_file"].notna() & (inc["model_file"].astype(str).str.strip()!="")]
    if len(inc)!=expected_n:
        raise ValueError(f"{final_map}: expected {expected_n} included species, found {len(inc)}")
    inc=inc.copy()
    inc["species"]=inc["species"].map(norm_species)
    if inc["species"].duplicated().any():
        dup=inc.loc[inc["species"].duplicated(keep=False),"species"].tolist()
        raise ValueError(f"Duplicate normalized species keys in included map: {dup}")
    inc[["species","model_file"]].to_csv(out,index=False)
    return len(inc)

def validate_community_output(cdir, expression_cols, expected_mapped):
    cdir=Path(cdir)
    mp=cdir/"species_to_model_mapping.csv"
    if not mp.exists():
        raise RuntimeError(f"Missing community mapping output: {mp}")
    m=pd.read_csv(mp)
    n=int(m["match_status"].isin(["matched","matched_provisional"]).sum())
    if n != int(expected_mapped):
        raise RuntimeError(
            f"Community mapping QC failed in {cdir}: mapped {n}, expected {expected_mapped}."
        )
    for cond in expression_cols:
        tax=cdir/f"taxonomy_{cond}.csv"
        cov=cdir/f"mapping_coverage_{cond}.json"
        raw=cdir/f"community_solution_fluxes_{cond}.csv"
        ex=cdir/f"community_exchange_fluxes_{cond}.csv"
        missing=[str(p) for p in (tax,cov,raw,ex) if not p.exists()]
        if missing:
            raise RuntimeError(
                f"MICOM community did not complete for {cond}. Missing: {missing}"
            )
        covj=json.loads(cov.read_text())
        retained=float(covj.get("retained_activity_weighted_abundance_before_renormalization",0))
        if retained <= 0:
            raise RuntimeError(f"Retained activity coverage is zero for {cond}.")
        print(
            f"[COMMUNITY QC] {cond}: modeled={covj.get('n_species_modeled')}, "
            f"retained_activity={retained:.3%}, active_exchange_file={ex.name}"
        )

def scenario_list(tradeoffs,scales):
    out=[
        {"name":"primary","tradeoff":0.5,"objective":"biomass","scale":0.1,"portal_mode":"forced"},
        {"name":"objective_atpm","tradeoff":0.5,"objective":"atpm","scale":0.1,"portal_mode":"forced"},
        {"name":"portal_mode_soft","tradeoff":0.5,"objective":"biomass","scale":0.1,"portal_mode":"soft"},
        {"name":"portal_mode_availability","tradeoff":0.5,"objective":"biomass","scale":0.1,"portal_mode":"availability"},
    ]
    for s in scales:
        if abs(s-0.1)>1e-12:
            out.append({"name":f"portal_scale_{tag(s)}","tradeoff":0.5,
                        "objective":"biomass","scale":s,"portal_mode":"forced"})
    for t in tradeoffs:
        if abs(t-0.5)>1e-12:
            out.append({"name":f"micom_tradeoff_{tag(t)}","tradeoff":t,
                        "objective":"biomass","scale":0.1,"portal_mode":"forced"})
    seen=set(); uniq=[]
    for x in out:
        k=(x["tradeoff"],x["objective"],x["scale"],x["portal_mode"])
        if k not in seen:
            seen.add(k); uniq.append(x)
    return uniq

def summarize(path):
    d=pd.read_csv(path)
    vm="variance_explained_microbiome"; vd="variance_explained_diet"; dom="dominant_driver"
    measurable=pd.to_numeric(d[vm],errors="coerce").fillna(0)>0
    return {
        "n_reactions":len(d),
        "n_measurable_microbiome":int(measurable.sum()),
        "n_microbiome_dominant":int(d[dom].astype(str).str.startswith("Microbiome").sum()),
        "n_stable":int(d[dom].astype(str).eq("Stable").sum()),
        "mean_ve_microbiome":float(pd.to_numeric(d[vm],errors="coerce").mean()),
        "mean_ve_diet":float(pd.to_numeric(d[vd],errors="coerce").mean()),
    }, d

def compare(ref,q,name):
    rows=[]; q=q.reindex(ref.index)
    for metric in ["delta_diet","delta_microbiome","variance_explained_diet","variance_explained_microbiome"]:
        if metric in ref.columns and metric in q.columns:
            a=pd.to_numeric(ref[metric],errors="coerce"); b=pd.to_numeric(q[metric],errors="coerce")
            ok=a.notna()&b.notna()
            rho=float(spearmanr(a[ok],b[ok]).correlation) if ok.sum()>=3 else np.nan
            rows.append({"scenario":name,"metric":metric,"spearman_vs_primary":rho})
    if "dominant_driver" in ref.columns and "dominant_driver" in q.columns:
        rows.append({"scenario":name,"metric":"dominant_driver","spearman_vs_primary":np.nan,
                     "categorical_agreement":float((ref["dominant_driver"].astype(str)==q["dominant_driver"].astype(str)).mean())})
    return rows

def host_scenario(sc,sdir,community_dir,args,out):
    portal=Path(community_dir)/"portal_metabolites_for_hepatic_model.json"
    for cond in [args.baseline_condition,args.treatment_condition]:
        cdir=sdir/f"condition_{cond}"
        expected=cdir/f"{cond}_flux_comparison.csv"
        cmd=[sys.executable,args.hepatic_runner,
             "--source_script",args.hepatic_source_script,
             "--hepatic_model",args.hepatic_model,
             "--expression_data",args.hepatic_expression,
             "--portal_metabolites",str(portal),
             "--condition",cond,
             "--results_dir",str(cdir),
             "--objective_mode",sc["objective"],
             "--pfba_fraction",str(args.pfba_fraction),
             "--solver",args.solver,
             "--diet_bounds",args.diet_bounds,
             "--portal_scaling",str(sc["scale"]),
             "--portal_mode",sc["portal_mode"],
             "--gene_mapping",args.gene_mapping]
        run(cmd,out/"logs"/f"{sc['name']}_{cond}.log",
            resume=args.resume,expected=expected,dry=args.dry_run)

    adir=sdir/"attribution"
    expected=adir/"flux_attribution_analysis.csv"
    cmd=[sys.executable,args.attribution_script,
         "--hepatic_results_dir",str(sdir),
         "--baseline_condition",args.baseline_condition,
         "--treatment_condition",args.treatment_condition,
         "--results_dir",str(adir)]
    run(cmd,out/"logs"/f"{sc['name']}_attribution.log",
        resume=args.resume,expected=expected,dry=args.dry_run)
    return expected

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--community_script",required=True)
    ap.add_argument("--hepatic_runner",required=True)
    ap.add_argument("--hepatic_source_script",required=True)
    ap.add_argument("--attribution_script",required=True)
    ap.add_argument("--metatranscriptome",required=True)
    ap.add_argument("--agora_dir",required=True)
    ap.add_argument("--species_map_primary",required=True)
    ap.add_argument("--species_map_sensitivity",required=True)
    ap.add_argument("--hepatic_model",required=True)
    ap.add_argument("--hepatic_expression",required=True)
    ap.add_argument("--gene_mapping",required=True)
    ap.add_argument("--diet_bounds",required=True)
    ap.add_argument("--expression_cols",default="ND_SCD,DD_HFD")
    ap.add_argument("--baseline_condition",default="ND_SCD")
    ap.add_argument("--treatment_condition",default="DD_HFD")
    ap.add_argument("--tradeoffs",default="0.3,0.5,0.7")
    ap.add_argument("--portal_scales",default="0.05,0.10,0.20")
    ap.add_argument("--pfba_fraction",type=float,default=1.0)
    ap.add_argument("--solver",default="gurobi")
    ap.add_argument("--output",default="RQ4_final_full_v2")
    ap.add_argument("--resume",action="store_true")
    ap.add_argument("--dry-run",action="store_true")
    args=ap.parse_args()

    out=Path(args.output); out.mkdir(parents=True,exist_ok=True)
    primary_inc=out/"species_map_primary_included_only.csv"
    sens_inc=out/"species_map_mapping_sensitivity_included_only.csv"
    print(f"[INFO] primary species={included_map(args.species_map_primary,primary_inc,27)}")
    print(f"[INFO] mapping-sensitivity species={included_map(args.species_map_sensitivity,sens_inc,28)}")

    tradeoffs=[float(x) for x in args.tradeoffs.split(",")]
    scales=[float(x) for x in args.portal_scales.split(",")]
    exprcols=[x.strip() for x in args.expression_cols.split(",")]
    if 0.5 not in tradeoffs: tradeoffs.append(0.5)

    community={}
    for t in sorted(set(tradeoffs)):
        cdir=out/f"community_primary_tradeoff_{tag(t)}"
        expected=cdir/"portal_metabolites_for_hepatic_model.json"
        cmd=[sys.executable,args.community_script,
             "--metatranscriptome",args.metatranscriptome,
             "--agora_dir",args.agora_dir,
             "--expression_cols",*exprcols,
             "--medium_json",args.diet_bounds,
             "--hepatic_model",args.hepatic_model,
             "--species_model_map",str(primary_inc),
             "--strict_species_map",
             "--expected_mapped_species","27",
             "--tradeoff",str(t),"--pfba","--solver",args.solver,
             "--results_dir",str(cdir)]
        run(cmd,out/"logs"/f"community_primary_t{tag(t)}.log",
            resume=args.resume,expected=expected,dry=args.dry_run)
        if not args.dry_run:
            validate_community_output(cdir,exprcols,27)
        community[t]=cdir

    scenarios=scenario_list(tradeoffs,scales)
    pd.DataFrame(scenarios).to_csv(out/"scenario_design.csv",index=False)

    if args.dry_run:
        # Also print mapping sensitivity community command.
        cdir=out/"community_mapping_sensitivity_tradeoff_0p5"
        cmd=[sys.executable,args.community_script,
             "--metatranscriptome",args.metatranscriptome,
             "--agora_dir",args.agora_dir,
             "--expression_cols",*exprcols,
             "--medium_json",args.diet_bounds,
             "--hepatic_model",args.hepatic_model,
             "--species_model_map",str(sens_inc),
             "--strict_species_map","--tradeoff","0.5","--pfba","--solver",args.solver,
             "--results_dir",str(cdir)]
        print("\n> "+" ".join(map(str,cmd)))
        print("[DRY RUN] no models executed.")
        return

    rows=[]; tables={}
    for sc in scenarios:
        sdir=out/"host_scenarios"/sc["name"]
        attr=host_scenario(sc,sdir,community[sc["tradeoff"]],args,out)
        sm,tab=summarize(attr)
        tables[sc["name"]]=tab.set_index("reaction_id")
        rows.append({**sc,"species_map":"primary_27",**sm})

    map_comm=out/"community_mapping_sensitivity_tradeoff_0p5"
    expected=map_comm/"portal_metabolites_for_hepatic_model.json"
    cmd=[sys.executable,args.community_script,
         "--metatranscriptome",args.metatranscriptome,
         "--agora_dir",args.agora_dir,
         "--expression_cols",*exprcols,
         "--medium_json",args.diet_bounds,
         "--hepatic_model",args.hepatic_model,
         "--species_model_map",str(sens_inc),
         "--strict_species_map","--expected_mapped_species","28",
         "--tradeoff","0.5","--pfba","--solver",args.solver,
         "--results_dir",str(map_comm)]
    run(cmd,out/"logs"/"community_mapping_sensitivity.log",
        resume=args.resume,expected=expected,dry=False)
    validate_community_output(map_comm,exprcols,28)

    sc={"name":"mapping_sensitivity_28_species","tradeoff":0.5,
        "objective":"biomass","scale":0.1,"portal_mode":"forced"}
    sdir=out/"host_scenarios"/sc["name"]
    attr=host_scenario(sc,sdir,map_comm,args,out)
    sm,tab=summarize(attr)
    tables[sc["name"]]=tab.set_index("reaction_id")
    rows.append({**sc,"species_map":"mapping_sensitivity_28",**sm})

    summary=pd.DataFrame(rows)
    summary.to_csv(out/"rq4_sensitivity_summary.csv",index=False)

    ref=tables["primary"]
    stability=[]
    for name,q in tables.items():
        stability.extend(compare(ref,q,name))
    pd.DataFrame(stability).to_csv(out/"rq4_reaction_stability_vs_primary.csv",index=False)

    primary=summary[summary["name"]=="primary"].iloc[0]
    ready=[]
    for _,r in summary.iterrows():
        ready.append({
            "scenario":r["name"],
            "species_map":r["species_map"],
            "n_measurable_microbiome":r["n_measurable_microbiome"],
            "delta_n_measurable_vs_primary":r["n_measurable_microbiome"]-primary["n_measurable_microbiome"],
            "n_microbiome_dominant":r["n_microbiome_dominant"],
            "delta_n_microbiome_dominant_vs_primary":r["n_microbiome_dominant"]-primary["n_microbiome_dominant"],
            "mean_ve_microbiome":r["mean_ve_microbiome"],
        })
    pd.DataFrame(ready).to_csv(out/"rq4_interpretation_summary.csv",index=False)

    portals=[]
    for t,cdir in community.items():
        pfile=Path(cdir)/"portal_metabolite_production.csv"
        p=pd.read_csv(pfile) if pfile.exists() else pd.DataFrame()
        portals.append({"scenario":f"primary_tradeoff_{t}","tradeoff":t,
                        "species_map":"primary_27","n_portal_rows":len(p)})
    pfile=map_comm/"portal_metabolite_production.csv"
    p=pd.read_csv(pfile) if pfile.exists() else pd.DataFrame()
    portals.append({"scenario":"mapping_sensitivity_28_species","tradeoff":0.5,
                    "species_map":"mapping_sensitivity_28","n_portal_rows":len(p)})
    pd.DataFrame(portals).to_csv(out/"rq4_portal_and_mapping_summary.csv",index=False)

    print("\n[DONE] RQ4 full final v2")
    print(summary.to_string(index=False))

if __name__=="__main__":
    main()

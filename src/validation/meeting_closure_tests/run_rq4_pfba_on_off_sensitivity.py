#!/usr/bin/env python3
"""
RQ4 pFBA sensitivity at the frozen primary settings.

Compares the frozen primary (community pFBA ON + hepatic pFBA ON) against:
 A. community pFBA OFF + hepatic pFBA ON
 B. community pFBA ON  + hepatic pFBA OFF
 C. community pFBA OFF + hepatic pFBA OFF

This distinguishes the MICOM flux-selection effect from the hepatic pFBA effect.
"""
from __future__ import annotations
import argparse, subprocess, sys, json
from pathlib import Path
import pandas as pd, numpy as np
from scipy.stats import spearmanr

def run(c,cwd):
    print("\n> "+" ".join(map(str,c)))
    p=subprocess.run([str(x) for x in c],cwd=str(cwd))
    if p.returncode: raise RuntimeError(f"failed ({p.returncode})")

def threshold_set(path):
    d=pd.read_csv(path)
    if d["microbiome_attributable"].dtype != bool:
        mask=d["microbiome_attributable"].astype(str).str.lower().isin(["true","1","yes"])
    else:
        mask=d["microbiome_attributable"]
    return set(d.loc[mask,"reaction_id"])

def compare(ref_attr, ref_fc, alt_attr, alt_fc):
    ref=pd.read_csv(ref_attr).set_index("reaction_id")
    alt=pd.read_csv(alt_attr).set_index("reaction_id")
    idx=ref.index.intersection(alt.index)
    x=pd.to_numeric(ref.loc[idx,"delta_microbiome"],errors="coerce")
    y=pd.to_numeric(alt.loc[idx,"delta_microbiome"],errors="coerce")
    ok=x.notna()&y.notna()
    rho=float(spearmanr(x[ok],y[ok]).correlation) if ok.sum()>=3 else np.nan
    rset=threshold_set(ref_fc); aset=threshold_set(alt_fc)
    return {
        "threshold_n":len(aset),
        "overlap_with_primary":len(rset&aset),
        "jaccard_with_primary":len(rset&aset)/len(rset|aset) if rset|aset else np.nan,
        "primary_set_retained_pct":100*len(rset&aset)/len(rset) if rset else np.nan,
        "delta_microbiome_spearman_vs_primary":rho,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project_root",default=".")
    ap.add_argument("--rq4_source_dir",required=True)
    ap.add_argument("--frozen_rq4",default="reviewer_results/RQ4_final_full_v3")
    ap.add_argument("--metatranscriptome",default="data/Meta_GSE104913.csv")
    ap.add_argument("--agora_dir",required=True)
    ap.add_argument("--hepatic_model",default="iMM1415.json")
    ap.add_argument("--hepatic_expression",default="data/GSE182668/male-C57BL6J-GSE182668_HFD_SCD_gene_expression.csv")
    ap.add_argument("--gene_mapping",default="mouse_entrez_to_symbol.csv")
    ap.add_argument("--diet_bounds",default="expanded_diet_bounds_flat.json")
    ap.add_argument("--solver",default="gurobi")
    ap.add_argument("--output",default="reviewer_results/meeting_closure/RQ4_pfba_factorial_sensitivity")
    a=ap.parse_args()

    root=Path(a.project_root).resolve(); src=Path(a.rq4_source_dir).resolve()
    frozen=(root/a.frozen_rq4).resolve(); out=(root/a.output).resolve()
    out.mkdir(parents=True,exist_ok=True)
    sp=frozen/"species_map_primary_included_only.csv"
    if not sp.exists(): raise FileNotFoundError(sp)

    # Frozen community-pFBA portal.
    portal_on=frozen/"community_primary_tradeoff_0p5/portal_metabolites_for_hepatic_model.json"
    if not portal_on.exists(): raise FileNotFoundError(portal_on)

    # Build community-FBA (pFBA off) portal once.
    comm_off=out/"community_pfba_off"
    portal_off=comm_off/"portal_metabolites_for_hepatic_model.json"
    if not portal_off.exists():
        run([sys.executable,src/"rq4_microbiome_community_modeling_final_v3.py",
             "--metatranscriptome",(root/a.metatranscriptome).resolve(),
             "--agora_dir",Path(a.agora_dir).resolve(),
             "--expression_cols","ND_SCD","DD_HFD",
             "--medium_json",(root/a.diet_bounds).resolve(),
             "--hepatic_model",(root/a.hepatic_model).resolve(),
             "--species_model_map",sp,"--strict_species_map","--expected_mapped_species","27",
             "--tradeoff","0.5","--no-pfba","--solver",a.solver,"--results_dir",comm_off],root)

    scenarios=[
        ("communityFBA_hostpFBA",portal_off,True),
        ("communitypFBA_hostFBA",portal_on,False),
        ("communityFBA_hostFBA",portal_off,False),
    ]

    ref_attr=frozen/"host_scenarios/primary/attribution/flux_attribution_analysis.csv"
    ref_fc=frozen/"host_scenarios/primary/condition_DD_HFD/DD_HFD_flux_comparison.csv"
    rows=[{
        "scenario":"primary_communitypFBA_hostpFBA",
        "community_pfba":True,"host_pfba":True,
        "threshold_n":len(threshold_set(ref_fc)),
        "overlap_with_primary":len(threshold_set(ref_fc)),
        "jaccard_with_primary":1.0,
        "primary_set_retained_pct":100.0,
        "delta_microbiome_spearman_vs_primary":1.0,
    }]

    for name,portal,use_host_pfba in scenarios:
        hostroot=out/name
        for cond in ["ND_SCD","DD_HFD"]:
            cdir=hostroot/f"condition_{cond}"
            cmd=[sys.executable,src/"rq4_hepatic_integration_final_v2.py",
                 "--source_script",src/"rq4_hepatic_integration_CORRECTED_v14_source.py",
                 "--hepatic_model",(root/a.hepatic_model).resolve(),
                 "--expression_data",(root/a.hepatic_expression).resolve(),
                 "--portal_metabolites",portal,"--condition",cond,"--results_dir",cdir,
                 "--objective_mode","biomass","--pfba_fraction","1.0","--solver",a.solver,
                 "--diet_bounds",(root/a.diet_bounds).resolve(),
                 "--gene_mapping",(root/a.gene_mapping).resolve(),
                 "--portal_scaling","0.1","--portal_mode","forced"]
            if not use_host_pfba:
                cmd.append("--no-pfba")
            run(cmd,root)

        attr=hostroot/"attribution"
        run([sys.executable,src/"rq4_attribution_analysis_source.py",
             "--hepatic_results_dir",hostroot,"--baseline_condition","ND_SCD",
             "--treatment_condition","DD_HFD","--results_dir",attr],root)

        vals=compare(
            ref_attr,ref_fc,
            attr/"flux_attribution_analysis.csv",
            hostroot/"condition_DD_HFD/DD_HFD_flux_comparison.csv"
        )
        rows.append({
            "scenario":name,
            "community_pfba": "communitypFBA" in name,
            "host_pfba": use_host_pfba,
            **vals,
        })

    res=pd.DataFrame(rows)
    res.to_csv(out/"RQ4_pFBA_factorial_summary.csv",index=False)
    (out/"RQ4_pFBA_factorial_summary.json").write_text(
        json.dumps(res.to_dict(orient="records"),indent=2),encoding="utf-8")
    print("\n"+res.to_string(index=False))

if __name__=="__main__":
    main()

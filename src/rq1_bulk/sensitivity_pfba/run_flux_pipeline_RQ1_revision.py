#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RQ1 revision rerun orchestrator.

Runs the modeling ONCE:
  1) replicate-level three-layer E-Flux + pFBA -> feasible raw fluxes
  2) dataset-effect adjustment -> visualization/QC matrix
  3) dataset-aware comprehensive inference on raw fluxes, using corrected matrix only for PCA
  4) optional interactive descriptive visualization on corrected matrix

This avoids the old duplicated "uncorrected branch / corrected branch" re-modeling.
"""
from __future__ import annotations
import argparse, csv, datetime as dt, hashlib, json, re, shlex, subprocess, sys
from pathlib import Path

try:
    import yaml
except Exception:
    yaml=None

def load_config(path):
    p=Path(path)
    if p.suffix.lower()==".json":
        return json.loads(p.read_text(encoding="utf-8"))
    if p.suffix.lower() in {".yaml",".yml"}:
        if yaml is None: raise RuntimeError("PyYAML is required for YAML configs.")
        return yaml.safe_load(p.read_text(encoding="utf-8"))
    raise ValueError("Config must be JSON or YAML.")

def cli(flag,value):
    if value is None or value is False:return []
    if value is True:return [f"--{flag}"]
    return [f"--{flag}",str(value)]

def run(cmd,log,cwd,dry=False):
    text=" ".join(shlex.quote(str(x)) for x in cmd)
    print(text)
    if dry:return
    with open(log,"a",encoding="utf-8") as f:
        f.write("\n$ "+text+"\n")
        r=subprocess.run(cmd,cwd=cwd,stdout=f,stderr=subprocess.STDOUT,text=True)
    if r.returncode!=0:
        try:
            lines=Path(log).read_text(encoding="utf-8",errors="replace").splitlines()
            tail="\n".join(lines[-80:])
            print("\n--- child-process log tail ---")
            print(tail)
            print("--- end log tail ---\n")
        except Exception:
            pass
        raise RuntimeError(f"Command failed ({r.returncode}): {text}")

def sha256(path,block=1024*1024):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        while True:
            b=f.read(block)
            if not b:break
            h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser(description="Run revised RQ1 pFBA pipeline.")
    ap.add_argument("--config",default="run_flux_config_RQ1_pfba.json")
    ap.add_argument("--dry-run",action="store_true")
    ap.add_argument("--skip_batch",action="store_true")
    ap.add_argument("--skip_interactive",action="store_true")
    a=ap.parse_args()
    cfg=load_config(a.config)
    root=Path(cfg.get("project_root",".")).resolve()
    out=root/Path(cfg["results_dir"])
    out.mkdir(parents=True,exist_ok=True)
    logs=out/"logs"; logs.mkdir(exist_ok=True)
    stamp=dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    log=logs/f"pipeline_{stamp}.log"

    py=str(cfg.get("python_executable",sys.executable))
    scripts=cfg["scripts"]
    required=[
        cfg["input_expression"],cfg["model_file"],cfg["diet_bounds_json"],cfg["mapping_file"],
        scripts["flux"],scripts["batch"],scripts["comprehensive"],scripts["combined"],scripts["fva_summary"],
    ]
    for rel in required:
        p=root/rel
        if not p.exists():raise FileNotFoundError(p)

    # Cheap preflight: validate the expected RQ1 sample counts from the CSV header
    # before launching any genome-scale optimizations. This prevents a silent subset
    # of cohorts from being modeled because of column-detection regressions.
    expected=cfg.get("expected_group_counts",{})
    if expected:
        with open(root/cfg["input_expression"],"r",encoding="utf-8-sig",newline="") as fh:
            header=next(csv.reader(fh))
        observed={}
        for g,nexp in expected.items():
            rx=re.compile(rf"^{re.escape(str(g))}(?:_|$)",re.I)
            observed[g]=sum(bool(rx.search(str(c))) for c in header)
            if observed[g] != int(nexp):
                raise ValueError(
                    f"RQ1 input header count mismatch for {g}: expected {nexp}, observed {observed[g]}. "
                    "Check the merged expression file before running."
                )
        print(f"[PRECHECK] RQ1 sample counts verified: {observed} (total={sum(observed.values())})")

    model_out=out/"modeling"
    common={
        "model_file":cfg["model_file"],
        "results_dir":str(model_out),
        "test_type":cfg.get("test_type","t-test"),
        "diet_bounds_json":cfg["diet_bounds_json"],
        "eflux_quantile":cfg.get("eflux_quantile",0.95),
        "eflux_floor":cfg.get("eflux_floor",0.1),
        "eflux_cap":cfg.get("eflux_cap",1000),
        "objective_id":cfg.get("objective_id","BIOMASS_mm_1_no_glygln"),
        "objective_sense":cfg.get("objective_sense","max"),
        "transporter_strategy":cfg.get("transporter_strategy","either"),
        "edge_abs_diff_threshold":cfg.get("edge_abs_diff_threshold",0.2),
        "mapping_file":cfg["mapping_file"],
        "solve_mode":cfg.get("solve_mode","pfba"),
        "pfba_fraction":cfg.get("pfba_fraction",1.0),
        "solver":cfg.get("solver","gurobi"),
        "fva_mode":cfg.get("fva_mode","representative"),
        "fva_fraction":cfg.get("fva_fraction",1.0),
        "fva_processes":cfg.get("fva_processes",1),
    }
    if cfg.get("fva_targets_file"):common["fva_targets_file"]=cfg["fva_targets_file"]
    cmd=[py,scripts["flux"],cfg["input_expression"]]
    for k,v in common.items():cmd+=cli(k,v)
    cmd+=["--write_replicates_long"]
    if cfg.get("aggregate",False):
        raise ValueError("Primary RQ1 revision must preserve biological replicates; aggregate=true is not allowed.")
    run(cmd,log,root,a.dry_run)

    raw=model_out/"flux_analysis"/"reaction_flux_comparison_extended.csv"
    corrected=out/"visualization"/"reaction_flux_comparison_extended_batch_adjusted.csv"
    if not a.dry_run and not raw.exists():raise FileNotFoundError(raw)

    if not a.skip_batch:
        corrected.parent.mkdir(parents=True,exist_ok=True)
        cmd=[py,scripts["batch"],"--input",str(raw),"--output_csv",str(corrected),
             "--outdir",str(out/"batch_correction_qc")]
        run(cmd,log,root,a.dry_run)
        pca_input=corrected
    else:
        pca_input=raw

    comp_out=out/"comprehensive_analysis"
    cmd=[py,scripts["comprehensive"],"--input",str(raw),"--pca_input",str(pca_input),
         "--output",str(comp_out),"--permutations",str(cfg.get("permutations",999)),
         "--fdr",str(cfg.get("fdr",0.05)),"--effect_size",str(cfg.get("effect_size",0.5))]
    run(cmd,log,root,a.dry_run)

    if cfg.get("fva_mode","none") == "representative":
        fva_dir=model_out/"fva"
        cmd=[py,scripts["fva_summary"],"--fva_dir",str(fva_dir),
             "--baseline","SCD","--targets","HFD,KD,WD",
             "--output",str(out/"fva_robustness")]
        run(cmd,log,root,a.dry_run)

    if not a.skip_interactive:
        cmd=[py,scripts["combined"],"--input",str(pca_input),
             "--output_dir",str(out/"descriptive_interactive")]
        run(cmd,log,root,a.dry_run)

    if not a.dry_run:
        manifest={
            "timestamp":dt.datetime.now().isoformat(),
            "pipeline":"RQ1_revision_pfba",
            "config":cfg,
            "outputs":{"raw_feasible_flux":str(raw),"batch_adjusted_visualization_flux":str(corrected),
                       "comprehensive":str(comp_out)},
            "script_hashes":{s:sha256(root/s) for s in scripts.values()},
            "input_hashes":{
                cfg["input_expression"]:sha256(root/cfg["input_expression"]),
                cfg["model_file"]:sha256(root/cfg["model_file"]),
                cfg["diet_bounds_json"]:sha256(root/cfg["diet_bounds_json"]),
                cfg["mapping_file"]:sha256(root/cfg["mapping_file"]),
            }
        }
        (out/f"manifest_{stamp}.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
        print(f"[OK] RQ1 pipeline complete: {out}")
    return 0

if __name__=="__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
Run all still-open computational action items from the advisor meeting.
Requires the final Complete_Multiscale_Hepatic_Pipeline package plus biological inputs.
"""
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path

def run(c,cwd,dry=False):
    print("\n> "+" ".join(map(str,c)))
    if dry:return
    p=subprocess.run([str(x) for x in c],cwd=str(cwd))
    if p.returncode: raise RuntimeError(f"failed ({p.returncode})")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--project_root",default=".")
    ap.add_argument("--pipeline_root",required=True,help="Complete_Multiscale_Hepatic_Pipeline directory")
    ap.add_argument("--agora_dir",required=True)
    ap.add_argument("--expression",default="GSEMERGED_SCD_HFD_KD_WD_gene_expression.csv")
    ap.add_argument("--model",default="iMM1415.json")
    ap.add_argument("--diet_bounds",default="expanded_diet_bounds_flat.json")
    ap.add_argument("--mapping",default="mouse_entrez_to_symbol.csv")
    ap.add_argument("--positive_controls",default=None)
    ap.add_argument("--frozen_rq1",default="RERUN_RQ1_PFBA")
    ap.add_argument("--frozen_rq4",default="reviewer_results/RQ4_final_full_v3")
    ap.add_argument("--rq4_metatranscriptome",default="data/Meta_GSE104913.csv")
    ap.add_argument("--rq4_hepatic_expression",default="data/GSE182668/male-C57BL6J-GSE182668_HFD_SCD_gene_expression.csv")
    ap.add_argument("--solver",default="gurobi")
    ap.add_argument("--skip_rq4_pfba_toggle",action="store_true")
    ap.add_argument("--dry-run",action="store_true")
    a=ap.parse_args()

    root=Path(a.project_root).resolve(); pipe=Path(a.pipeline_root).resolve()
    here=Path(__file__).resolve().parent
    res=root/"reviewer_results/meeting_closure"; res.mkdir(parents=True,exist_ok=True)
    pc=Path(a.positive_controls).resolve() if a.positive_controls else pipe/"configs/positive_controls.json"

    # 1. Direct three-layer FBA vs pFBA, same GSE101657 samples/settings.
    run([sys.executable,here/"run_direct_fba_vs_pfba.py",
         "--project_root",root,
         "--modeling_script",pipe/"src/rq1_bulk/sensitivity_pfba/map_fixv5_multigroupsv8_layered_manuscript_rerun.py",
         "--expression",a.expression,"--model",a.model,"--diet_bounds",a.diet_bounds,"--mapping",a.mapping,
         "--objective","BIOMASS_mm_1_no_glygln","--solver",a.solver,
         "--output","reviewer_results/meeting_closure/FBA_pFBA_GSE101657"],root,a.dry_run)

    # 2. Cap sensitivity missing from earlier P90/P95/P99/floor sweep.
    capdir=res/"cap_sensitivity_GSE101657"
    run([sys.executable,pipe/"src/benchmarking/sweep_threelayer_pipeline_v2.py",
         "--flux_script",pipe/"src/rq1_bulk/sensitivity_pfba/map_fixv5_multigroupsv8_layered_manuscript_rerun.py",
         "--model",(root/a.model).resolve(),"--rnaseq",(root/a.expression).resolve(),
         "--dataset","GSE101657","--conditions","SCD,HFD,KD","--baseline","SCD",
         "--mapping",(root/a.mapping).resolve(),"--diet_bounds",(root/a.diet_bounds).resolve(),
         "--positive_controls",pc,"--quantiles","0.95","--floors","0.1",
         "--caps","100,1000,10000","--objectives","BIOMASS_mm_1_no_glygln,ATPM",
         "--solver",a.solver,"--results_dir",capdir],root,a.dry_run)
    if not a.dry_run:
        run([sys.executable,here/"summarize_cap_sensitivity.py","--dir",capdir],root,False)

    # 3. ATPM FVA/pFVA of the original RQ1 reaction claims.
    run([sys.executable,pipe/"src/rq1_bulk/sensitivity_pfba/run_rq1_cohort_representative_fva.py",
         "--expression",(root/a.expression).resolve(),"--model",(root/a.model).resolve(),
         "--diet_bounds",(root/a.diet_bounds).resolve(),"--mapping",(root/a.mapping).resolve(),
         "--modeling_script",pipe/"src/rq1_bulk/sensitivity_pfba/map_fixv5_multigroupsv8_layered_manuscript_rerun.py",
         "--rq1_results",(root/a.frozen_rq1).resolve(),
         "--output",res/"RQ1_ATPM_FVA",
         "--objective_id","ATPM","--solver",a.solver,"--fva_fraction","1.0",
         "--parsimony_eps","0.01,0.05","--tolerance","1e-7","--processes","1"],root,a.dry_run)

    # 4. Stronger-than-required RQ4 pFBA-on/off sensitivity.
    if not a.skip_rq4_pfba_toggle:
        run([sys.executable,here/"run_rq4_pfba_on_off_sensitivity.py",
             "--project_root",root,"--rq4_source_dir",pipe/"src/rq4_microbiome",
             "--frozen_rq4",a.frozen_rq4,"--metatranscriptome",a.rq4_metatranscriptome,
             "--agora_dir",a.agora_dir,"--hepatic_model",a.model,
             "--hepatic_expression",a.rq4_hepatic_expression,
             "--gene_mapping",a.mapping,"--diet_bounds",a.diet_bounds,
             "--solver",a.solver],root,a.dry_run)

    print("\n[DONE] Remaining meeting-action tests completed.")
if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Verify final frozen headline counts after a full pipeline run."""
from pathlib import Path
import argparse, json, sys

EXPECTED={
    "RQ1_HFD_anchor_n":20,
    "RQ1_shared_all_three_SCD_referenced_n":27,
    "RQ2_universal_magnitude_core_n":16,
    "RQ3_P99_biomass_WDChow_responsive_unique_reactions":562,
    "RQ4_primary_HFD_threshold_attributable_n":47,
}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--cross-summary",default="reviewer_results/final_cross_rq/cross_RQ_summary_final.json")
    a=ap.parse_args()
    p=Path(a.cross_summary)
    if not p.exists():
        raise FileNotFoundError(p)
    d=json.loads(p.read_text())
    bad=[]
    for k,v in EXPECTED.items():
        obs=d.get(k)
        print(f"{k}: observed={obs}, expected={v}")
        if obs!=v: bad.append((k,obs,v))
    if bad:
        print("[FAIL] frozen headline counts differ.")
        raise SystemExit(2)
    print("[PASS] final frozen headline counts reproduced.")

if __name__=="__main__":
    main()

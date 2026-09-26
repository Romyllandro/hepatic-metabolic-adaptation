#!/usr/bin/env python3
from __future__ import annotations
import argparse, re
from pathlib import Path
import pandas as pd

BFRAG_QUERY = "Bacteroides fragilis NCTC 9343 strain ATCC 25285"
BFRAG_MODEL = "Bacteroides_fragilis_NCTC_9343.mat"
CAMPY_QUERY = "Campylobacter jejuni subsp jejuni NCTC 11168 = ATCC 700819"
CAMPY_MODEL = "Campylobacter_jejuni_subsp_jejuni_NCTC_11168_BN148.mat"

def norm(x):
    return re.sub(r"[^a-z0-9]+","_",str(x).lower()).strip("_")

def locate_row(df, query):
    q=norm(query)
    hits=df[df["species"].astype(str).map(norm)==q]
    if len(hits)!=1:
        raise ValueError(f"Expected exactly one row for {query!r}; found {len(hits)}")
    return hits.index[0]

def validate_models(df, agora_dir):
    missing=[]
    for _,r in df[df["included"].astype(bool)].iterrows():
        p=Path(agora_dir)/str(r["model_file"])
        if not p.exists():
            missing.append(str(p))
    if missing:
        raise FileNotFoundError("Included AGORA model(s) not found:\n  " + "\n  ".join(missing))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--reviewed_map",required=True)
    ap.add_argument("--agora_dir",required=True)
    ap.add_argument("--output_dir",default="rq4_final_maps")
    ap.add_argument("--expected_exact",type=int,default=26)
    ap.add_argument("--expected_primary",type=int,default=27)
    ap.add_argument("--expected_sensitivity",type=int,default=28)
    args=ap.parse_args()

    out=Path(args.output_dir); out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(args.reviewed_map)
    required={"species","model_file","included","mapping_status","mapping_basis"}
    miss=required-set(df.columns)
    if miss:
        raise ValueError(f"Reviewed map missing columns: {sorted(miss)}")

    df["included"]=df["included"].astype(str).str.lower().isin(["true","1","yes"])
    exact=df[df["included"] & df["mapping_basis"].astype(str).str.lower().eq("exact")].copy()
    if len(exact)!=args.expected_exact:
        raise ValueError(
            f"Expected {args.expected_exact} exact mappings, found {len(exact)}."
        )

    primary=df.copy()
    primary["included"]=False
    primary["model_file"]=""
    primary["model_path"]=""
    primary["mapping_status"]="excluded_from_primary"
    primary["manual_review"]="complete"

    for idx,r in exact.iterrows():
        primary.loc[idx,"included"]=True
        primary.loc[idx,"model_file"]=r["model_file"]
        primary.loc[idx,"model_path"]=str((Path(args.agora_dir)/str(r["model_file"])).resolve())
        primary.loc[idx,"mapping_status"]="validated_exact"
        primary.loc[idx,"mapping_basis"]="exact"

    bidx=locate_row(primary,BFRAG_QUERY)
    bpath=Path(args.agora_dir)/BFRAG_MODEL
    if not bpath.exists():
        raise FileNotFoundError(f"Verified B. fragilis synonym model not found: {bpath}")
    primary.loc[bidx,"included"]=True
    primary.loc[bidx,"model_file"]=BFRAG_MODEL
    primary.loc[bidx,"model_path"]=str(bpath.resolve())
    primary.loc[bidx,"mapping_status"]="validated_strain_synonym"
    primary.loc[bidx,"mapping_basis"]="verified_strain_synonym_NCTC9343_ATCC25285"
    primary.loc[bidx,"manual_review"]="complete"
    primary.loc[bidx,"notes"]="Primary RQ4: verified NCTC 9343 / ATCC 25285 strain synonym."

    cidx=locate_row(primary,CAMPY_QUERY)
    primary.loc[cidx,"included"]=False
    primary.loc[cidx,"model_file"]=""
    primary.loc[cidx,"model_path"]=""
    primary.loc[cidx,"mapping_status"]="excluded_strain_variant_proxy"
    primary.loc[cidx,"mapping_basis"]="NCTC11168_BN148_variant_not_primary"
    primary.loc[cidx,"manual_review"]="complete"
    primary.loc[cidx,"notes"]="Primary RQ4 excludes BN148 because it is a documented NCTC 11168 strain variant, not an exact identity match."

    if int(primary["included"].sum())!=args.expected_primary:
        raise ValueError(f"Primary map count={int(primary['included'].sum())}, expected {args.expected_primary}.")
    validate_models(primary,args.agora_dir)

    sensitivity=primary.copy()
    cpath=Path(args.agora_dir)/CAMPY_MODEL
    if not cpath.exists():
        raise FileNotFoundError(f"Campylobacter sensitivity model not found: {cpath}")
    sensitivity.loc[cidx,"included"]=True
    sensitivity.loc[cidx,"model_file"]=CAMPY_MODEL
    sensitivity.loc[cidx,"model_path"]=str(cpath.resolve())
    sensitivity.loc[cidx,"mapping_status"]="included_strain_variant_proxy_sensitivity"
    sensitivity.loc[cidx,"mapping_basis"]="strain_variant_proxy_NCTC11168_BN148"
    sensitivity.loc[cidx,"manual_review"]="complete"
    sensitivity.loc[cidx,"notes"]="Mapping sensitivity only: BN148 is a documented NCTC 11168 variant; not treated as an exact primary mapping."

    if int(sensitivity["included"].sum())!=args.expected_sensitivity:
        raise ValueError(f"Sensitivity map count={int(sensitivity['included'].sum())}, expected {args.expected_sensitivity}.")
    validate_models(sensitivity,args.agora_dir)

    p1=out/"rq4_species_map_primary_27.csv"
    p2=out/"rq4_species_map_mapping_sensitivity_28.csv"
    primary.to_csv(p1,index=False)
    sensitivity.to_csv(p2,index=False)

    pd.DataFrame([
        {"species":BFRAG_QUERY,"primary_decision":"include","sensitivity_decision":"include",
         "model_file":BFRAG_MODEL,"basis":"verified strain synonym"},
        {"species":CAMPY_QUERY,"primary_decision":"exclude","sensitivity_decision":"include",
         "model_file":CAMPY_MODEL,"basis":"strain-variant proxy; sensitivity only"}
    ]).to_csv(out/"rq4_species_map_decisions.csv",index=False)

    qc=[
        f"Input reviewed map: {Path(args.reviewed_map).resolve()}",
        f"Rows/species in reviewed map: {len(df)}",
        f"Validated exact mappings recovered: {len(exact)}",
        f"Primary strict map included: {int(primary['included'].sum())}",
        f"Mapping-sensitivity map included: {int(sensitivity['included'].sum())}",
        "",
        f"Primary synonym: {BFRAG_QUERY} -> {BFRAG_MODEL}",
        f"Sensitivity-only proxy: {CAMPY_QUERY} -> {CAMPY_MODEL}",
        "",
        "All included AGORA model files exist."
    ]
    (out/"rq4_species_map_final_QC.txt").write_text("\n".join(qc),encoding="utf-8")
    print("\n".join(qc))

if __name__=="__main__":
    main()

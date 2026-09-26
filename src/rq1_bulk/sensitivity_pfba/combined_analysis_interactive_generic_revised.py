#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Descriptive cross-dataset visualization for a batch-adjusted flux matrix.

No reaction-level hypothesis tests are performed here. This script is intentionally
restricted to PCA, centroid distances, hierarchical clustering, and optional Plotly
HTML. Primary inference belongs in comprehensive_flux_analysis_revised.py.
"""
from __future__ import annotations
import argparse, re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import linkage, dendrogram

try:
    import plotly.express as px
    _PLOTLY=True
except Exception:
    _PLOTLY=False

ALIASES={
 "SCD":{"SCD","SC","CD","CHOW","ND","CONTROL","LFD","NCD","CTRL"},
 "HFD":{"HFD","HF","HIGHFAT"},
 "KD":{"KD","KETO","KETOGENIC"},
 "WD":{"WD","WESTERN","WESTERNDIET"},
}
def group_of(s):
    toks=re.split(r"[^A-Z0-9]+",str(s).upper())
    for t in toks:
        for c,a in ALIASES.items():
            if t==c or t in a:return c
    return str(s).split("_",1)[0].upper()
def dataset_of(s):
    m=re.search(r"GSE\d+",str(s),re.I)
    if m:return m.group(0).upper()
    m=re.search(r"GSM\d+",str(s),re.I)
    return m.group(0).upper() if m else "UNKNOWN"
def samples(df,suffix):
    bad=re.compile(r"(MeanFlux|StdFlux|SEMFlux|_N)$",re.I)
    return [str(c) for c in df.columns if str(c).endswith(suffix) and not bad.search(str(c))]

def main():
    p=argparse.ArgumentParser(description="Descriptive PCA/clustering of batch-adjusted flux profiles.")
    p.add_argument("--input","-i",required=True)
    p.add_argument("--output_dir","-o",default="descriptive_batch_adjusted")
    p.add_argument("--sample_suffix",default="_Flux")
    p.add_argument("--no_interactive",action="store_true")
    a=p.parse_args()
    out=Path(a.output_dir); out.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(a.input); sc=samples(df,a.sample_suffix)
    if not sc: raise ValueError("No sample columns detected.")
    meta=pd.DataFrame({"SampleID":sc,"Group":[group_of(c) for c in sc],"Dataset":[dataset_of(c) for c in sc]})
    if (meta.Dataset=="UNKNOWN").any(): raise ValueError("Unparsed dataset IDs in sample names.")
    X=df[sc].apply(pd.to_numeric,errors="coerce").fillna(0).to_numpy(float).T
    Xz=StandardScaler().fit_transform(X)
    ncomp=min(10,Xz.shape[0]-1,Xz.shape[1])
    pc=PCA(n_components=ncomp,random_state=42); Z=pc.fit_transform(Xz)
    for i in range(min(5,ncomp)):meta[f"PC{i+1}"]=Z[:,i]
    meta.to_csv(out/"pca_scores_metadata.csv",index=False)

    groups=sorted(meta.Group.unique())
    full_cent={g:Xz[(meta.Group==g).to_numpy()].mean(0) for g in groups}
    pc_cent={g:Z[(meta.Group==g).to_numpy(),:2].mean(0) for g in groups}
    Df=np.array([[np.linalg.norm(full_cent[a]-full_cent[b]) for b in groups] for a in groups])
    Dp=np.array([[np.linalg.norm(pc_cent[a]-pc_cent[b]) for b in groups] for a in groups])
    pd.DataFrame(Df,index=groups,columns=groups).to_csv(out/"centroid_distances_full_zscore.csv")
    pd.DataFrame(Dp,index=groups,columns=groups).to_csv(out/"centroid_distances_PC12.csv")

    fig,ax=plt.subplots(figsize=(7,6))
    for g in groups:
        m=(meta.Group==g).to_numpy()
        ax.scatter(Z[m,0],Z[m,1],s=50,alpha=.8,label=g)
    ax.set_xlabel(f"PC1 ({pc.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pc.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("Batch-adjusted flux PCA (descriptive)")
    ax.legend(); fig.tight_layout(); fig.savefig(out/"PCA_by_diet.png",dpi=300); plt.close(fig)

    L=linkage(Xz,method="average",metric="euclidean")
    fig,ax=plt.subplots(figsize=(13,6))
    dendrogram(L,labels=meta.SampleID.tolist(),leaf_rotation=90,leaf_font_size=6,ax=ax)
    ax.set_title("Hierarchical clustering of batch-adjusted flux profiles")
    fig.tight_layout(); fig.savefig(out/"hierarchical_clustering.png",dpi=300); plt.close(fig)

    if _PLOTLY and not a.no_interactive:
        idf=meta.copy()
        fig=px.scatter(idf,x="PC1",y="PC2",color="Group",symbol="Dataset",hover_name="SampleID",
                       title="Batch-adjusted flux PCA")
        fig.write_html(out/"PCA_interactive.html",include_plotlyjs=True)
    print(f"[OK] Descriptive outputs: {out}")
    return 0
if __name__=="__main__":
    raise SystemExit(main())

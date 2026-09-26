#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Revised downstream analysis for RQ1.

Primary inference:
    original feasible sample-level pFBA fluxes
    + diet effect estimated while adjusting for dataset.

Visualization:
    optional batch-adjusted flux matrix for PCA/clustering/centroid displays.

This separation prevents statistically batch-adjusted vectors (which need not
satisfy S*v=0) from becoming the primary reaction-level inferential data.
"""

from __future__ import annotations
import argparse, json, math, re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.stats import hypergeom
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm


RANDOM_SEED = 42

GROUP_ALIASES = {
    "SCD": {"SCD","SC","CD","CHOW","ND","CONTROL","STDCHOW","STANDARDCHOW","LFD","NCD","CTRL","CONTROL"},
    "HFD": {"HFD","HF","HIGHFAT","HIGH_FAT","HIGH-FAT"},
    "KD": {"KD","KETO","KETOGENIC","KETO_DIET","KETO-DIET"},
    "WD": {"WD","WESTERN","WESTERNDIET"},
}


def canonical_group(text):
    up = str(text).upper()
    toks = re.split(r"[^A-Z0-9]+", up)
    for tok in toks:
        for canon, aliases in GROUP_ALIASES.items():
            if tok == canon or tok in aliases:
                return canon
    return str(text).split("_",1)[0].upper()


def parse_dataset(text):
    # Prefer study-level GSE identifiers over sample-level GSM identifiers.
    m = re.search(r"GSE\d+", str(text), flags=re.I)
    if m:
        return m.group(0).upper()
    m = re.search(r"GSM\d+", str(text), flags=re.I)
    return m.group(0).upper() if m else "UNKNOWN"


def detect_sample_columns(df, suffix="_Flux"):
    ignore = re.compile(r"(MeanFlux|StdFlux|SEMFlux|_N)$", re.I)
    return [str(c) for c in df.columns if str(c).endswith(suffix) and not ignore.search(str(c))]


def load_flux_table(path, suffix="_Flux"):
    df = pd.read_csv(path)
    required = ["ReactionID","ReactionName","Subsystem"]
    miss = [c for c in required if c not in df.columns]
    if miss:
        raise ValueError(f"Missing required columns: {miss}")
    sample_cols = detect_sample_columns(df, suffix)
    if not sample_cols:
        raise ValueError("No sample-level *_Flux columns detected.")
    meta = pd.DataFrame({
        "SampleID": sample_cols,
        "Group": [canonical_group(c) for c in sample_cols],
        "Dataset": [parse_dataset(c) for c in sample_cols],
    })
    if (meta["Dataset"]=="UNKNOWN").any():
        bad=meta.loc[meta["Dataset"]=="UNKNOWN","SampleID"].head().tolist()
        raise ValueError(f"Could not parse GSE/GSM dataset IDs from sample columns: {bad}")
    X = df[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float).T
    return df, sample_cols, meta, X


def cohen_d(x, y):
    x=np.asarray(x,float); y=np.asarray(y,float)
    x=x[np.isfinite(x)]; y=y[np.isfinite(y)]
    if len(x)<2 or len(y)<2:
        return np.nan
    diff=float(np.mean(x)-np.mean(y))
    dof=len(x)+len(y)-2
    pv=((len(x)-1)*np.var(x,ddof=1)+(len(y)-1)*np.var(y,ddof=1))/max(dof,1)
    if pv<=0:
        return 0.0 if abs(diff)<1e-15 else float(np.sign(diff)*np.inf)
    return diff/np.sqrt(pv)


def comparison_order(groups, control="SCD"):
    groups=sorted(set(groups))
    out=[]
    if control in groups:
        out += [(g,control) for g in groups if g!=control]
    non=[g for g in groups if g!=control]
    out += list(combinations(non,2))
    return out


def shared_datasets(meta, g1, g2):
    ds=[]
    for d, sub in meta.groupby("Dataset"):
        gs=set(sub["Group"])
        if g1 in gs and g2 in gs:
            ds.append(d)
    return sorted(ds)


def bh(p):
    p=np.asarray(p,float)
    q=np.full_like(p,np.nan)
    ok=np.isfinite(p)
    if ok.any():
        q[ok]=multipletests(p[ok], method="fdr_bh")[1]
    return q


def zscore_features(X):
    X=np.asarray(X,float)
    mu=np.nanmean(X,axis=0)
    sd=np.nanstd(X,axis=0,ddof=1)
    sd[~np.isfinite(sd) | (sd==0)] = 1.0
    return (np.nan_to_num(X,nan=mu) - mu)/sd


def pseudo_f_from_distance(D, labels):
    labels=np.asarray(labels)
    N=len(labels); groups=np.unique(labels); k=len(groups)
    if N<=k or k<2:
        return np.nan, np.nan
    A=D**2
    SST=A.sum()/N
    SSW=0.0
    for g in groups:
        idx=np.where(labels==g)[0]; ng=len(idx)
        if ng>1:
            SSW += A[np.ix_(idx,idx)].sum()/ng
    SSB=max(0.0,SST-SSW)
    F=(SSB/(k-1))/(SSW/(N-k)) if SSW>0 else np.inf
    R2=SSB/SST if SST>0 else 0.0
    return float(F),float(R2)


def restricted_permanova(X, labels, strata, n_perm=999, seed=42):
    X=zscore_features(X)
    D=squareform(pdist(X,metric="euclidean"))
    labels=np.asarray(labels,dtype=object)
    strata=np.asarray(strata,dtype=object)
    F,R2=pseudo_f_from_distance(D,labels)
    rng=np.random.default_rng(seed)
    count=0
    uniq_strata=np.unique(strata)
    for _ in range(n_perm):
        perm=labels.copy()
        for s in uniq_strata:
            idx=np.where(strata==s)[0]
            perm[idx]=rng.permutation(perm[idx])
        Fp,_=pseudo_f_from_distance(D,perm)
        if np.isfinite(Fp) and Fp>=F:
            count+=1
    p=(count+1)/(n_perm+1)
    return {"F":F,"R2":R2,"p":float(p),"N":len(labels),"n_strata":len(uniq_strata)}


def pca_analysis(X, meta, outdir, prefix="batch_adjusted"):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    Xz=zscore_features(X)
    ncomp=min(20,Xz.shape[0]-1,Xz.shape[1])
    pca=PCA(n_components=ncomp,random_state=RANDOM_SEED)
    Z=pca.fit_transform(Xz)

    # Full standardized-flux-space centroid distances (recommended for Fig 2b).
    groups=sorted(meta["Group"].unique())
    cent_full={g:Xz[(meta["Group"]==g).to_numpy()].mean(axis=0) for g in groups}
    cent_pc={g:Z[(meta["Group"]==g).to_numpy(),:2].mean(axis=0) for g in groups}
    full=np.zeros((len(groups),len(groups)))
    pc12=np.zeros_like(full)
    for i,g1 in enumerate(groups):
        for j,g2 in enumerate(groups):
            full[i,j]=np.linalg.norm(cent_full[g1]-cent_full[g2])
            pc12[i,j]=np.linalg.norm(cent_pc[g1]-cent_pc[g2])
    pd.DataFrame(full,index=groups,columns=groups).to_csv(outdir/f"{prefix}_centroid_distances_full_zscore.csv")
    pd.DataFrame(pc12,index=groups,columns=groups).to_csv(outdir/f"{prefix}_centroid_distances_PC12.csv")

    scores=pd.DataFrame({"SampleID":meta["SampleID"],"Group":meta["Group"],"Dataset":meta["Dataset"]})
    for i in range(min(5,Z.shape[1])):
        scores[f"PC{i+1}"]=Z[:,i]
    scores.to_csv(outdir/f"{prefix}_pca_scores.csv",index=False)

    fig,ax=plt.subplots(figsize=(7,6))
    for g in groups:
        m=(meta["Group"]==g).to_numpy()
        ax.scatter(Z[m,0],Z[m,1],s=50,alpha=.8,label=g)
        c=cent_pc[g]
        ax.scatter([c[0]],[c[1]],marker="D",s=90,edgecolor="black")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("PCA of batch-adjusted flux profiles" if "batch" in prefix else "PCA of flux profiles")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir/f"{prefix}_PCA.png",dpi=300)
    plt.close(fig)
    return {
        "pc1":float(pca.explained_variance_ratio_[0]),
        "pc2":float(pca.explained_variance_ratio_[1]),
        "centroid_full_file":str(outdir/f"{prefix}_centroid_distances_full_zscore.csv"),
    }


def build_design(submeta, g1):
    treat=(submeta["Group"].values==g1).astype(float)
    ds=pd.get_dummies(pd.Categorical(submeta["Dataset"]),drop_first=True,dtype=float).to_numpy()
    X=np.column_stack([np.ones(len(submeta)),treat,ds])
    return X


def reaction_statistics(raw_df, meta, fdr=0.05, effect=0.5):
    results={}
    per_dataset={}
    comps=comparison_order(meta["Group"].unique())
    for g1,g2 in comps:
        dsets=shared_datasets(meta,g1,g2)
        cname=f"{g1}_vs_{g2}"
        if not dsets:
            print(f"[SKIP] {cname}: no dataset contains both groups.")
            continue
        mask=meta["Dataset"].isin(dsets) & meta["Group"].isin([g1,g2])
        smeta=meta.loc[mask].reset_index(drop=True)
        cols=smeta["SampleID"].tolist()
        Y=raw_df[cols].apply(pd.to_numeric,errors="coerce").to_numpy(dtype=float)  # rxn x samples
        design=build_design(smeta,g1)
        if np.linalg.matrix_rank(design)<design.shape[1]:
            raise ValueError(f"Rank-deficient inferential design for {cname}.")
        rows=[]
        for i,row in raw_df.iterrows():
            y=Y[i,:]
            good=np.isfinite(y)
            if good.sum()<4:
                continue
            Xg=design[good,:]; yg=y[good]
            if np.linalg.matrix_rank(Xg)<Xg.shape[1]:
                continue
            fit=sm.OLS(yg,Xg).fit(cov_type="HC3")
            beta=float(fit.params[1]); p=float(fit.pvalues[1])
            v1=yg[smeta.loc[good,"Group"].values==g1]
            v2=yg[smeta.loc[good,"Group"].values==g2]
            rows.append({
                "ReactionID":row["ReactionID"],"ReactionName":row["ReactionName"],
                "Subsystem":row["Subsystem"],"Contrast":cname,
                "Treatment":g1,"Reference":g2,
                "SharedDatasets":";".join(dsets),"N_datasets":len(dsets),
                "N_treatment":len(v1),"N_reference":len(v2),
                f"{g1}_mean":float(np.mean(v1)),f"{g2}_mean":float(np.mean(v2)),
                "AdjustedMeanDiff":beta,
                "Cohen_d":cohen_d(v1,v2),
                "p_value":p,
            })
        sdf=pd.DataFrame(rows)
        if len(sdf):
            sdf["q_value"]=bh(sdf["p_value"].values)
            sdf["Significant"]=(sdf["q_value"]<=fdr) & (sdf["Cohen_d"].abs()>=effect)
            sdf["Direction"]=np.where(sdf["AdjustedMeanDiff"]>0,"Up","Down")
            sdf=sdf.sort_values(["q_value","p_value","ReactionID"])
        results[cname]=sdf

        # Dataset-specific effects for rank-product and consistency.
        dsrows=[]
        for d in dsets:
            md=meta[(meta["Dataset"]==d)&(meta["Group"].isin([g1,g2]))]
            c1=md.loc[md["Group"]==g1,"SampleID"].tolist()
            c2=md.loc[md["Group"]==g2,"SampleID"].tolist()
            for _,r in raw_df.iterrows():
                a=pd.to_numeric(r[c1],errors="coerce").to_numpy(float)
                b=pd.to_numeric(r[c2],errors="coerce").to_numpy(float)
                a=a[np.isfinite(a)]; b=b[np.isfinite(b)]
                if len(a)<1 or len(b)<1:
                    continue
                if len(a)>=2 and len(b)>=2:
                    with np.errstate(all="ignore"):
                        _,pv=stats.ttest_ind(a,b,equal_var=False,nan_policy="omit")
                else:
                    pv=np.nan
                dsrows.append({
                    "ReactionID":r["ReactionID"],"Dataset":d,"Contrast":cname,
                    "MeanDiff":float(np.mean(a)-np.mean(b)),
                    "AbsMeanDiff":float(abs(np.mean(a)-np.mean(b))),
                    "Cohen_d":cohen_d(a,b),"Welch_p":float(pv) if np.isfinite(pv) else np.nan,
                    "N_treatment":len(a),"N_reference":len(b)
                })
        per_dataset[cname]=pd.DataFrame(dsrows)
        print(f"[STATS] {cname}: datasets={len(dsets)}, reactions={len(sdf)}, significant={int(sdf['Significant'].sum()) if len(sdf) else 0}")
    return results,per_dataset


def rank_product(per_dataset):
    out={}
    for cname,df in per_dataset.items():
        if df.empty or df["Dataset"].nunique()<2:
            continue
        wide=df.pivot(index="ReactionID",columns="Dataset",values="AbsMeanDiff")
        ranks=wide.rank(axis=0,ascending=False,method="average")
        norm=ranks.divide(ranks.count(axis=0),axis=1)
        rp=np.exp(np.log(norm.clip(lower=1e-12)).mean(axis=1))
        dwide=df.pivot(index="ReactionID",columns="Dataset",values="MeanDiff")
        direction_agreement=np.maximum((dwide>0).mean(axis=1),(dwide<0).mean(axis=1))
        med_abs=wide.median(axis=1)
        rdf=pd.DataFrame({
            "ReactionID":rp.index,"RankProduct":rp.values,
            "DirectionAgreement":direction_agreement.reindex(rp.index).values,
            "MedianAbsMeanDiff":med_abs.reindex(rp.index).values,
            "N_datasets":wide.notna().sum(axis=1).reindex(rp.index).values,
        }).sort_values("RankProduct")
        rdf["RankProductRank"]=np.arange(1,len(rdf)+1)
        out[cname]=rdf
    return out


def subsystem_and_enrichment(stats_by_contrast, fdr=0.05):
    subsys_out={}; enrich_out={}
    for cname,sdf in stats_by_contrast.items():
        if sdf.empty:
            continue
        gr=[]
        for ss,g in sdf.groupby("Subsystem",dropna=False):
            sig=g[g["Significant"]]
            gr.append({
                "Subsystem":ss,"N_reactions":len(g),"N_significant":len(sig),
                "N_up_significant":int((sig["AdjustedMeanDiff"]>0).sum()),
                "N_down_significant":int((sig["AdjustedMeanDiff"]<0).sum()),
                "MeanAdjustedDiff":float(g["AdjustedMeanDiff"].mean()),
                "MeanAbsCohen_d":float(g["Cohen_d"].replace([np.inf,-np.inf],np.nan).abs().mean()),
            })
        subsys_out[cname]=pd.DataFrame(gr).sort_values("N_significant",ascending=False)

        sig=sdf[sdf["Significant"]]
        M=len(sdf); N=len(sig); er=[]
        if N:
            for ss,g in sdf.groupby("Subsystem",dropna=False):
                n=len(g); k=int(sig["Subsystem"].eq(ss).sum())
                if k==0: continue
                p=float(hypergeom.sf(k-1,M,n,N))
                er.append({"Subsystem":ss,"k_significant":k,"n_total":n,
                           "expected":n*N/M,"enrichment_ratio":k/(n*N/M),"p_value":p})
        edf=pd.DataFrame(er)
        if len(edf):
            edf["q_value"]=bh(edf["p_value"].values)
            edf["Significant"]=edf["q_value"]<=fdr
            edf=edf.sort_values(["q_value","enrichment_ratio"],ascending=[True,False])
        enrich_out[cname]=edf
    return subsys_out,enrich_out


def write_cytoscape(stats_by_contrast, outdir):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    for cname,sdf in stats_by_contrast.items():
        sig=sdf[sdf["Significant"]]
        if sig.empty: continue
        edges=pd.DataFrame({
            "source":sig["ReactionID"],"target":sig["Subsystem"],
            "interaction":"in_subsystem","adjusted_mean_diff":sig["AdjustedMeanDiff"],
            "cohen_d":sig["Cohen_d"],"q_value":sig["q_value"],
        })
        edges.to_csv(outdir/f"edges_{cname}.csv",index=False)
        nodes_rxn=pd.DataFrame({
            "id":sig["ReactionID"],"label":sig["ReactionName"],"type":"reaction",
            "subsystem":sig["Subsystem"],"effect":sig["AdjustedMeanDiff"],"q_value":sig["q_value"]
        })
        ss=sorted(sig["Subsystem"].dropna().astype(str).unique())
        nodes_ss=pd.DataFrame({"id":ss,"label":ss,"type":"subsystem"})
        pd.concat([nodes_rxn,nodes_ss],ignore_index=True).to_csv(outdir/f"nodes_{cname}.csv",index=False)


def plot_volcano(stats_by_contrast,outdir):
    outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
    for cname,sdf in stats_by_contrast.items():
        if sdf.empty: continue
        q=np.clip(sdf["q_value"].to_numpy(float),1e-300,1)
        x=sdf["AdjustedMeanDiff"].to_numpy(float)
        sig=sdf["Significant"].to_numpy(bool)
        fig,ax=plt.subplots(figsize=(7,6))
        ax.scatter(x[~sig],-np.log10(q[~sig]),s=10,alpha=.45)
        if sig.any():
            ax.scatter(x[sig],-np.log10(q[sig]),s=16,alpha=.8)
        ax.set_xlabel("Dataset-adjusted mean flux difference")
        ax.set_ylabel("-log10(BH-FDR)")
        ax.set_title(cname.replace("_vs_"," vs "))
        fig.tight_layout()
        fig.savefig(outdir/f"volcano_{cname}.png",dpi=300)
        plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description="Dataset-aware RQ1 downstream analysis on feasible pFBA fluxes.")
    ap.add_argument("--input","-i",required=True,help="Original feasible pFBA wide flux CSV.")
    ap.add_argument("--pca_input",default=None,help="Optional batch-adjusted CSV used only for descriptive PCA/centroids.")
    ap.add_argument("--output","-o",default="RQ1_comprehensive_analysis")
    ap.add_argument("--sample_suffix",default="_Flux")
    ap.add_argument("--permutations",type=int,default=999)
    ap.add_argument("--fdr",type=float,default=0.05)
    ap.add_argument("--effect_size",type=float,default=0.5)
    args=ap.parse_args()

    out=Path(args.output)
    csvdir=out/"csv_outputs"; figdir=out/"figures"; netdir=out/"cytoscape_networks"; repdir=out/"reports"
    for d in [csvdir,figdir,netdir,repdir]: d.mkdir(parents=True,exist_ok=True)

    raw,sample_cols,meta,Xraw=load_flux_table(args.input,args.sample_suffix)
    meta.to_csv(csvdir/"sample_metadata.csv",index=False)
    print(pd.crosstab(meta["Dataset"],meta["Group"]))

    # Descriptive PCA/centroid analysis on batch-adjusted representation when supplied.
    pca_source="raw"
    Xpca=Xraw
    if args.pca_input:
        corr,corr_cols,cmeta,Xcorr=load_flux_table(args.pca_input,args.sample_suffix)
        if set(corr_cols)!=set(sample_cols):
            raise ValueError("pca_input sample columns do not match raw input.")
        corr_idx={c:i for i,c in enumerate(corr_cols)}
        Xpca=np.vstack([Xcorr[corr_idx[c],:] for c in sample_cols])
        pca_source="batch_adjusted"
    pca_res=pca_analysis(Xpca,meta,csvdir,prefix=pca_source)

    # Restricted-permutation PERMANOVA on original feasible fluxes (primary)
    perma=[]
    overall=restricted_permanova(Xraw,meta["Group"].values,meta["Dataset"].values,args.permutations,RANDOM_SEED)
    overall.update({"Contrast":"OVERALL","Data":"raw_feasible"})
    perma.append(overall)
    pair=[]
    for g1,g2 in comparison_order(meta["Group"].unique()):
        dsets=shared_datasets(meta,g1,g2)
        if not dsets: continue
        m=meta["Dataset"].isin(dsets)&meta["Group"].isin([g1,g2])
        rr=restricted_permanova(Xraw[m.to_numpy(),:],meta.loc[m,"Group"].values,meta.loc[m,"Dataset"].values,
                                args.permutations,RANDOM_SEED)
        rr.update({"Contrast":f"{g1}_vs_{g2}","Data":"raw_feasible","SharedDatasets":";".join(dsets)})
        pair.append(rr)
    pdf=pd.DataFrame(pair)
    if len(pdf):
        pdf["q_value"]=bh(pdf["p"].values)
    pd.DataFrame(perma).to_csv(csvdir/"permanova_overall_restricted_raw.csv",index=False)
    pdf.to_csv(csvdir/"permanova_pairwise_restricted_raw.csv",index=False)

    # Optional corrected PERMANOVA only as sensitivity/descriptive comparison.
    if args.pca_input:
        pairc=[]
        oc=restricted_permanova(Xpca,meta["Group"].values,meta["Dataset"].values,args.permutations,RANDOM_SEED)
        pd.DataFrame([{**oc,"Contrast":"OVERALL","Data":"batch_adjusted"}]).to_csv(
            csvdir/"permanova_overall_restricted_batch_adjusted.csv",index=False)
        for g1,g2 in comparison_order(meta["Group"].unique()):
            dsets=shared_datasets(meta,g1,g2)
            if not dsets: continue
            m=meta["Dataset"].isin(dsets)&meta["Group"].isin([g1,g2])
            rr=restricted_permanova(Xpca[m.to_numpy(),:],meta.loc[m,"Group"].values,meta.loc[m,"Dataset"].values,
                                    args.permutations,RANDOM_SEED)
            rr.update({"Contrast":f"{g1}_vs_{g2}","SharedDatasets":";".join(dsets)})
            pairc.append(rr)
        cdf=pd.DataFrame(pairc)
        if len(cdf): cdf["q_value"]=bh(cdf["p"].values)
        cdf.to_csv(csvdir/"permanova_pairwise_restricted_batch_adjusted.csv",index=False)

    stats_by, per_ds=reaction_statistics(raw,meta,args.fdr,args.effect_size)
    for cname,sdf in stats_by.items():
        sdf.to_csv(csvdir/f"reaction_stats_{cname}.csv",index=False)
    for cname,ddf in per_ds.items():
        ddf.to_csv(csvdir/f"dataset_specific_effects_{cname}.csv",index=False)

    rp=rank_product(per_ds)
    for cname,rdf in rp.items():
        rdf.to_csv(csvdir/f"rank_product_{cname}.csv",index=False)

    subsys,enrich=subsystem_and_enrichment(stats_by,args.fdr)
    for cname,df in subsys.items(): df.to_csv(csvdir/f"subsystem_analysis_{cname}.csv",index=False)
    for cname,df in enrich.items(): df.to_csv(csvdir/f"pathway_enrichment_{cname}.csv",index=False)

    write_cytoscape(stats_by,netdir)
    plot_volcano(stats_by,figdir)

    summary={
        "input_raw":str(args.input),"pca_input":str(args.pca_input) if args.pca_input else None,
        "n_reactions":len(raw),"n_samples":len(meta),
        "groups":meta["Group"].value_counts().to_dict(),
        "datasets":meta["Dataset"].value_counts().to_dict(),
        "pca_source":pca_source,"pca":pca_res,
        "inference":"Dataset-adjusted OLS with HC3 robust SE on original feasible pFBA fluxes; only datasets containing both contrasted groups are used.",
        "fdr":args.fdr,"effect_size_threshold":args.effect_size,
        "contrasts":{
            c:{
                "n_tested":len(sdf),
                "n_significant":int(sdf["Significant"].sum()) if len(sdf) else 0,
                "shared_datasets":sdf["SharedDatasets"].iloc[0] if len(sdf) else ""
            } for c,sdf in stats_by.items()
        },
        "non_estimable_contrasts":[
            f"{a}_vs_{b}" for a,b in comparison_order(meta["Group"].unique())
            if len(shared_datasets(meta,a,b))==0
        ]
    }
    (repdir/"analysis_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    with open(repdir/"README_RESULTS.txt","w",encoding="utf-8") as f:
        f.write("PRIMARY INFERENCE uses original feasible pFBA fluxes with dataset adjustment.\n")
        f.write("PCA/centroid displays may use the batch-adjusted representation supplied by --pca_input.\n")
        f.write("KD_vs_WD is skipped if no dataset contains both diets.\n")
        f.write("Rank-product is computed only when >=2 datasets contain both groups.\n")
    print(f"[OK] Revised comprehensive analysis complete: {out}")
    return 0

if __name__=="__main__":
    raise SystemExit(main())

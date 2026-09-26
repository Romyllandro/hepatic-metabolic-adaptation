#!/usr/bin/env python3
"""
Regenerate revised manuscript Figure 2 (RQ1) and Figure 3 (RQ2).

Required inputs
---------------
1. RERUN_RQ1_PFBA.zip (or extracted directory)
2. RERUN_RQ2_PFBA_FVA.zip (or extracted directory)
3. targeted_parsimonious_FVA.zip (or extracted directory)

Outputs
-------
Figure2_RQ1_revised.png
Figure3_RQ2_revised.png
plus component panels used to assemble each figure.

The script uses the frozen revised definitions:
- RQ1 primary reaction significance from dataset-aware pFBA inference.
- RQ2 magnitude responsiveness at |mean HFD-SCD| >= 0.20.
- Targeted RQ2 parsimonious-FVA robustness at 1% and 5%.
"""
import argparse, io, json, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from PIL import Image, ImageDraw

class Source:
    def __init__(self,path):
        self.path=Path(path); self.is_zip=self.path.suffix.lower()==".zip"
        self.z=zipfile.ZipFile(self.path) if self.is_zip else None
    def _find(self,suffix):
        if self.is_zip:
            m=[n for n in self.z.namelist() if n.endswith(suffix)]
            if not m: raise FileNotFoundError(suffix)
            return sorted(m,key=len)[0]
        m=list(self.path.rglob(Path(suffix).name))
        if not m: raise FileNotFoundError(suffix)
        return m[0]
    def csv(self,suffix):
        x=self._find(suffix)
        return pd.read_csv(io.BytesIO(self.z.read(x))) if self.is_zip else pd.read_csv(x)
    def json(self,suffix):
        x=self._find(suffix)
        if self.is_zip: return json.loads(self.z.read(x))
        return json.loads(Path(x).read_text())
    def close(self):
        if self.z: self.z.close()

def save(fig,path,dpi=180):
    fig.savefig(path,dpi=dpi,bbox_inches="tight"); plt.close(fig)

def trim(im,pad=6):
    bg=Image.new(im.mode,im.size,"white")
    a=np.asarray(im).astype(int); b=np.asarray(bg).astype(int)
    diff=Image.fromarray(np.max(np.abs(a-b),axis=2).astype(np.uint8))
    box=diff.getbbox()
    if not box:return im
    l,t,r,b=box
    return im.crop((max(0,l-pad),max(0,t-pad),min(im.width,r+pad),min(im.height,b+pad)))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--rq1",required=True)
    ap.add_argument("--rq2",required=True)
    ap.add_argument("--rq2-pfva",required=True)
    ap.add_argument("--output",default="revised_figures")
    a=ap.parse_args()
    out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
    r1=Source(a.rq1); r2=Source(a.rq2); pfs=Source(a.rq2_pfva)

    # ---------------- FIGURE 2 ----------------
    pca=r1.csv("batch_adjusted_pca_scores.csv")
    summary=r1.json("analysis_summary.json")
    dist=r1.csv("batch_adjusted_centroid_distances_full_zscore.csv")
    dist=dist.set_index(dist.columns[0]).loc[["SCD","HFD","KD","WD"],["SCD","HFD","KD","WD"]]
    perma=r1.csv("permanova_overall_restricted_raw.csv").iloc[0]
    contrasts=["HFD_vs_SCD","KD_vs_SCD","WD_vs_SCD"]
    st={c:r1.csv(f"reaction_stats_{c}.csv").set_index("ReactionID") for c in contrasts}
    shared=sorted(set.intersection(*[set(d.index[d["Significant"]==True]) for d in st.values()]))

    # a
    fig=plt.figure(figsize=(5,4)); ax=fig.add_axes([.14,.15,.8,.75])
    for g in ["SCD","HFD","KD","WD"]:
        d=pca[pca.Group==g]; ax.scatter(d.PC1,d.PC2,s=24,alpha=.75,label=g)
    ax.axhline(0,lw=.6,alpha=.35); ax.axvline(0,lw=.6,alpha=.35)
    ax.set_xlabel(f"PC1 ({summary['pca']['pc1']*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({summary['pca']['pc2']*100:.1f}% variance)")
    ax.set_title("a | Batch-adjusted global flux space",loc="left",fontweight="bold")
    ax.legend(frameon=False,fontsize=8,ncol=2)
    ax.text(.02,.02,f"Restricted PERMANOVA on raw feasible pFBA\n$R^2$={perma['R2']:.3f}, p={perma['p']:.3f}",
            transform=ax.transAxes,fontsize=7,va="bottom")
    p2a=out/"Fig2a_PCA.png"; save(fig,p2a)

    # b1
    fig=plt.figure(figsize=(4,3.6)); ax=fig.add_axes([.18,.18,.70,.68])
    arr=dist.values.astype(float); im=ax.imshow(arr)
    ax.set_xticks(range(4),dist.columns); ax.set_yticks(range(4),dist.index)
    ax.set_title("b | Batch-adjusted centroid distance",loc="left",fontweight="bold",fontsize=10)
    for i in range(4):
        for j in range(4): ax.text(j,i,f"{arr[i,j]:.1f}",ha="center",va="center",fontsize=8)
    fig.colorbar(im,ax=ax,fraction=.05,pad=.04,label="Euclidean distance")
    p2b1=out/"Fig2b1_centroid.png"; save(fig,p2b1)

    # b2
    order=["SCD","HFD","KD","WD"]; count=np.full((4,4),np.nan)
    pairs={("SCD","HFD"):61,("SCD","KD"):50,("SCD","WD"):199,("HFD","KD"):2,("HFD","WD"):286}
    for (x,y),v in pairs.items():
        i,j=order.index(x),order.index(y); count[i,j]=count[j,i]=v
    fig=plt.figure(figsize=(4,3.6)); ax=fig.add_axes([.18,.18,.70,.68]); im=ax.imshow(np.ma.masked_invalid(count))
    ax.set_xticks(range(4),order); ax.set_yticks(range(4),order)
    ax.set_title("Significant reactions (estimable contrasts)",loc="left",fontweight="bold",fontsize=10)
    for i in range(4):
        for j in range(4):
            if i==j: txt="—"
            elif {order[i],order[j]}=={"KD","WD"}: txt="NE"
            else: txt=f"{int(count[i,j])}" if np.isfinite(count[i,j]) else "NE"
            ax.text(j,i,txt,ha="center",va="center",fontsize=8)
    fig.colorbar(im,ax=ax,fraction=.05,pad=.04,label="FDR < 0.05")
    p2b2=out/"Fig2b2_counts.png"; save(fig,p2b2)

    # c
    core=pd.DataFrame([(rid,st["HFD_vs_SCD"].loc[rid,"Subsystem"]) for rid in shared],columns=["ReactionID","Subsystem"])
    core["is_rgroup"]=core.Subsystem.eq("R Group Synthesis")
    core=core.sort_values(["is_rgroup","Subsystem","ReactionID"]).reset_index(drop=True)
    fig=plt.figure(figsize=(10,5.4)); ax=fig.add_axes([.18,.10,.78,.82]); y=np.arange(len(core))
    marks={"HFD_vs_SCD":"o","KD_vs_SCD":"s","WD_vs_SCD":"^"}; labs={"HFD_vs_SCD":"HFD","KD_vs_SCD":"KD","WD_vs_SCD":"WD"}
    offs={"HFD_vs_SCD":-.18,"KD_vs_SCD":0,"WD_vs_SCD":.18}
    for c in contrasts:
        vals=[st[c].loc[r,"AdjustedMeanDiff"] for r in core.ReactionID]
        ax.scatter(vals,y+offs[c],s=18,marker=marks[c],label=labs[c])
    ax.axvline(0,lw=.8); ax.set_yticks(y,[r+"†" if q else r for r,q in zip(core.ReactionID,core.is_rgroup)],fontsize=6.5)
    ax.invert_yaxis(); ax.set_xlabel("Dataset-adjusted mean flux difference vs SCD (model flux units)")
    ax.set_title("c | Shared pFBA-associated reaction signature (27 reactions; all directionally concordant)",loc="left",fontweight="bold")
    ax.legend(frameon=False,ncol=3,fontsize=8,loc="lower right")
    ax.text(.99,.01,"† iMM1415 R-group bookkeeping/pseudo-reaction",transform=ax.transAxes,ha="right",fontsize=6.5)
    p2c=out/"Fig2c_shared27.png"; save(fig,p2c)

    # d
    selected=["Transport, Extracellular","Extracellular exchange","Transport, Mitochondrial","Transport, Peroxisomal",
              "Glycolysis/Gluconeogenesis","Citric Acid Cycle","Folate Metabolism","Fatty acid activation",
              "Glycerophospholipid Metabolism","Nucleotides","Oxidative Phosphorylation","Arginine and Proline Metabolism"]
    ds=[]
    for c,label in zip(contrasts,["HFD","KD","WD"]):
        q=r1.csv(f"subsystem_analysis_{c}.csv").set_index("Subsystem")
        fig=plt.figure(figsize=(3.4,4.3)); ax=fig.add_axes([.42,.12,.54,.78]); yy=np.arange(len(selected))
        up=[q.loc[s,"N_up_significant"] if s in q.index else 0 for s in selected]
        down=[q.loc[s,"N_down_significant"] if s in q.index else 0 for s in selected]
        ax.barh(yy,[-v for v in down],height=.7,label="Down"); ax.barh(yy,up,height=.7,label="Up"); ax.axvline(0,lw=.7)
        ax.set_yticks(yy,selected if label=="HFD" else [""]*len(selected),fontsize=6); ax.invert_yaxis()
        ax.set_title(label,fontweight="bold"); ax.set_xlabel("Significant reactions",fontsize=7)
        if label=="WD": ax.legend(frameon=False,fontsize=7,loc="lower right")
        p=out/f"Fig2d_{label}.png"; save(fig,p); ds.append(p)

    W,H=1428,1252; can=Image.new("RGB",(W,H),"white")
    A=trim(Image.open(p2a).convert("RGB")); B1=trim(Image.open(p2b1).convert("RGB")); B2=trim(Image.open(p2b2).convert("RGB"))
    A.thumbnail((470,340)); B1.thumbnail((380,330)); B2.thumbnail((380,330))
    can.paste(A,(25,10)); can.paste(B1,(515,10)); can.paste(B2,(1010,10))
    C=trim(Image.open(p2c).convert("RGB")); C.thumbnail((1360,440)); can.paste(C,((W-C.width)//2,350))
    ImageDraw.Draw(can).text((30,805),"d | Pathway-level direction of significant dietary flux responses",fill="black")
    x=20
    for p in ds:
        im=trim(Image.open(p).convert("RGB")); im.thumbnail((450,405)); can.paste(im,(x,835)); x+=465
    can.save(out/"Figure2_RQ1_revised.png",dpi=(300,300))

    # ---------------- FIGURE 3 ----------------
    per=r2.csv("RQ2_per_strain_summary.csv")
    cross=r2.csv("RQ2_cross_strain_reaction_summary.csv")
    core16=r2.csv("RQ2_universal_magnitude_core.csv")
    long=r2.csv("RQ2_all_strain_reaction_stats_long.csv")
    sens=r2.csv("RQ2_threshold_sensitivity.csv")
    pf=pfs.csv("RQ2_targeted_pfva_cross_strain_summary.csv")
    disp={"129S1SvImJ":"129S1/SvImJ","C57BL6J":"C57BL/6J","AJ":"A/J","CASTEiJ":"CAST/EiJ",
          "DBA2J":"DBA/2J","PWKPhJ":"PWK/PhJ","WSBEiJ":"WSB/EiJ","NODShiLtJ":"NOD/ShiLtJ","NZOHlLtJ":"NZO/HlLtJ"}
    so=["129S1SvImJ","AJ","C57BL6J","CASTEiJ","DBA2J","NODShiLtJ","NZOHlLtJ","PWKPhJ","WSBEiJ"]

    d=per.set_index("Strain").loc[so].reset_index()
    fig=plt.figure(figsize=(6.2,4.1)); ax=fig.add_axes([.25,.17,.72,.72]); xx=np.arange(len(d)); w=.38
    ax.bar(xx-w/2,d.N_magnitude_responsive,width=w,label="Magnitude-responsive |Δv|≥0.20")
    sp=d.N_statistically_supported.astype(float).values; ni=so.index("NZOHlLtJ"); sp[ni]=np.nan
    ax.bar(xx+w/2,sp,width=w,label="FDR<0.05 & |d|≥0.5"); ax.text(ni+w/2,1,"NE",ha="center",fontsize=7)
    ax.set_xticks(xx,[disp[s] for s in so],rotation=45,ha="right",fontsize=7); ax.set_ylabel("Number of reactions")
    ax.set_title("a | Strain-specific HFD response",loc="left",fontweight="bold"); ax.legend(frameon=False,fontsize=7)
    p3a=out/"Fig3a_strains.png"; save(fig,p3a)

    hist=cross.loc[cross.N_magnitude_responsive>0,"N_magnitude_responsive"].value_counts().sort_index()
    fig=plt.figure(figsize=(4.7,3)); ax=fig.add_axes([.15,.20,.80,.68]); x=np.arange(1,10); vals=[int(hist.get(i,0)) for i in x]
    ax.bar(x,vals)
    for xi,v in zip(x,vals): ax.text(xi,v+.8,str(v),ha="center",fontsize=7)
    ax.set_xticks(x); ax.set_xlabel("Number of strains with |Δv|≥0.20"); ax.set_ylabel("Reactions")
    ax.set_title("b | Conservation spectrum (105-reaction union)",loc="left",fontweight="bold",fontsize=9)
    p3b1=out/"Fig3b1_spectrum.png"; save(fig,p3b1)

    fig=plt.figure(figsize=(4.7,2.7)); ax=fig.add_axes([.16,.22,.79,.65])
    ax.plot(sens.theta,sens.N_union,marker="o",label="Union (≥1 strain)")
    ax.plot(sens.theta,sens.N_universal,marker="o",label="Universal (9/9)"); ax.axvline(.20,lw=.8)
    ax.set_xlabel("Magnitude threshold θ"); ax.set_ylabel("Number of reactions")
    ax.set_title("Threshold sensitivity",loc="left",fontweight="bold",fontsize=9); ax.legend(frameon=False,fontsize=7)
    p3b2=out/"Fig3b2_threshold.png"; save(fig,p3b2)

    mat=[]
    for rid in core16.ReactionID:
        mat.append([float(long[(long.ReactionID==rid)&(long.Strain==s)].HFD_minus_SCD.iloc[0]) for s in so])
    mat=np.asarray(mat); m=np.nanmax(np.abs(mat))
    fig=plt.figure(figsize=(9.7,5.2)); ax=fig.add_axes([.18,.13,.68,.79])
    im=ax.imshow(mat,aspect="auto",norm=TwoSlopeNorm(vmin=-m,vcenter=0,vmax=m))
    ax.set_xticks(np.arange(9),[disp[s] for s in so],rotation=45,ha="right",fontsize=7)
    ax.set_yticks(np.arange(len(core16)),core16.ReactionID,fontsize=7)
    ax.set_title("c | Universal magnitude core: Δ flux (HFD−SCD) and targeted parsimonious-FVA robustness",loc="left",fontweight="bold",fontsize=10)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]): ax.text(j,i,f"{mat[i,j]:+.2f}",ha="center",va="center",fontsize=5.3)
    q1=pf[pf.ParsimonyEpsilon==.01].set_index("ReactionID"); q5=pf[pf.ParsimonyEpsilon==.05].set_index("ReactionID")
    ax.set_xlim(-.5,11.8); ax.text(9.65,-1.25,"pFVA\n1%",ha="center",fontweight="bold",fontsize=7); ax.text(10.8,-1.25,"pFVA\n5%",ha="center",fontweight="bold",fontsize=7)
    for i,rid in enumerate(core16.ReactionID):
        ax.text(9.65,i,f"{int(q1.loc[rid,'N_match_primary'])}/9",ha="center",va="center",fontsize=6.5)
        ax.text(10.8,i,f"{int(q5.loc[rid,'N_match_primary'])}/9",ha="center",va="center",fontsize=6.5)
    fig.colorbar(im,ax=ax,fraction=.035,pad=.03,label="HFD−SCD flux")
    p3c=out/"Fig3c_core.png"; save(fig,p3c)

    W,H=1429,1011; can=Image.new("RGB",(W,H),"white")
    A=trim(Image.open(p3a).convert("RGB")); A.thumbnail((720,370))
    B1=trim(Image.open(p3b1).convert("RGB")); B2=trim(Image.open(p3b2).convert("RGB")); B1.thumbnail((650,220)); B2.thumbnail((650,190))
    can.paste(A,(15,15)); can.paste(B1,(760,10)); can.paste(B2,(760,225))
    C=trim(Image.open(p3c).convert("RGB")); C.thumbnail((1385,560)); can.paste(C,((W-C.width)//2,430))
    can.save(out/"Figure3_RQ2_revised.png",dpi=(300,300))
    r1.close(); r2.close(); pfs.close()
    print(out/"Figure2_RQ1_revised.png")
    print(out/"Figure3_RQ2_revised.png")

if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Capture the exact local environment used for a reproducibility release."""
from pathlib import Path
import argparse, json, platform, subprocess, sys, datetime

def run(args):
    try:
        return subprocess.check_output(args,text=True,stderr=subprocess.STDOUT)
    except Exception as e:
        return f"ERROR: {e}\n"

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",default="environment_capture")
    a=ap.parse_args()
    out=Path(a.output); out.mkdir(parents=True,exist_ok=True)
    meta={
        "timestamp":datetime.datetime.now().isoformat(),
        "python":sys.version,
        "executable":sys.executable,
        "platform":platform.platform(),
        "machine":platform.machine(),
    }
    (out/"platform.json").write_text(json.dumps(meta,indent=2),encoding="utf-8")
    (out/"pip_freeze.txt").write_text(run([sys.executable,"-m","pip","freeze"]),encoding="utf-8")
    (out/"conda_environment.yml").write_text(run(["conda","env","export"]),encoding="utf-8")
    (out/"gurobi_version.txt").write_text(run([sys.executable,"-c","import gurobipy as g; print(g.gurobi.version())"]),encoding="utf-8")
    (out/"cobra_version.txt").write_text(run([sys.executable,"-c","import cobra; print(cobra.__version__)"]),encoding="utf-8")
    (out/"micom_version.txt").write_text(run([sys.executable,"-c","import micom; print(micom.__version__)"]),encoding="utf-8")
    print(f"[DONE] wrote exact environment snapshot to {out.resolve()}")

if __name__=="__main__":
    main()

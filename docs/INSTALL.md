# Installation

## 1. Python environment

The production analyses were run on Windows 10/11 with Python 3.9 in a conda environment named
`metabolic-modeling`. Any recent Python (3.9–3.11) should work.

```bash
conda env create -f environment.yml
conda activate hepatic-flux
```

or, in an existing environment:

```bash
pip install -r requirements.txt
```

| Package | Used for |
|---|---|
| cobra (COBRApy) | model loading, FBA/pFBA/FVA |
| gurobipy | LP solver (licence required, see below) |
| micom | microbial community modeling (microbiome stage only) |
| riptide | RIPTiDe comparator (E. coli benchmark only) |
| pandas, numpy, scipy, statsmodels, scikit-learn | data handling and statistics |
| matplotlib, seaborn, pillow, openpyxl, pyyaml | figures, Excel tables, manifests |

## 2. Gurobi

Every optimisation stage uses Gurobi through COBRApy/optlang. Free academic licences are available at
https://www.gurobi.com/academia/. After installing, confirm it works:

```bash
python -c "import gurobipy; m = gurobipy.Model(); print('Gurobi OK', gurobipy.gurobi.version())"
```

### Why the solver matters

FBA problems for genome-scale models often have many optimal flux vectors that share one objective value.
Different solvers, or different Gurobi versions and settings, can return different vertices of that
optimal set. We checked this directly: the A/J strain with GLPK gave the same biomass optimum
(0.04036) as Gurobi, but 64 instead of 96 reactions with |Δv| ≥ 0.20.

* To reproduce the **exact** published numbers, use Gurobi. Its deterministic settings are in
  `src/common/deterministic_solver_v2.py`: Threads=1, Method=1 (dual simplex), NumericFocus=3, Seed=1.
* Other solvers (GLPK, CPLEX) are supported through `--solver`. Expect differences at the level of
  individual reactions, but not in objective values.

This solution dependence is the reason the manuscript reports pFBA and FVA sensitivity analyses.

## 3. Record your environment

Before archiving your own run, capture the exact package versions:

```bash
python src/utils/capture_environment.py --output environment_capture
```

## 4. Check the installation

```bash
cp pipeline_config.example.json pipeline_config.json
python run_pipeline.py --config pipeline_config.json --stages preflight
python run_pipeline.py --config pipeline_config.json --stages verify --use-reference
```

The second command should end with `34/34 checks passed`.

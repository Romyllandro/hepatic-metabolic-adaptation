# Meeting action-item closure tests

Most advisor-meeting action items are already completed in the frozen revision.
This package targets the few remaining empirical gaps so that the response can say
they were **tested**, rather than only justified.

## Still-open tests this package closes

1. **Direct FBA vs pFBA** on the same GSE101657 SCD/HFD/KD samples, same three-layer
   constraints, same biomass-associated objective.
2. **E-Flux cap sensitivity** at cap = 100, 1000, 10000, holding P95/floor 0.1 fixed,
   under biomass and ATPM.
3. **ATPM FVA/pFVA**, because the existing FVA was performed under the primary
   biomass-associated objective and the meeting explicitly asked whether an
   alternative objective remains robust to alternative optima.
4. **RQ4 MICOM pFBA on/off** at the primary settings. This is stronger than the
   original meeting request (which asked for literature support plus an RQ1
   FBA-vs-pFBA test), but it directly tests whether the pFBA choice changes the
   microbiome→liver conclusion.

## Command

Assuming:
- your biological project is `C:\Users\rolan\manuscript_run\pn_pipeline`
- you extracted `Complete_Multiscale_Hepatic_Pipeline_CODE.zip` to
  `C:\Users\rolan\Complete_Multiscale_Hepatic_Pipeline`

copy this folder into the project as `meeting_tests`, then run:

```bat
python meeting_tests\run_all_remaining_meeting_tests.py ^
  --project_root C:\Users\rolan\manuscript_run\pn_pipeline ^
  --pipeline_root C:\Users\rolan\Complete_Multiscale_Hepatic_Pipeline ^
  --agora_dir C:\Users\rolan\AGORA2_mat ^
  --solver gurobi
```

First inspect the commands with:

```bat
python meeting_tests\run_all_remaining_meeting_tests.py ^
  --project_root C:\Users\rolan\manuscript_run\pn_pipeline ^
  --pipeline_root C:\Users\rolan\Complete_Multiscale_Hepatic_Pipeline ^
  --agora_dir C:\Users\rolan\AGORA2_mat ^
  --solver gurobi ^
  --dry-run
```

## Outputs

```text
reviewer_results\meeting_closure\
  FBA_pFBA_GSE101657\
    fba\
    pfba\
    comparison\
      FBA_vs_pFBA_summary.json
      reaction_level_FBA_vs_pFBA.csv

  cap_sensitivity_GSE101657\
    cap_sensitivity_summary.csv
    ...

  RQ1_ATPM_FVA\
    ...

  RQ4_no_pfba_sensitivity\
    RQ4_pFBA_on_off_summary.json
    ...
```

## Interpretation thresholds

Do not define "consistent" as numerically identical. Evaluate:
- Spearman response correlation;
- direction agreement;
- top-reaction Jaccard;
- preservation of pathway/cell-type hierarchy;
- FVA interval direction;
- whether the primary conclusion changes.

Suggested language:
- **highly consistent**: direction >=95% and rho >=0.8, without loss of the core conclusion;
- **moderately consistent / quantitatively sensitive**: direction >=90% but rho or top-set
  overlap noticeably lower;
- **assumption-sensitive**: major rank/set changes that alter reaction-level claims;
- **inconsistent**: central direction/pathway/cell-type conclusion reverses.

Those labels are reporting conventions, not formal statistical cutoffs.


## v2 hotfix note

The first closure package wrote the temporary cohort matrix as:

`GSE101657_SCD_HFD_expression.csv`

but the production RQ1 model intentionally requires a parser-compatible filename
ending in:

`_<GROUPS>_gene_expression.csv`

v2 fixes this to a name such as:

`GSE101657_SCD_HFD_KD_gene_expression.csv`

and evaluates both HFD-vs-SCD and KD-vs-SCD when those conditions are present.

v2 also strengthens the RQ4 pFBA sensitivity into a small factorial test:
- community pFBA off / host pFBA on,
- community pFBA on / host pFBA off,
- both pFBA off,
all compared with the frozen pFBA/pFBA primary result.

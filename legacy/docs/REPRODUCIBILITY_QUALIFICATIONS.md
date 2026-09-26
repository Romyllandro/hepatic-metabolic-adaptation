# Reproducibility qualifications

This package represents the **final corrected analysis path**, not every historical
script used during development.

Notably excluded from the production path:
- early FBA-only variants;
- old per-reaction flux-batch-correction inferential workflow;
- benchmark versions that matched 0 expression genes;
- the invalid RQ4 v2 strict-map handoff that mapped 0/57 species;
- old cross-RQ summaries built before RQ3/RQ4 were frozen.

External resources are not bundled when redistribution/licensing is unclear.
The input manifest documents what must be supplied.

The environment files provide a practical reconstruction but are not an exact
lockfile because the exact local package versions used in the user's workstation
were not available in this packaging runtime. Run `capture_environment.py` on the
production workstation and archive that output for byte-level environment provenance.

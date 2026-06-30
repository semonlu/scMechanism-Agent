# Full Workflow Log

- Result directory: `submission/example/results`
- Report: `submission/example/manuscript_report.md`
- Environment setup/check scripts: `scripts/env_setup/*.ps1`
- Latest integrated environment check: `logs/environment_check_20260630_203000`
- Seurat core workflow: `submission/example/scripts/03_run_gse223751_seurat.R`
- Course modules runner: `submission/example/scripts/04_run_course_modules.ps1`
- Report generator: `scripts/write_analysis_report.py`
- Known-bad older pseudotime artifacts removed: old root cell file, duplicate old sessionInfo file, and PDF inspection scratch images.
- Current Monocle3 report figure uses multi-state `marker_support_label` input rather than a single Fibroblasts-only subset.

# Implementation Plan

The skill uses a layered workflow:

1. Parse the clinical or biological question.
2. Diagnose public accession, supplementary file list, or local file format.
3. Ask for missing sample metadata, group design, organism, tissue, and comparison.
4. Check R/Python environment requirements before Seurat V5 course-derived execution.
5. Select or render reproducible R/Python scripts.
6. Review output folders and tables.
7. Draft cautious Methods, Results, figure legends, limitations, and validation suggestions.
8. Validate the skill with 5 fixed platform cases.

Reusable components:

- `SKILL.md`: routing and high-level workflow.
- `agents/`: role-specific reasoning steps.
- `references/`: format rules, QC rules, pseudotime rules, CellChat limits, environment setup, and platform checklist.
- `scripts/`: deterministic helpers for diagnosis, planning, rendering, validation, report writing, and upload packaging.
- `examples/validation_input_output_comparison.md`: 5-case platform validation report.


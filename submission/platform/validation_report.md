# Five-Case Validation Report

The canonical validation report is:

- `examples/validation_input_output_comparison.md`

Platform reviewers should see five case types:

1. clinical research question.
2. 10x MEX supplementary file diagnosis.
3. h5ad/RDS processed object diagnosis.
4. FASTQ/SRA raw-data boundary handling.
5. uploaded marker/DEG/enrichment/CellChat result interpretation.

The purpose is not to rerun a complete single-cell pipeline five times. The purpose is to prove that the Skill can be invoked on the platform and can stably complete the decision loop:

```text
judge -> plan -> generate code route -> review result -> write report
```


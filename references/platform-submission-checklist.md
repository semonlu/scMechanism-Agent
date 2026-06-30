# Medical AI Skill Platform Submission Checklist

Use this checklist before uploading the skill to the Medical AI Skill development platform.

Source requirement: `D:/单细胞测试/医疗AI Skill开发平台_参赛者操作手册.pdf`.

## Required Folder Shape

The uploadable skill folder must include:

- `SKILL.md`: the top-level skill instructions.
- `agents/*.md`: role files that carry the domain workflow.
- `references/`: supporting rules, checklists, and reference material.
- `examples/`: example prompts and validation evidence.

This skill may also include:

- `scripts/`: deterministic helpers for environment checks, workflow rendering, result validation, report writing, and packaging.
- `templates/`: reusable R/Python/Markdown templates.

Do not include in the platform upload package:

- `.git/`, Git metadata, or local editor files.
- raw expression matrices, FASTQ/SRA files, Seurat RDS objects, h5ad/loom objects, or private clinical metadata.
- local environment check folders such as `seuratv5_environment_check/`.
- generated cache folders such as `__pycache__/`.
- heavy full-run submission outputs unless the platform specifically requests them.

## Required Validation Evidence

Before submission, prepare at least 5 familiar, conclusion-clear validation cases.

For each case, record:

- user input or task prompt.
- expected skill behavior.
- observed skill output.
- pass/fail judgment.
- what file or rule would be changed if the output fails.

The canonical comparison report for this skill is:

- `examples/validation_input_output_comparison.md`

## Required Local Checks

Run:

```powershell
python scripts/validate_platform_skill.py --skill-root .
python scripts/validate_full_workflow.py --project-root . --example-root submission/example --out-md submission/example/full_workflow_validation.md
```

When building an upload package, run:

```powershell
python scripts/package_platform_skill.py --skill-root . --out dist/scMechanism-Agent-skill.zip
```

The package script uses an allowlist so local data and full analysis artifacts are not uploaded accidentally.


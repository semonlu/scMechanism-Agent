# Medical AI Skill Platform Submission Checklist

Use this checklist before uploading the skill to the Medical AI Skill development platform.

Source requirement: `D:/单细胞测试/医疗AI Skill开发平台_参赛者操作手册.pdf`.

## Required Folder Shape

The uploadable skill folder must include:

- `SKILL.md`: the top-level skill instructions.
- `README.md`: a short project-facing overview for reviewers.
- `environment.yml`: lightweight Python/Scanpy reproducibility entrypoint.
- `agents/*.md`: role files that carry the domain workflow.
- `references/`: supporting rules, checklists, and reference material.
- `examples/`: example prompts and validation evidence.
- `submission/platform/`: concise contest-facing materials for demonstration.
- `SECURITY_AND_PRIVACY.md`: data, privacy, and safety boundaries.

This skill may also include:

- `scripts/`: deterministic helpers for environment checks, workflow rendering, result validation, report writing, and packaging.
- `templates/`: reusable R/Python/Markdown templates.

Do not include in the platform upload package:

- `.git/`, Git metadata, or local editor files.
- raw expression matrices, FASTQ/SRA files, Seurat RDS objects, h5ad/loom objects, or private clinical metadata.
- local environment check folders such as `seuratv5_environment_check/`.
- generated cache folders such as `__pycache__/`.
- heavy full-run example outputs unless the platform specifically requests them.

For this repository, package `submission/platform/` for contest demonstration. Keep the heavier local full-run example under `submission/example/` in GitHub/local development, not in the platform upload zip.

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

## Recommended Contest Materials

Prepare these reviewer-facing files:

- `scMechanism-Agent-skill.zip`
- Skill introduction.
- clinical pain point.
- implementation plan.
- business value.
- 5-case input-output comparison report.
- platform demonstration script.
- security and privacy statement.
- GitHub link.

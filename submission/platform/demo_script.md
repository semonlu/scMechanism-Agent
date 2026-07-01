# Platform Demonstration Script

Use this sequence during the contest demonstration.

## 1. Upload and Enable

1. Open the Medical AI Skill platform.
2. Go to Skill Management.
3. Upload `dist/scMechanism-Agent-skill.zip`.
4. Enable the uploaded Skill.

## 2. Show Skill Structure

Show that the uploaded Skill contains:

- `SKILL.md`
- `agents/`
- `references/`
- `templates/`
- `examples/`
- `submission/platform/`
- `SECURITY_AND_PRIVACY.md`

## 3. Run 5 Validation Cases

Run the five prompts in:

- `examples/validation_input_output_comparison.md`

For each case, record:

- input prompt.
- expected output.
- platform observed output summary or screenshot.
- pass/fail judgment.
- triggered agent or rule files when visible in detailed mode.

## 4. Explain the Closed Loop

Use this narrative:

```text
I encoded single-cell analysis experience into a reusable Skill.
The Skill can diagnose input data, plan analysis, generate code, check results, and draft reports.
The 5 validation cases show stable behavior across clinical question parsing, 10x format diagnosis, h5ad/RDS diagnosis, FASTQ/SRA boundary handling, and uploaded result interpretation.
```

## 5. Safety Statement

State that the Skill supports research and hypothesis generation only. It does not make clinical diagnoses or treatment recommendations and does not upload private matrices or patient data.


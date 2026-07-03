# scMechanism-Agent Architecture

## Purpose

scMechanism-Agent is a local-first single-cell research workflow controller. It
helps a user diagnose public or local scRNA-seq inputs, choose defensible
analysis modules, render runnable scripts, validate outputs, and write cautious
mechanism-oriented reports. Heavy computation remains local and auditable.

## Layers

### Skill Control Layer

`SKILL.md` is the routing contract. It decides whether a request is about input
diagnosis, planning, code generation, environment preparation, result review,
or reporting. The `agents/*.md` files provide specialized role instructions for
those modes.

### Reference Knowledge Layer

`references/` contains stable rules that should not be hidden inside prompts:
supported GEO formats, QC expectations, Seurat/Scanpy rules, annotation rules,
CellChat/pseudotime interpretation limits, environment setup, and workflow
contracts.

### Execution Script Layer

`scripts/` contains deterministic local tools:

- Python controllers for diagnosis, planning, rendering, validation, packaging,
  and report generation.
- `scripts/env_setup/` for cross-platform environment checks and installers.
- `scripts/course_adapted/` for runnable Seurat V5 modules.
- `scripts/course_source/` for traceable course-derived source references.

### Template Layer

`templates/` contains reusable Scanpy and report templates. Use templates for
portable or lightweight workflows. Use `scripts/course_adapted/` for
course-traced Seurat workflows.

### Validation And Evidence Layer

Validation is script-driven. `scripts/validate_result_bundle.py` checks local
analysis folders, while `scripts/validate_full_workflow.py` can validate a
full example directory when one is present locally. Large raw data, 10x
matrices, h5ad files, and Seurat RDS objects are intentionally excluded from
source control.

## Environment Strategy

The environment is split into profiles:

- `environment.yml` creates a cross-platform conda base with Python, Scanpy, and
  R base.
- `scripts/env_setup/install_r_packages.R --profile minimal` installs the R
  packages required by the core Seurat/SingleR workflow.
- `--profile extended` adds optional CellChat, Monocle3, and enrichment support.
- `--profile course` preserves the broad historical course package installer.

`scripts/env_setup/check_environment.py` is the preferred verification entrypoint
on Windows, Linux, and macOS.

## Data Flow

1. Diagnose input files or accessions.
2. Generate an analysis plan with data-readiness and biological-risk notes.
3. Render course-adapted R scripts or Scanpy templates with explicit paths.
4. Run local analysis outside the agent runtime.
5. Validate result folders and expected outputs.
6. Write reports that separate observations, statistical inference, and
   mechanism hypotheses.
7. Keep generated data, analysis objects, and private metadata out of source
   control.

## Safety Boundaries

The project supports clinical research hypothesis generation, not diagnosis or
treatment. Cell-cell communication, pseudotime, CNV, enrichment, and virtual
perturbation outputs must be described as computational inference unless
independently validated.

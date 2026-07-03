# Cross-Platform Environment Setup

This project supports Windows, Linux, and macOS through one minimal conda
environment plus profile-based R package installation.

## Profiles

| Profile | Purpose | Includes |
|---|---|---|
| `python` | Scanpy/helper-only checks | Python imports from `environment.yml` |
| `minimal` | Submission/example main workflow | Python stack, R base, Seurat, SingleR, celldex, plotting/table packages |
| `extended` | Optional course modules | `minimal` plus CellChat, Monocle3, clusterProfiler, organism annotation packages |
| `course` | Full course archive compatibility | The broad historical Seurat V5 course package set |

Use `minimal` by default. Use `extended` only when running
`scripts/course_adapted/02_marker_enrichment_from_seurat.R`,
`03_cellchat_from_seurat.R`, or `04_monocle3_from_seurat.R` on a local Seurat
object.

## Recommended Fresh Install

From the repository root:

```bash
conda env create -f environment.yml
conda activate scmechanism-agent
Rscript scripts/env_setup/install_r_packages.R --profile minimal
python scripts/env_setup/check_environment.py --profile minimal
```

To update an existing environment:

```bash
python scripts/env_setup/install_minimal_env.py --env-name scmechanism-agent
```

The installer runs `conda env update`, installs the `minimal` R profile, and
then verifies the environment.

## Optional Extended Modules

```bash
Rscript scripts/env_setup/install_r_packages.R --profile extended
python scripts/env_setup/check_environment.py --profile extended
```

The extended profile may take longer and can require compilers. Install it only
for workflows that actually call CellChat, Monocle3, or enrichment modules.

## Platform Notes

### Windows

- Install Miniconda or Mambaforge first.
- Install R 4.3 or newer, or use the `r-base` provided by `environment.yml`.
- Install Rtools matching the R major/minor version when R packages compile
  from source.
- PowerShell scripts in `scripts/env_setup/*.ps1` remain available for legacy
  Windows/course installs, but the Python entrypoints above are the preferred
  cross-platform path.

### Linux

- Use a recent conda distribution.
- If R source packages compile, ensure system build tools are available
  (`gcc`, `g++`, `gfortran`, `make`, development headers).
- Prefer the minimal profile on servers where package installs are managed by
  administrators.

### macOS

- Apple Silicon and Intel macOS are both supported by conda-forge for the
  Python/R base layer.
- If R packages compile from source, install Xcode Command Line Tools.
- Bioconductor data packages can be slow over unstable networks; rerunning the
  installer is safe because installed packages are skipped unless `--force` is
  used.

## Verification Commands

```bash
python scripts/env_setup/env_setup_selftest.py
python scripts/env_setup/check_environment.py --profile minimal
```

If the Seurat workflow script stops at a missing 10x folder, that means the
local data are absent, not that the environment is broken. The repository keeps
large data and RDS objects out of source control.

When a full local example directory is available, validate it with:

```bash
python scripts/validate_full_workflow.py --project-root . --example-root /path/to/example --out-md /tmp/scmechanism_full_workflow_validation.md
```

# Requirements Inventory

## Environment Profiles

Use `references/environment/cross-platform-setup.md` as the primary
Windows/Linux/macOS setup guide.

The recommended default is the `minimal` profile:

- `environment.yml` installs Python, Scanpy-compatible packages, and R base.
- `Rscript scripts/env_setup/install_r_packages.R --profile minimal` installs
  Seurat, SingleR, celldex, SingleCellExperiment, and plotting/table packages
  required by the core Seurat/SingleR workflow scripts.
- `python scripts/env_setup/check_environment.py --profile minimal` verifies
  the local environment.

Use `extended` only for optional downstream modules:

- marker enrichment: `clusterProfiler`, `org.Mm.eg.db`, `org.Hs.eg.db`
- CellChat: `CellChat`
- Monocle3: `monocle3`

Use `course` only when reproducing the broad historical Seurat V5 course stack.

## External Software

The course bundle includes installers in `1.软件安装/`. The skill requires R 4.3 or newer and Python 3.10 or newer:

- `R-4.4.2-win.exe`
- `rtools44-6335-6327.exe`
- `JAGS-4.3.1.exe`
- `RStudio-2024.12.0-467.exe`

Recommended installed paths:

- R: any R 4.3 or newer installation
- Rscript: the corresponding `bin\Rscript.exe`
- Rtools: Rtools43 for R 4.3, Rtools44 for R 4.4, or the matching newer Rtools release
- Conda: `C:\ProgramData\miniconda3`
- Git gzip: `C:\Program Files\Git\usr\bin\gzip.exe`

## R Package Sources

The complete installer uses four groups:

1. CRAN packages, including the Seurat core stack, plotting, IO, and utility packages.
2. Bioconductor packages, including annotation, enrichment, single-cell QC, CNV, and deconvolution dependencies.
3. GitHub packages explicitly listed in the course install script.
4. Local tar.gz packages bundled in `2.R包安装/`.

## Packages Imported by Course Scripts

The scripts import these R packages:

```text
AnnoProbe, ape, assertthat, BiocParallel, car, CCA, CellChat, celldex,
circlize, clusterProfiler, clustree, ComplexHeatmap, copykat, COSG,
cowplot, data.table, decontX, dendextend, doParallel, DOSE, DoubletFinder,
dplyr, e1071, fastSave, foreach, future, future.apply, ggalluvial, ggExtra,
ggforce, ggplot2, ggpubr, ggsci, ggthemes, glmGamPoi, grid, gridExtra,
GSEABase, GSVA, harmony, hdf5r, hdWGCNA, igraph, infercnv, irGSEA, limma,
Matrix, miscTools, monocle, monocle3, MuSiC, NMF, openai, org.Hs.eg.db,
parallel, patchwork, phylogram, plyr, preprocessCore, presto, qs, R.utils,
randomcoloR, remotes, reshape2, ridge, scales, scater, scattermore,
scCustomize, scDblFinder, SCINA, scPred, SCpubr, Seurat, SeuratObject,
SingleCellExperiment, SingleR, starTracer, stringr, tcltk, tidyverse,
tricycle, UCell, vioplot, WGCNA
```

Base/recommended packages such as `grid`, `parallel`, `tcltk`, and `Matrix` are still checked but may already ship with R or standard binary repositories.

## GitHub Packages

The course install script lists:

```text
GfellerLab/EPIC
Moonerss/CIBERSORT
dviraran/xCell
ebecht/MCPcounter, subdir Source
cansysbio/ConsensusTME
cit-bioinfo/mMCP-counter
omnideconv/immunedeconv
Hy4m/linkET
chuiqin/irGSEA
jinworks/CellChat
jokergoo/circlize
navinlabcode/copykat
chris-mcginnis-ucsf/DoubletFinder
immunogenomics/presto
campbio/decontX
powellgenomicslab/scPred at 9f407b7436f40d44224a5976a94cc6815c6e837f
genecell/COSGR
JerryZhang-1222/starTracer
cole-trapnell-lab/monocle3
smorabit/hdWGCNA at dev
xuranw/MuSiC
jackbibby1/SCPA
```

It also pins:

```text
crossmatch 1.3.1
multicross 2.1.0
```

## Course-Local R Archives

Install these from `2.R包安装/` when present:

```text
openai.tar.gz
hdWGCNA-0.3.00.tar.gz
```

`hdWGCNA-0.3.00.tar.gz` can be used as a local fallback if GitHub install fails or if the user wants the course-provided version.

## Python Packages

The course is mostly R, but reticulate/API/scverse-compatible modules benefit from a local Conda environment:

```text
python>=3.10
pandas
numpy
scipy
matplotlib
seaborn
scikit-learn
scanpy
anndata
scrublet
celltypist
gprofiler-official
umap-learn
igraph
leidenalg
openai
requests
jupyter
ipykernel
```

Use `seuratv5-course-py` as the default Conda environment name.

For platform-independent local work, this repository also provides:

```text
environment.yml
```

Create it with:

```powershell
conda env create -f environment.yml
conda activate scmechanism-agent
Rscript scripts/env_setup/install_r_packages.R --profile minimal
python scripts/env_setup/check_environment.py --profile minimal
```

This file is the cross-platform conda base layer. R package installation is
profile-based and handled by `scripts/env_setup/install_r_packages.R`.

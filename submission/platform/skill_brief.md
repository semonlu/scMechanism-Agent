# Skill Introduction

`scMechanism Agent` is a Medical AI Skill for public and local single-cell transcriptomics mechanism research.

It helps researchers complete a reproducible loop:

```text
clinical or biological question
  -> public/local data format diagnosis
  -> analysis route planning
  -> Seurat/Scanpy/CellChat/Monocle3 script selection
  -> result quality review
  -> cautious mechanism interpretation
  -> Methods/Results/report drafting
```

The skill is designed for GEO/SRA/local scRNA-seq workflows, especially 10x MEX, 10x H5, h5ad, Seurat RDS, marker tables, enrichment tables, CellChat results, and pseudotime outputs.

It does not replace a local bioinformatics runtime. Heavy computation can run locally, on a server, or on an approved analysis environment; the platform demonstration focuses on whether the Skill can be invoked, can choose the correct workflow, and can produce stable rule-based outputs.


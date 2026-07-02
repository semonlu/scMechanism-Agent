# GEO/SRA Input Diagnosis

## Detected Formats
- AnnData h5ad
- R/Seurat object

## Accessions
`{'GEO Series': [], 'GEO Sample': [], 'SRA Run': [], 'SRA Project': []}`

## Direct Downstream Analysis
`yes`

## Needs Fastq Rebuild
`False`

## Recommended Reading
- Use scanpy.read_h5ad; verify raw counts/layers before differential expression.
- Use readRDS/load and inspect class, assays, layers, reductions, and metadata.

## Major Risks
- Sample metadata/group design not detected.
- Processed objects may not retain raw counts needed for DE, CellChat, CNV, or deconvolution.

## Next Step
`Generate an analysis plan only after organism, metadata columns, and research comparison are known.`

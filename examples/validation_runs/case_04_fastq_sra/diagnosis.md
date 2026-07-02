# GEO/SRA Input Diagnosis

## Detected Formats
- FASTQ/SRA raw reads

## Accessions
`{'GEO Series': [], 'GEO Sample': [], 'SRA Run': ['SRR12345678', 'SRR12345679'], 'SRA Project': []}`

## Direct Downstream Analysis
`no`

## Needs Fastq Rebuild
`True`

## Recommended Reading
- Do not run downstream Seurat/Scanpy directly; first rebuild matrices with Cell Ranger, STARsolo, kallisto-bustools, or alevin-fry.

## Major Risks
- Sample metadata/group design not detected.
- Raw-read reconstruction requires reference genome, chemistry, and compute resources.

## Next Step
`Generate an analysis plan only after organism, metadata columns, and research comparison are known.`

# Security and Privacy

This skill is for public or locally owned single-cell research workflows and hypothesis generation. It is not a clinical diagnostic, treatment, or patient-management system.

## Data Handling

- Do not upload private patient-level clinical records, protected health information, API keys, tokens, or institutional secrets.
- Do not include raw FASTQ/SRA data, expression matrices, h5ad/loom files, Seurat RDS objects, or large private result tables in the repository.
- Keep heavy computation and private data on the user's local machine, institutional server, or approved compute environment.
- Use the Skill to guide format diagnosis, workflow planning, script generation, result review, and report drafting.

## Biological and Clinical Limits

- Treat CellChat, pseudotime, enrichment, CNV, deconvolution, and virtual perturbation outputs as computational inference.
- Do not present public single-cell reanalysis as direct clinical evidence for an individual patient.
- Require independent experimental, clinical, or cohort validation before making mechanistic or translational claims.

## External Services

- Do not send private matrices or identifiable clinical metadata to external LLM/API services.
- If an MCP or online search tool is used, use it for public metadata, documentation, literature, and accession lookup rather than private data transfer.

## Repository Hygiene

Keep generated local data, result folders, large matrices, credentials, and private metadata out of version control. Share only scripts, templates, references, and small non-sensitive examples.

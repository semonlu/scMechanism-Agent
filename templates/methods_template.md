# Methods Template

Public single-cell RNA-seq data were obtained from {{DATA_SOURCE}}. Input files were diagnosed as {{INPUT_FORMAT}}. Data were processed locally using {{PIPELINE}}. Cells were filtered after inspection of detected genes, total counts, mitochondrial percentage, and sample-level quality metrics. Normalization, variable feature detection, scaling, dimensionality reduction, neighborhood graph construction, clustering, and UMAP visualization were performed with recorded parameters in the analysis scripts.

Cell-type annotation was assigned using canonical marker genes and, when available, reference-based annotation tools. Differential expression and marker analyses were conducted within the documented comparison design. Functional enrichment was performed using organism-appropriate databases. Cell-cell communication and pseudotime analyses were run only when cell labels, sample design, and data structure met the assumptions described in the method-status tables.

All analyses were treated as exploratory public-data reanalysis and were not used for clinical diagnosis or treatment decision-making.

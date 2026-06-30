# Result Quality Check

Overall quality judgment: basically usable after manual biological review.

## Detected Files

- metadata: 1
- markers: 9
- deg: 0
- enrichment: 4
- cellchat: 5
- pseudotime: 14
- umap: 10

## Monocle3 Pseudotime Sanity Check

- trajectory cells: 5,000
- finite pseudotime cells: 5,000
- unique finite pseudotime values: 4,537
- cell subset: marker-supported mesenchymal trajectory states
- cell states: chondrocyte 820; cycling 762; osteoblast 292; tenocyte_fibroblast 3,126
- root rule: E15 cells
- stage distribution in sampled trajectory object: E15 1,813; P1 1,388; P7 915; P14 301; P28 583

## Main Risks

- No obvious structural output gaps were detected.
- The dataset is a mouse shoulder enthesis development dataset, not a frozen-shoulder disease-control cohort.
- SingleR labels and trajectory interpretation still require marker-level and tissue-context review before manuscript use.

## Manuscript Readiness

Use the outputs as a reproducible example and exploratory analysis package. Formal manuscript claims need replicated design, annotation review, statistical modeling, and manually selected biological contrasts.

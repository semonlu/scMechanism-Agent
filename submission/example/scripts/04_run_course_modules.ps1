param(
  [string]$Rscript = "Rscript",
  [string]$ProjectRoot = "",
  [string]$ExampleRoot = ""
)

$ErrorActionPreference = "Stop"

if (-not $ExampleRoot) {
  $ExampleRoot = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
}
if (-not $ProjectRoot) {
  $ProjectRoot = Split-Path -Parent (Split-Path -Parent $ExampleRoot)
}

$rendered = Join-Path $ExampleRoot "rendered_course_modules"
$results = Join-Path $ExampleRoot "results"
$coreObject = Join-Path $results "gse223751_seurat\objects\gse223751_processed_seurat.rds"
$courseResults = Join-Path $results "course_modules"

New-Item -ItemType Directory -Force -Path $rendered | Out-Null

$renderTemplate = Join-Path $ProjectRoot "scripts\render_template.py"

function Convert-ToRPath {
  param([string]$Path)
  return ($Path -replace '\\', '/')
}

$coreObjectR = Convert-ToRPath $coreObject
$markerOutputR = Convert-ToRPath (Join-Path $courseResults 'marker_enrichment')
$cellchatOutputR = Convert-ToRPath (Join-Path $courseResults 'cellchat')
$monocleOutputR = Convert-ToRPath (Join-Path $courseResults 'monocle3')

python $renderTemplate `
  --template (Join-Path $ProjectRoot "scripts\course_adapted\02_marker_enrichment_from_seurat.R") `
  --out (Join-Path $rendered "02_marker_enrichment_from_seurat.R") `
  --define "SEURAT_RDS=$coreObjectR" `
  --define "OUTPUT_DIR=$markerOutputR" `
  --define "CLUSTER_COL=singleR_label" `
  --define "ORGANISM=mouse" `
  --define "TOP_N=20"

python $renderTemplate `
  --template (Join-Path $ProjectRoot "scripts\course_adapted\03_cellchat_from_seurat.R") `
  --out (Join-Path $rendered "03_cellchat_from_seurat.R") `
  --define "SEURAT_RDS=$coreObjectR" `
  --define "OUTPUT_DIR=$cellchatOutputR" `
  --define "CELLTYPE_COL=singleR_label" `
  --define "ORGANISM=mouse" `
  --define "MIN_CELLS=50"

python $renderTemplate `
  --template (Join-Path $ProjectRoot "scripts\course_adapted\04_monocle3_from_seurat.R") `
  --out (Join-Path $rendered "04_monocle3_from_seurat.R") `
  --define "SEURAT_RDS=$coreObjectR" `
  --define "OUTPUT_DIR=$monocleOutputR" `
  --define "CELLTYPE_COL=marker_support_label" `
  --define "SUBSET_QUERY=marker_support_label %in% c('tenocyte_fibroblast', 'chondrocyte', 'osteoblast', 'cycling')" `
  --define "ROOT_CELLS_FILE=" `
  --define "ROOT_QUERY=grepl('^E15$', stage)" `
  --define "MAX_CELLS=5000" `
  --define "BALANCE_COL=marker_support_label" `
  --define "USE_SEURAT_UMAP=true" `
  --define "RANDOM_SEED=20260630"

& $Rscript (Join-Path $rendered "02_marker_enrichment_from_seurat.R")
& $Rscript (Join-Path $rendered "03_cellchat_from_seurat.R")
& $Rscript (Join-Path $rendered "04_monocle3_from_seurat.R")

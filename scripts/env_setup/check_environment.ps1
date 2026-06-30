param(
  [string]$CourseRoot = "",
  [string]$Rscript = "",
  [string]$CondaEnv = "seuratv5-course-py"
)

$ErrorActionPreference = "Stop"
$reportRoot = Join-Path (Get-Location) "seuratv5_environment_check"
$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$reportDir = Join-Path $reportRoot $timestamp
New-Item -ItemType Directory -Force -Path $reportDir | Out-Null

function Add-ReportLine {
  param([string]$File, [string]$Line)
  Add-Content -Path (Join-Path $reportDir $File) -Value $Line -Encoding UTF8
}

function Test-CommandOrPath {
  param([string]$Name, [string]$Path)
  $exists = $false
  $resolved = ""
  if ($Path -and (Test-Path -LiteralPath $Path)) {
    $exists = $true
    $resolved = (Resolve-Path -LiteralPath $Path).Path
  } else {
    $cmd = Get-Command $Name -ErrorAction SilentlyContinue
    if ($cmd) {
      $exists = $true
      $resolved = $cmd.Source
    }
  }
  Add-ReportLine -File "tools.tsv" -Line "$Name`t$exists`t$resolved"
  if ($exists) {
    Write-Host ("OK   {0}: {1}" -f $Name, $resolved)
  } else {
    Write-Host ("MISS {0}" -f $Name)
  }
}

Add-ReportLine -File "tools.tsv" -Line "tool`tavailable`tpath"
Test-CommandOrPath -Name "Rscript" -Path $Rscript
Test-CommandOrPath -Name "R" -Path ""
Test-CommandOrPath -Name "Rtools make" -Path "E:\rtools44\usr\bin\make.exe"
Test-CommandOrPath -Name "gzip" -Path "C:\Program Files\Git\usr\bin\gzip.exe"
Test-CommandOrPath -Name "conda" -Path "C:\ProgramData\miniconda3\Scripts\conda.exe"

foreach ($jags in @(
  "C:\Program Files\JAGS\JAGS-4.3.1\x64\bin\jags-terminal.exe",
  "C:\Program Files\JAGS\JAGS-4.3.1\bin\jags-terminal.exe"
)) {
  if (Test-Path -LiteralPath $jags) {
    Add-ReportLine -File "tools.tsv" -Line "JAGS`tTrue`t$jags"
    Write-Host "OK   JAGS: $jags"
    break
  }
}

if ($CourseRoot -and (Test-Path -LiteralPath $CourseRoot)) {
  Add-ReportLine -File "course.tsv" -Line "item`tavailable`tpath"
  Add-ReportLine -File "course.tsv" -Line "course_root`tTrue`t$CourseRoot"
  foreach ($fileName in @(
    "R-*.exe",
    "rtools44-6335-6327.exe",
    "rtools43*.exe",
    "JAGS-4.3.1.exe",
    "monocle_xxdchange.tar.gz",
    "openai.tar.gz",
    "xxdAIcelltype.tar.gz",
    "hdWGCNA-0.3.00.tar.gz"
  )) {
    $found = Get-ChildItem -LiteralPath $CourseRoot -Recurse -File -Filter $fileName -ErrorAction SilentlyContinue | Select-Object -First 1
    $foundPath = if ($found) { $found.FullName } else { "" }
    Add-ReportLine -File "course.tsv" -Line "$fileName`t$($found -ne $null)`t$foundPath"
  }
} else {
  Add-ReportLine -File "course.tsv" -Line "item`tavailable`tpath"
  Add-ReportLine -File "course.tsv" -Line "course_root`tFalse`t$CourseRoot"
}

$rscriptPath = $Rscript
if (-not $rscriptPath -or -not (Test-Path -LiteralPath $rscriptPath)) {
  $cmd = Get-Command Rscript -ErrorAction SilentlyContinue
  if ($cmd) { $rscriptPath = $cmd.Source }
}

if ($rscriptPath -and (Test-Path -LiteralPath $rscriptPath)) {
  $rVersionText = & $rscriptPath -e "cat(as.character(getRversion()))"
  Add-ReportLine -File "tools.tsv" -Line "R version`t$([version]$rVersionText -ge [version]'4.3.0')`t$rVersionText"
  if ([version]$rVersionText -lt [version]"4.3.0") {
    Write-Host "MISS R version: detected $rVersionText, requires >= 4.3"
  } else {
    Write-Host "OK   R version: $rVersionText (requires >= 4.3)"
  }
  $rCheck = @'
packages <- c(
  "AnnoProbe","ape","assertthat","BiocParallel","car","CCA","CellChat","celldex",
  "circlize","clusterProfiler","clustree","ComplexHeatmap","copykat","COSG",
  "cowplot","data.table","decontX","dendextend","doParallel","DOSE",
  "DoubletFinder","dplyr","e1071","fastSave","foreach","future","future.apply",
  "ggalluvial","ggExtra","ggforce","ggplot2","ggpubr","ggsci","ggthemes",
  "glmGamPoi","GSEABase","GSVA","harmony","hdf5r","hdWGCNA","igraph",
  "infercnv","irGSEA","limma","Matrix","miscTools","monocle","monocle3",
  "MuSiC","NMF","openai","org.Hs.eg.db","patchwork","phylogram","plyr",
  "preprocessCore","presto","qs","R.utils","randomcoloR","remotes","reshape2",
  "ridge","scales","scater","scattermore","scCustomize","scDblFinder","SCINA",
  "scPred","SCpubr","Seurat","SeuratObject","SingleCellExperiment","SingleR",
  "starTracer","stringr","tidyverse","tricycle","UCell","vioplot","WGCNA",
  "xxdAIcelltype"
)
out <- data.frame(package = packages, installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE))
out$version <- vapply(packages, function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else ""
}, character(1))
utils::write.table(out, file = commandArgs(trailingOnly = TRUE)[1], sep = "\t", quote = FALSE, row.names = FALSE)
cat("R package check complete: ", sum(out$installed), "/", nrow(out), " installed\n", sep = "")
missing <- out$package[!out$installed]
if (length(missing)) cat("Missing R packages: ", paste(missing, collapse = ", "), "\n", sep = "")
'@
  $rCheckPath = Join-Path $reportDir "check_r_packages.R"
  Set-Content -Path $rCheckPath -Value $rCheck -Encoding ASCII
  & $rscriptPath $rCheckPath (Join-Path $reportDir "r_packages.tsv")
} else {
  Write-Host "Skipping R package check because Rscript was not found."
}

$conda = Get-Command conda -ErrorAction SilentlyContinue
if (-not $conda -and (Test-Path -LiteralPath "C:\ProgramData\miniconda3\Scripts\conda.exe")) {
  $condaPath = "C:\ProgramData\miniconda3\Scripts\conda.exe"
} elseif ($conda) {
  $condaPath = $conda.Source
} else {
  $condaPath = ""
}

if ($condaPath) {
  $envList = & $condaPath env list
  $envExists = ($envList | Select-String -SimpleMatch $CondaEnv) -ne $null
  if (-not $envExists) {
    Write-Host "Skipping Python package check because conda env was not found: $CondaEnv"
    Write-Host "Create it with install_python_env.ps1 or install_environment.ps1 -InstallPythonEnv."
    Write-Host ""
    Write-Host "Environment check report: $reportDir"
    exit 0
  }

  $pyCheck = @'
import importlib.util
packages = {
    "pandas": "pandas",
    "numpy": "numpy",
    "scipy": "scipy",
    "matplotlib": "matplotlib",
    "seaborn": "seaborn",
    "sklearn": "scikit-learn",
    "scanpy": "scanpy",
    "anndata": "anndata",
    "scrublet": "scrublet",
    "celltypist": "celltypist",
    "gprofiler": "gprofiler-official",
    "umap": "umap-learn",
    "igraph": "igraph",
    "leidenalg": "leidenalg",
    "openai": "openai",
    "requests": "requests",
    "jupyter": "jupyter",
    "ipykernel": "ipykernel",
}
for module, name in packages.items():
    print(f"{name}\t{importlib.util.find_spec(module) is not None}")
'@
  $pyCheckPath = Join-Path $reportDir "check_python_packages.py"
  Set-Content -Path $pyCheckPath -Value $pyCheck -Encoding ASCII
  & $condaPath run -n $CondaEnv python $pyCheckPath | Tee-Object -FilePath (Join-Path $reportDir "python_packages.tsv")
} else {
  Write-Host "Skipping Python package check because conda was not found."
}

Write-Host ""
Write-Host "Environment check report: $reportDir"

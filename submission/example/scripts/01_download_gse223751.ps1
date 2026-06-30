param(
  [string]$OutDir = ".\data\raw",
  [string]$Url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223751/suppl/GSE223751_RAW.tar"
)

$ErrorActionPreference = "Stop"
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$outFile = Join-Path $OutDir "GSE223751_RAW.tar"

if (Test-Path $outFile) {
  Write-Host "Found existing file: $outFile"
  Get-Item $outFile | Select-Object FullName, Length, LastWriteTime
  exit 0
}

Write-Host "Downloading $Url"
curl.exe -L --fail --retry 3 --retry-delay 5 -o $outFile $Url
Get-Item $outFile | Select-Object FullName, Length, LastWriteTime

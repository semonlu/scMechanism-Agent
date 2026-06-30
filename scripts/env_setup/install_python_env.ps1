param(
  [string]$CondaEnv = "seuratv5-course-py",
  [string]$PythonVersion = "3.10",
  [switch]$DryRun
)

$ErrorActionPreference = "Stop"

if ([version]$PythonVersion -lt [version]"3.10") {
  throw "Python modules require Python >= 3.10. Requested: $PythonVersion"
}

$packages = @(
  "pandas",
  "numpy",
  "scipy",
  "matplotlib",
  "seaborn",
  "scikit-learn",
  "scanpy",
  "anndata",
  "scrublet",
  "celltypist",
  "gprofiler-official",
  "umap-learn",
  "igraph",
  "leidenalg",
  "openai",
  "requests",
  "jupyter",
  "ipykernel"
)

function Invoke-Logged {
  param([string]$Exe, [string[]]$Arguments)
  $cmd = "$Exe $($Arguments -join ' ')"
  if ($DryRun) {
    Write-Host "[DryRun] $cmd"
    return
  }
  Write-Host $cmd
  & $Exe @Arguments
  if ($LASTEXITCODE -ne 0) {
    throw "Command failed: $cmd"
  }
}

$conda = Get-Command conda -ErrorAction SilentlyContinue
if (-not $conda) {
  $candidate = "C:\ProgramData\miniconda3\Scripts\conda.exe"
  if (Test-Path -LiteralPath $candidate) {
    $condaPath = $candidate
  } else {
    throw "conda was not found. Install Miniconda/Anaconda or add conda to PATH."
  }
} else {
  $condaPath = $conda.Source
}

$envExists = $false
if (-not $DryRun) {
  $envList = & $condaPath env list
  $envExists = ($envList | Select-String -SimpleMatch $CondaEnv) -ne $null
}

if (-not $envExists) {
  Invoke-Logged -Exe $condaPath -Arguments @("create", "-n", $CondaEnv, "python=$PythonVersion", "-y")
} else {
  Write-Host "Conda environment already exists: $CondaEnv"
}

Invoke-Logged -Exe $condaPath -Arguments @("run", "-n", $CondaEnv, "python", "-m", "pip", "install", "--upgrade", "pip", "setuptools", "wheel")
Invoke-Logged -Exe $condaPath -Arguments (@("run", "-n", $CondaEnv, "python", "-m", "pip", "install") + $packages)
Invoke-Logged -Exe $condaPath -Arguments @("run", "-n", $CondaEnv, "python", "-m", "ipykernel", "install", "--user", "--name", $CondaEnv, "--display-name", "Python ($CondaEnv)")

Write-Host "Python environment is ready: $CondaEnv"

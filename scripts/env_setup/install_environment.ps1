param(
  [string]$CourseRoot = "",
  [string]$RHome = "",
  [string]$CondaEnv = "seuratv5-course-py",
  [switch]$InstallRPackages,
  [switch]$InstallPythonEnv,
  [switch]$DryRun
)

$ErrorActionPreference = "Stop"
$EnvSetupDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$SkillRoot = Split-Path -Parent (Split-Path -Parent $EnvSetupDir)

function Write-Step {
  param([string]$Message)
  Write-Host ""
  Write-Host "==> $Message"
}

function Resolve-CourseRoot {
  param([string]$InputPath)
  if ($InputPath -and (Test-Path -LiteralPath $InputPath)) {
    return (Resolve-Path -LiteralPath $InputPath).Path
  }

  $candidates = Get-ChildItem -LiteralPath $SkillRoot -Recurse -File -Filter "monocle_xxdchange.tar.gz" -ErrorAction SilentlyContinue |
    ForEach-Object { Split-Path -Parent (Split-Path -Parent $_.FullName) }

  if ($candidates.Count -gt 0) {
    return $candidates[0]
  }

  return ""
}

function Add-UserPathEntry {
  param([string]$PathEntry)
  if (-not $PathEntry -or -not (Test-Path -LiteralPath $PathEntry)) {
    Write-Host "Skip missing PATH entry: $PathEntry"
    return
  }

  $currentUserPath = [Environment]::GetEnvironmentVariable("Path", "User")
  if (-not $currentUserPath) { $currentUserPath = "" }
  $parts = $currentUserPath -split ";" | Where-Object { $_ -ne "" }
  $trimChars = [char[]]@("\", "/")
  $alreadyPresent = $parts | Where-Object { $_.TrimEnd($trimChars) -ieq $PathEntry.TrimEnd($trimChars) }

  if (-not $alreadyPresent) {
    $newPath = if ($currentUserPath) { "$PathEntry;$currentUserPath" } else { $PathEntry }
    if ($DryRun) {
      Write-Host "[DryRun] Add to user PATH: $PathEntry"
    } else {
      [Environment]::SetEnvironmentVariable("Path", $newPath, "User")
      Write-Host "Added to user PATH: $PathEntry"
    }
  } else {
    Write-Host "PATH already contains: $PathEntry"
  }

  $envParts = $env:Path -split ";" | Where-Object { $_ -ne "" }
  $envAlready = $envParts | Where-Object { $_.TrimEnd($trimChars) -ieq $PathEntry.TrimEnd($trimChars) }
  if (-not $envAlready) {
    $env:Path = "$PathEntry;$env:Path"
  }
}

function Resolve-RscriptPath {
  param([string]$InputRHome)

  if ($InputRHome) {
    $candidate = Join-Path $InputRHome "bin\Rscript.exe"
    if (Test-Path -LiteralPath $candidate) {
      return (Resolve-Path -LiteralPath $candidate).Path
    }
  }

  $cmd = Get-Command Rscript -ErrorAction SilentlyContinue
  if ($cmd) {
    return $cmd.Source
  }

  $candidates = @()
  foreach ($root in @("E:\", "C:\Program Files\R")) {
    if (Test-Path -LiteralPath $root) {
      $candidates += Get-ChildItem -LiteralPath $root -Directory -ErrorAction SilentlyContinue |
        Where-Object { $_.Name -match "^R-[0-9]+\.[0-9]+" } |
        Sort-Object Name -Descending |
        ForEach-Object { Join-Path $_.FullName "bin\Rscript.exe" }
    }
  }

  foreach ($candidate in $candidates) {
    if (Test-Path -LiteralPath $candidate) {
      return (Resolve-Path -LiteralPath $candidate).Path
    }
  }

  return ""
}

function Test-RVersion {
  param([string]$RscriptPath)
  if (-not $RscriptPath -or -not (Test-Path -LiteralPath $RscriptPath)) {
    return
  }
  $versionText = & $RscriptPath -e "cat(as.character(getRversion()))"
  if ([version]$versionText -lt [version]"4.3.0") {
    throw "R >= 4.3 is required, but detected R $versionText at $RscriptPath"
  }
  Write-Host "Detected R $versionText at $RscriptPath"
}

function Invoke-Script {
  param([string]$File, [string[]]$Arguments)
  $display = "powershell -NoProfile -ExecutionPolicy Bypass -File `"$File`" $($Arguments -join ' ')"
  if ($DryRun) {
    Write-Host "[DryRun] $display"
    return
  }
  & powershell -NoProfile -ExecutionPolicy Bypass -File $File @Arguments
}

$CourseRoot = Resolve-CourseRoot -InputPath $CourseRoot
if ($CourseRoot) {
  Write-Host "Course root: $CourseRoot"
} else {
  Write-Warning "Course root was not found. Pass -CourseRoot explicitly when installing local course archives."
}

$resolvedRscript = Resolve-RscriptPath -InputRHome $RHome
Test-RVersion -RscriptPath $resolvedRscript

Write-Step "Adding runtime directories to PATH"
if ($RHome) {
  Add-UserPathEntry -PathEntry (Join-Path $RHome "bin")
} elseif ($resolvedRscript) {
  Add-UserPathEntry -PathEntry (Split-Path -Parent $resolvedRscript)
} else {
  Write-Warning "Rscript was not found. Install R >= 4.3 or pass -RHome."
}

foreach ($root in @("E:\rtools44", "C:\rtools44", "E:\rtools43", "C:\rtools43")) {
  Add-UserPathEntry -PathEntry (Join-Path $root "usr\bin")
  Add-UserPathEntry -PathEntry (Join-Path $root "x86_64-w64-mingw32.static.posix\bin")
}

Add-UserPathEntry -PathEntry "C:\Program Files\Git\usr\bin"
Add-UserPathEntry -PathEntry "C:\ProgramData\miniconda3"
Add-UserPathEntry -PathEntry "C:\ProgramData\miniconda3\Scripts"
Add-UserPathEntry -PathEntry "C:\ProgramData\miniconda3\condabin"

foreach ($jags in @(
  "C:\Program Files\JAGS\JAGS-4.3.1\x64\bin",
  "C:\Program Files\JAGS\JAGS-4.3.1\bin"
)) {
  Add-UserPathEntry -PathEntry $jags
}

if ($InstallPythonEnv) {
  Write-Step "Installing Python environment"
  $args = @("-CondaEnv", $CondaEnv)
  if ($DryRun) { $args += "-DryRun" }
  Invoke-Script -File (Join-Path $EnvSetupDir "install_python_env.ps1") -Arguments $args
}

if ($InstallRPackages) {
  Write-Step "Installing R packages"
  $rscript = Resolve-RscriptPath -InputRHome $RHome
  if (-not (Test-Path -LiteralPath $rscript)) {
    throw "Rscript was not found. Install R >= 4.3, pass -RHome, or add Rscript to PATH."
  }
  Test-RVersion -RscriptPath $rscript
  $rArgs = @((Join-Path $EnvSetupDir "install_r_packages.R"))
  if ($CourseRoot) {
    $rArgs += "--course-root"
    $rArgs += $CourseRoot
  }
  if ($DryRun) {
    Write-Host "[DryRun] `"$rscript`" $($rArgs -join ' ')"
  } else {
    & $rscript @rArgs
  }
}

Write-Step "Checking environment"
$checkArgs = @("-CondaEnv", $CondaEnv)
if ($resolvedRscript) {
  $checkArgs += "-Rscript"
  $checkArgs += $resolvedRscript
}
if ($CourseRoot) {
  $checkArgs += "-CourseRoot"
  $checkArgs += $CourseRoot
}
Invoke-Script -File (Join-Path $EnvSetupDir "check_environment.ps1") -Arguments $checkArgs

Write-Host ""
Write-Host "Environment setup command completed. Open a new terminal so the updated user PATH is fully applied."

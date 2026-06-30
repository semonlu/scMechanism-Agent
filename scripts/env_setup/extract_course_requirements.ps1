param(
  [string]$CourseRoot = "",
  [string]$OutFile = "seuratv5_course_requirements.tsv"
)

$ErrorActionPreference = "Stop"

if (-not $CourseRoot -or -not (Test-Path -LiteralPath $CourseRoot)) {
  throw "Pass a valid -CourseRoot path."
}

$patterns = @(
  'library\(([^),]+)',
  'require\(([^),]+)',
  '::([A-Za-z][A-Za-z0-9_.]*)'
)

$packages = New-Object System.Collections.Generic.HashSet[string]
Get-ChildItem -LiteralPath $CourseRoot -Recurse -Include *.r,*.R -File | ForEach-Object {
  $text = Get-Content -LiteralPath $_.FullName -Raw
  foreach ($pattern in $patterns) {
    [regex]::Matches($text, $pattern) | ForEach-Object {
      $raw = $_.Groups[1].Value.Trim().Trim('"').Trim("'")
      if ($raw -match '^[A-Za-z][A-Za-z0-9_.]*$') {
        [void]$packages.Add($raw)
      }
    }
  }
}

$packages | Sort-Object | ForEach-Object { "$_" } | Set-Content -Path $OutFile -Encoding UTF8
Write-Host "Wrote package inventory: $OutFile"

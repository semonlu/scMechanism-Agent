# scMechanism Agent Local MCP Backend

This directory is an optional local extension. It lets the Medical AI Skill platform call a sandboxed backend running on your own computer through a Cloudflare Quick Tunnel.

The server is intentionally limited:

- no arbitrary shell command tool.
- no user-provided command execution.
- inputs are restricted to `mcp_server/workspace/inputs/`.
- outputs are restricted to `mcp_server/workspace/outputs/`.
- logs are written to `mcp_server/workspace/logs/`.
- platform connection tests should use `ping`, `list_available_pipelines`, `check_runtime_environment`, or GEO listing/downloading tools.
- only real-data workflow tools are exposed through MCP.

## Directory Layout

```text
mcp_server/
  local_mcp_server.py
  config.yaml
  README_LOCAL_MCP.md
  skill_runtime/
    SKILL.md
    agents/
    references/
    templates/
    scripts/
      course_source/
      course_adapted/
      env_setup/
  pipelines/
    run_scanpy_basic.py
    run_seurat_basic.R
    run_cellchat.R
    run_monocle3.R
  workspace/
    code/
    inputs/
    outputs/
    logs/
    runtime/
```

`skill_runtime/` is a copied local runtime mirror of the public Skill assets. It contains agents, references, workflow documents, templates, environment notes, cleaned course-source scripts, and runnable course-adapted scripts. The MCP server exposes this directory as read-only tools so the platform can inspect the workflow logic, while execution remains limited to whitelisted tools.

## Create A Conda Environment

```powershell
conda create -n scmechanism-mcp python=3.10 -y
conda activate scmechanism-mcp
```

The MCP server itself only needs the Python standard library. The Scanpy workflow should run in a scientific Python environment. This repository uses the existing course environment by default:

```powershell
conda activate seuratv5-course-py
python -m pip install scanpy anndata pandas numpy scipy matplotlib seaborn celltypist gseapy bbknn scanorama harmonypy==0.0.10 liana leidenalg igraph
```

By default the MCP server looks for:

```text
C:\ProgramData\miniconda3\envs\seuratv5-course-py\python.exe
```

You can override this without changing code:

```powershell
$env:SCMECHANISM_PYTHON = "C:\ProgramData\miniconda3\envs\seuratv5-course-py\python.exe"
python .\mcp_server\local_mcp_server.py --host 127.0.0.1 --port 8765
```

For Seurat execution, install R >= 4.3, add `Rscript` to PATH, and install Seurat packages using the main skill environment scripts:

```powershell
.\scripts\env_setup\install_environment.ps1
```

The local MCP can find `Rscript` from PATH, from the `RSCRIPT` environment variable, or from common Windows R locations such as `E:\R-4.4.2\bin\Rscript.exe`.

## Start The Local MCP Server

From the repository root:

```powershell
cd path\to\scMechanism-Agent
python .\mcp_server\local_mcp_server.py --host 127.0.0.1 --port 8765
```

Health endpoint:

```powershell
Invoke-RestMethod http://127.0.0.1:8765/health
```

List tools:

```powershell
Invoke-RestMethod http://127.0.0.1:8765/mcp
Invoke-RestMethod http://127.0.0.1:8765/mcp/
```

Both `/mcp` and `/mcp/` are registered explicitly and return the same JSON response. The server does not rely on trailing-slash redirects for MCP access.

## Platform Timeout Model

Many Medical AI Skill platforms terminate one MCP HTTP request after about 120 seconds. Real single-cell jobs often take much longer, so long-running work is exposed as asynchronous `start_*` tools:

```text
start_geo_download / start_extract_workspace_archive / start_seurat_basic / start_scanpy_basic / start_cellchat / start_monocle3
```

Each `start_*` call returns quickly with:

```text
status: submitted
job_id: <id>
log_file: mcp_server/workspace/logs/jobs/<id>.log
```

The platform should then call:

```text
get_job_status(job_id="<id>")
read_job_log(job_id="<id>")
list_result_files(project_id="<project>")
read_report(project_id="<project>")
```

Direct synchronous implementations still exist inside the backend for local developer diagnostics, but normal MCP discovery advertises the asynchronous versions so platform calls do not time out.

For any real analysis, also keep the data input synchronized:

```text
GEO data: start_geo_download -> start_extract_workspace_archive_if_needed -> validate_data_analysis_qc
Manual/platform input: register_input_dataset -> validate_data_analysis_qc
```

The backend writes `data_input_manifest.json` and `data_analysis_qc.md` under `workspace/outputs/<project_id>/`.

Call a JSON-RPC MCP-style tool:

```powershell
$body = @{
  jsonrpc = "2.0"
  id = 1
  method = "tools/call"
  params = @{
    name = "create_project"
    arguments = @{
      project_name = "case01"
    }
  }
} | ConvertTo-Json -Depth 6

Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8765/mcp -ContentType "application/json" -Body $body
```

List readable MCP resources from the workspace sandbox:

```powershell
$body = @{
  jsonrpc = "2.0"
  id = 10
  method = "resources/list"
  params = @{
    max_files = 200
  }
} | ConvertTo-Json -Depth 6

Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8765/mcp -ContentType "application/json" -Body $body
```

List real GEO supplementary files:

```powershell
$body = @{
  jsonrpc = "2.0"
  id = 2
  method = "tools/call"
  params = @{
    name = "list_geo_supplementary_files"
    arguments = @{
      gse_accession = "GSE176078"
    }
  }
} | ConvertTo-Json -Depth 6

Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8765/mcp -ContentType "application/json" -Body $body
```

Submit a real GEO supplementary-file download job into the sandbox:

```powershell
$body = @{
  jsonrpc = "2.0"
  id = 3
  method = "tools/call"
  params = @{
    name = "start_geo_download"
    arguments = @{
      project_id = "case01"
      gse_accession = "GSE176078"
      file_regex = "matrix|barcodes|features|genes|h5|h5ad|rds|tar|tgz"
      max_files = 10
    }
  }
} | ConvertTo-Json -Depth 6

Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8765/mcp -ContentType "application/json" -Body $body
```

The response returns immediately with a `job_id`. Poll the job until it finishes:

```powershell
$body = @{
  jsonrpc = "2.0"
  id = 4
  method = "tools/call"
  params = @{
    name = "get_job_status"
    arguments = @{
      job_id = "<returned-job-id>"
    }
  }
} | ConvertTo-Json -Depth 6

Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8765/mcp -ContentType "application/json" -Body $body
```

After downloading, extract archives with `start_extract_workspace_archive` when needed, inspect files under `mcp_server/workspace/inputs/case01/GSE176078/`, then run `start_seurat_basic` or `start_scanpy_basic` with the real input path.

## Expose With Cloudflare Quick Tunnel

Install `cloudflared` from Cloudflare, then run:

```powershell
cloudflared tunnel --url http://127.0.0.1:8765
```

Cloudflare prints a temporary HTTPS URL like:

```text
https://example-random.trycloudflare.com
```

Keep this terminal open while the Medical AI Skill platform needs access.

## Expose With Ngrok

If you use ngrok instead of Cloudflare, expose the same local port:

```powershell
ngrok http 8765
```

Ngrok prints a public HTTPS URL like:

```text
https://example.ngrok-free.dev
```

Verify that both MCP paths return HTTP 200 and not `307 Temporary Redirect`:

```powershell
curl.exe -i https://example.ngrok-free.dev/mcp
curl.exe -i https://example.ngrok-free.dev/mcp/
```

## Add To Medical AI Skill Platform MCP Management

1. Open the Medical AI Skill platform.
2. Go to MCP Management.
3. Add a new MCP server.
4. Use a clear name such as `scmechanism-local-mcp`.
5. Set the URL to the tunnel URL plus `/mcp`. Prefer the no-trailing-slash MCP endpoint:

```text
https://example-random.trycloudflare.com/mcp
```

For ngrok, use the same pattern:

```text
https://example.ngrok-free.dev/mcp
```

If the platform connection test fails, try the trailing-slash endpoint:

```text
https://example.ngrok-free.dev/mcp/
```

If the platform expects to discover endpoints from a root server URL, try the root tunnel URL:

```text
https://example.ngrok-free.dev
```

6. If the platform asks for a health URL, use:

```text
https://example-random.trycloudflare.com/health
```

7. Test connection.
8. Use only the whitelisted tools listed below.

## Available Tools

- `ping()`
- `list_available_pipelines()`
- `create_project(project_name: str)`
- `list_result_files(project_id: str)`
- `read_report(project_id: str)`
- `register_input_dataset(project_id: str, input_path: str, input_type: str = "", source_label: str = "", metadata_path: str = "")`
- `validate_data_analysis_qc(project_id: str, input_path: str, input_type: str = "", metadata_path: str = "", module: str = "analysis")`
- `get_job_status(job_id: str)`
- `read_job_log(job_id: str, max_bytes: int = 1000000)`
- `list_geo_supplementary_files(gse_accession: str, file_regex: str = "")`
- `start_geo_download(project_id: str, gse_accession: str, file_regex: str = "", max_files: int = 20, overwrite: bool = false, max_bytes_per_file: int = 2000000000)`
- `start_extract_workspace_archive(archive_path: str, output_dir: str = "", overwrite: bool = false, max_members: int = 20000, project_id: str = "archive_extract")`
- `list_skill_runtime_files(relative_path: str = "", max_files: int = 300)`
- `read_skill_runtime_file(relative_path: str, max_bytes: int = 1000000)`
- `check_runtime_environment(project_id: str = "environment_check")`
- `render_workflow_scripts(project_id: str, input_path: str, species: str = "human", input_type: str = "", metadata_path: str = "", sample_id: str = "", batch_col: str = "", condition_col: str = "", reference_name: str = "")`
- `start_scanpy_basic(project_id: str, input_path: str, species: str = "human", input_type: str = "", batch_col: str = "", metadata_path: str = "", sample_id: str = "", group_col: str = "", annotation_method: str = "marker_summary", timeout_seconds: int = 3600)`
- `start_seurat_basic(project_id: str, input_path: str, species: str = "human", input_type: str = "", metadata_path: str = "", sample_id: str = "", batch_col: str = "", condition_col: str = "", reference_name: str = "", run_annotation: bool = true, run_marker_enrichment: bool = true, timeout_seconds: int = 7200)`
- `propose_downstream_modules(project_id: str)`
- `validate_result_bundle(project_id: str)`
- `start_cellchat(project_id: str, approval_token: str, celltype_col: str = "singleR_label", species: str = "human", min_cells: int = 10, timeout_seconds: int = 7200)`
- `start_monocle3(project_id: str, approval_token: str, celltype_col: str = "singleR_label", subset_query: str = "", root_query: str = "", timeout_seconds: int = 7200)`
- `list_workspace_files(relative_path: str = "", max_files: int = 200)`
- `read_workspace_file(relative_path: str, max_bytes: int = 1000000)`
- `write_workspace_file(relative_path: str, content: str = "", content_base64: str = "", overwrite: bool = false)`
- `download_to_workspace(url: str, relative_path: str, overwrite: bool = false, max_bytes: int = 100000000)`
- `run_workspace_python(project_id: str, script_path: str, args: list[str] = [], timeout_seconds: int = 120)`

`start_scanpy_basic` submits a guarded Scanpy workflow for `h5ad`, `10x_mtx`, or `10x_h5` inputs under `workspace/inputs/`. It uses the analysis Python detected by `SCMECHANISM_PYTHON` or the default `seuratv5-course-py` environment. If no metadata is supplied for a 10x folder, it creates default `sample_id`, `group`, and `batch` columns so downstream marker and review steps still have a stable grouping column. If metadata is supplied, pass `metadata_path` under `workspace/inputs/`; row names, `cell`, or `barcode` are matched to cell barcodes. It writes QC metrics, UMAP figures, marker tables, annotation-evidence tables, optional enrichment, logs, and a processed h5ad under `workspace/outputs/<project_id>/scanpy_basic/`.

`start_seurat_basic` submits the copied course-adapted Seurat V5 workflow on real input under `workspace/inputs/`. It supports `10x_mtx`, `10x_nonstandard`, `10x_h5`, `rds`, and `csv`:

1. `01_seurat_v5_core_pipeline.R`
2. `05_singler_cell_annotation.R`
3. `02_marker_enrichment_from_seurat.R`
4. `validate_result_bundle.py`
5. `propose_downstream_modules.py`

It does not automatically run CellChat or Monocle3.

Before analysis starts, `start_seurat_basic` and `start_scanpy_basic` write `data_analysis_qc.md`. If a project manifest exists and the requested `input_path` does not match the downloaded/extracted/registered input, the job stops before running biological analysis.

## Required Downstream Approval Gate

CellChat and pseudotime/Monocle3 must not run immediately after clustering. The MCP first generates:

```text
workspace/outputs/<project_id>/downstream_proposal.md
```

The user should review the proposed microenvironment, target cell lineages, cell type column, root/state assumptions, and warnings. Only then call:

```text
start_cellchat(..., approval_token="APPROVED_CELLCHAT")
start_monocle3(..., approval_token="APPROVED_MONOCLE3")
```

Without the matching approval token, both tools return `approval_required` and do not run.

## Minimal Local Workflow

1. List/download GEO data or put local data under `mcp_server/workspace/inputs/`.
2. Create a project:

```text
create_project(project_name="case01")
```

3. Check packages:

```text
check_runtime_environment(project_id="case01")
```

4. Run Seurat or Scanpy:

```text
list_geo_supplementary_files(gse_accession="GSE176078")
start_geo_download(project_id="case01", gse_accession="GSE176078", file_regex="scRNASeq|RAW|matrix|h5ad|rds|tar")
get_job_status(job_id="<download-job-id>")
read_job_log(job_id="<download-job-id>")
start_extract_workspace_archive(archive_path="case01/GSE176078/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz", output_dir="case01/GSE176078/extracted", project_id="case01")
get_job_status(job_id="<extract-job-id>")
validate_data_analysis_qc(project_id="case01", input_path="case01/GSE176078/extracted_sample", input_type="10x_nonstandard")
start_seurat_basic(project_id="case01", input_path="case01/GSE176078/extracted_sample", input_type="10x_nonstandard", species="human")
start_scanpy_basic(project_id="case01", input_path="processed.h5ad", input_type="h5ad", species="human")
start_scanpy_basic(project_id="case02", input_path="sample_10x", input_type="10x_mtx", species="human", sample_id="sample_10x")
```

5. Read `downstream_proposal.md`.
6. If approved, run CellChat or Monocle3 with the required approval token.

## Workspace File Sandbox

The platform can save, read, and process files only inside:

```text
mcp_server/workspace/
```

The MCP server exposes workspace files in two ways:

- tools: `list_workspace_files`, `read_workspace_file`, `write_workspace_file`, `download_to_workspace`, and `run_workspace_python`.
- resources: JSON-RPC `resources/list` and `resources/read`, using `sandbox://workspace/<relative-path>` URIs.

Allowed workspace areas:

- `workspace/inputs/`: downloaded or uploaded input files.
- `workspace/code/`: Python scripts that the platform may write and run.
- `workspace/outputs/`: project result folders.
- `workspace/logs/`: MCP audit logs.

All file tools reject absolute paths, `..` traversal, and path hints such as Desktop, Downloads, browser cache, WeChat, hospital folders, API keys, or token files.

Download rules:

- only `http://` and `https://` URLs are accepted.
- localhost, private IPs, link-local IPs, reserved IPs, and `.local` hosts are rejected.
- downloads are capped by `max_bytes`.

Script-running rules:

- only Python scripts under `workspace/code/` can be executed.
- the server does not expose `run_command(cmd)`.
- script arguments cannot contain absolute paths or parent-directory traversal.
- a static safety scan blocks scripts that import or call shell/network/private-path primitives such as `subprocess`, `os.system`, `socket`, `requests`, `urllib`, `ctypes`, or `shutil.rmtree`.
- a runtime Python guard restricts file operations to `workspace/` plus the Python installation needed to import standard libraries.
- stdout/stderr and return code are returned as JSON and written into the selected project log folder.

## Stop The Service

Press `Ctrl+C` in the terminal running `local_mcp_server.py`.

Stop the Cloudflare tunnel by pressing `Ctrl+C` in the `cloudflared` terminal.

## Security Notes

- Do not place patient-private data in `workspace/inputs/`.
- Do not expose this server longer than necessary.
- Do not share the Cloudflare tunnel URL publicly.
- Do not add a tool that executes arbitrary shell commands.
- Do not read Desktop, Downloads, browser cache, WeChat files, hospital folders, API keys, or token files.
- Do not upload or process protected health information.
- Do not delete system files.

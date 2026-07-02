# scMechanism Agent Local MCP Backend

This directory is an optional local extension. It lets the Medical AI Skill platform call a lightweight backend running on your own computer through a Cloudflare Quick Tunnel.

The server is intentionally limited:

- no arbitrary shell command tool.
- no user-provided command execution.
- inputs are restricted to `mcp_server/workspace/inputs/`.
- outputs are restricted to `mcp_server/workspace/outputs/`.
- logs are written to `mcp_server/workspace/logs/`.
- demo output is synthetic and contains no patient data.

## Directory Layout

```text
mcp_server/
  local_mcp_server.py
  config.yaml
  README_LOCAL_MCP.md
  pipelines/
    run_demo_pipeline.py
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

## Create A Conda Environment

```powershell
conda create -n scmechanism-mcp python=3.10 -y
conda activate scmechanism-mcp
```

The first MCP version only needs the Python standard library for the server and demo pipeline. Install optional scientific packages only when you are ready to enable real local execution:

```powershell
pip install scanpy pandas numpy scipy matplotlib seaborn
```

For Seurat execution, install R >= 4.3, add `Rscript` to PATH, and install Seurat packages using the main skill environment scripts.

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

Call a JSON-RPC MCP-style tool:

```powershell
$body = @{
  jsonrpc = "2.0"
  id = 1
  method = "tools/call"
  params = @{
    name = "create_project"
    arguments = @{
      project_name = "demo_project"
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

Run demo pipeline:

```powershell
$body = @{
  jsonrpc = "2.0"
  id = 2
  method = "tools/call"
  params = @{
    name = "run_demo_pipeline"
    arguments = @{
      project_id = "demo_project"
    }
  }
} | ConvertTo-Json -Depth 6

Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8765/mcp -ContentType "application/json" -Body $body
```

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
- `run_demo_pipeline(project_id: str)`
- `list_result_files(project_id: str)`
- `read_report(project_id: str)`
- `run_scanpy_basic(project_id: str, input_path: str, species: str = "human")`
- `run_seurat_basic(project_id: str, input_path: str, species: str = "human")`
- `list_workspace_files(relative_path: str = "", max_files: int = 200)`
- `read_workspace_file(relative_path: str, max_bytes: int = 1000000)`
- `write_workspace_file(relative_path: str, content: str = "", content_base64: str = "", overwrite: bool = false)`
- `download_to_workspace(url: str, relative_path: str, overwrite: bool = false, max_bytes: int = 100000000)`
- `run_workspace_python(project_id: str, script_path: str, args: list[str] = [], timeout_seconds: int = 120)`

`run_scanpy_basic` and `run_seurat_basic` are guarded placeholders in the first version. They validate input location and dependencies, then return clear JSON errors if dependencies or inputs are missing.

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

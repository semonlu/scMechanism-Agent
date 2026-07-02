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
    inputs/
    outputs/
    logs/
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
```

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

## Add To Medical AI Skill Platform MCP Management

1. Open the Medical AI Skill platform.
2. Go to MCP Management.
3. Add a new MCP server.
4. Use a clear name such as `scmechanism-local-mcp`.
5. Set the URL to the Cloudflare tunnel URL plus `/mcp`, for example:

```text
https://example-random.trycloudflare.com/mcp
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

`run_scanpy_basic` and `run_seurat_basic` are guarded placeholders in the first version. They validate input location and dependencies, then return clear JSON errors if dependencies or inputs are missing.

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


#!/usr/bin/env python3
"""Local optional MCP-like HTTP backend for lightweight scRNA-seq tasks.

This server intentionally exposes only a small whitelist of tools. It never
executes user-provided shell commands and restricts user input files to
workspace/inputs and outputs to workspace/outputs.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import re
import shutil
import sys
import time
import uuid
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Callable
from urllib.parse import parse_qs, urlparse


SERVER_ROOT = Path(__file__).resolve().parent
WORKSPACE_ROOT = SERVER_ROOT / "workspace"
INPUTS_DIR = WORKSPACE_ROOT / "inputs"
OUTPUTS_DIR = WORKSPACE_ROOT / "outputs"
LOGS_DIR = WORKSPACE_ROOT / "logs"

PIPELINES = [
    "scanpy_basic",
    "seurat_basic",
    "cellchat",
    "monocle3",
    "demo_pipeline",
    "report_skeleton",
]

FORBIDDEN_PATH_HINTS = [
    "desktop",
    "downloads",
    "wechat",
    "weixin",
    "browser",
    "cache",
    "token",
    "api_key",
    "hospital",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def ensure_workspace() -> None:
    for path in [INPUTS_DIR, OUTPUTS_DIR, LOGS_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def safe_slug(value: str, default: str = "project") -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", (value or "").strip()).strip("._-")
    return slug[:80] or default


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(parent.resolve(strict=False))
        return True
    except ValueError:
        return False


def require_project_dir(project_id: str, create: bool = False) -> Path:
    project_slug = safe_slug(project_id, "project")
    project_dir = OUTPUTS_DIR / project_slug
    if not is_relative_to(project_dir, OUTPUTS_DIR):
        raise ValueError("Invalid project_id")
    if create:
        project_dir.mkdir(parents=True, exist_ok=True)
    if not project_dir.exists():
        raise FileNotFoundError(f"Project not found: {project_slug}")
    return project_dir


def validate_input_path(input_path: str) -> Path:
    if not input_path:
        raise ValueError("input_path is required")
    lowered = input_path.lower()
    if any(hint in lowered for hint in FORBIDDEN_PATH_HINTS):
        raise ValueError("input_path appears to reference a forbidden private/system location")
    candidate = Path(input_path)
    if not candidate.is_absolute():
        candidate = INPUTS_DIR / candidate
    candidate = candidate.resolve(strict=False)
    if not is_relative_to(candidate, INPUTS_DIR):
        raise ValueError("input_path must be inside mcp_server/workspace/inputs")
    return candidate


def log_event(tool: str, job_id: str, payload: dict[str, Any]) -> Path:
    ensure_workspace()
    log_file = LOGS_DIR / f"{safe_slug(tool)}-{job_id}.log"
    line = json.dumps({"time": utc_now(), "tool": tool, "job_id": job_id, **payload}, ensure_ascii=False)
    with log_file.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")
    return log_file


def response_payload(
    *,
    status: str,
    message: str,
    project_id: str | None = None,
    job_id: str | None = None,
    output_dir: Path | str | None = None,
    log_file: Path | str | None = None,
    warnings: list[str] | None = None,
    data: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "status": status,
        "project_id": project_id,
        "job_id": job_id,
        "output_dir": str(output_dir) if output_dir else None,
        "log_file": str(log_file) if log_file else None,
        "message": message,
        "warnings": warnings or [],
        "data": data or {},
    }


def tool_ping() -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    log_file = log_event("ping", job_id, {"status": "ok"})
    return response_payload(
        status="ok",
        project_id=None,
        job_id=job_id,
        output_dir=OUTPUTS_DIR,
        log_file=log_file,
        message="Local scMechanism MCP server is running.",
        data={
            "time_utc": utc_now(),
            "server_root": str(SERVER_ROOT),
            "workspace": str(WORKSPACE_ROOT),
            "python": sys.version.split()[0],
        },
    )


def tool_list_available_pipelines() -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    log_file = log_event("list_available_pipelines", job_id, {"pipelines": PIPELINES})
    return response_payload(
        status="ok",
        project_id=None,
        job_id=job_id,
        output_dir=OUTPUTS_DIR,
        log_file=log_file,
        message="Available whitelisted pipelines returned.",
        data={"pipelines": PIPELINES},
    )


def tool_create_project(project_name: str) -> dict[str, Any]:
    project_id = safe_slug(project_name, f"project_{int(time.time())}")
    job_id = uuid.uuid4().hex
    project_dir = require_project_dir(project_id, create=True)
    for sub in ["tables", "figures", "objects", "logs"]:
        (project_dir / sub).mkdir(exist_ok=True)
    log_file = log_event("create_project", job_id, {"project_id": project_id, "output_dir": str(project_dir)})
    return response_payload(
        status="ok",
        project_id=project_id,
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Project directory created under workspace/outputs.",
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def tool_run_demo_pipeline(project_id: str) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    project_dir = require_project_dir(project_id, create=True)
    tables_dir = project_dir / "tables"

    write_csv(
        tables_dir / "marker_genes.csv",
        [
            {"cluster": "0", "cell_type_hint": "Macrophage", "gene": "LYZ", "avg_log2FC": 2.1, "p_val_adj": 0.001},
            {"cluster": "1", "cell_type_hint": "Fibroblast", "gene": "COL1A1", "avg_log2FC": 2.4, "p_val_adj": 0.002},
            {"cluster": "2", "cell_type_hint": "T cell", "gene": "CD3D", "avg_log2FC": 1.8, "p_val_adj": 0.005},
        ],
    )
    write_csv(
        tables_dir / "deg_results.csv",
        [
            {"gene": "IL1B", "cell_type": "Macrophage", "comparison": "disease_vs_control", "log2FC": 1.5, "p_val_adj": 0.01},
            {"gene": "COL3A1", "cell_type": "Fibroblast", "comparison": "disease_vs_control", "log2FC": 1.2, "p_val_adj": 0.03},
        ],
    )
    write_csv(
        tables_dir / "enrichment_results.csv",
        [
            {"term": "inflammatory response", "p_adj": 0.01, "genes": "IL1B;LYZ"},
            {"term": "extracellular matrix organization", "p_adj": 0.02, "genes": "COL1A1;COL3A1"},
        ],
    )
    write_csv(
        tables_dir / "celltype_proportion.csv",
        [
            {"sample_id": "S1", "condition": "disease", "cell_type": "Macrophage", "proportion": 0.32},
            {"sample_id": "S2", "condition": "control", "cell_type": "Macrophage", "proportion": 0.18},
            {"sample_id": "S1", "condition": "disease", "cell_type": "Fibroblast", "proportion": 0.41},
            {"sample_id": "S2", "condition": "control", "cell_type": "Fibroblast", "proportion": 0.28},
        ],
    )
    report = project_dir / "report_skeleton.md"
    report.write_text(
        """# Demo scRNA-seq Mechanism Report Skeleton

## Summary

This demo result is synthetic and contains no patient data. It is only used to test platform-to-local MCP connectivity.

## Observed Patterns

- Macrophage markers and IL1B suggest an inflammatory signal.
- Fibroblast markers and extracellular-matrix enrichment suggest stromal remodeling.
- Cell type proportions show a larger macrophage/fibroblast fraction in the demo disease sample.

## Safety

These are generated demo files, not clinical evidence.
""".strip()
        + "\n",
        encoding="utf-8",
    )
    run_log = project_dir / "run_log.txt"
    run_log.write_text(f"{utc_now()} demo pipeline completed for {project_id}\n", encoding="utf-8")
    log_file = log_event("run_demo_pipeline", job_id, {"project_id": project_id, "status": "ok"})
    return response_payload(
        status="ok",
        project_id=project_id,
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Demo pipeline generated lightweight result files.",
        data={"files": [p.name for p in sorted(tables_dir.glob("*.csv"))] + [report.name, run_log.name]},
    )


def tool_list_result_files(project_id: str) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    project_dir = require_project_dir(project_id, create=False)
    files = [str(p.relative_to(project_dir)).replace("\\", "/") for p in sorted(project_dir.rglob("*")) if p.is_file()]
    log_file = log_event("list_result_files", job_id, {"project_id": project_id, "n_files": len(files)})
    return response_payload(
        status="ok",
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Result files listed.",
        data={"files": files},
    )


def tool_read_report(project_id: str) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    project_dir = require_project_dir(project_id, create=False)
    report = project_dir / "report_skeleton.md"
    if not report.exists():
        log_file = log_event("read_report", job_id, {"project_id": project_id, "status": "missing_report"})
        return response_payload(
            status="error",
            project_id=safe_slug(project_id),
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message="report_skeleton.md not found. Run run_demo_pipeline first or create the report.",
            warnings=["No report file exists for this project."],
        )
    text = report.read_text(encoding="utf-8")
    log_file = log_event("read_report", job_id, {"project_id": project_id, "status": "ok"})
    return response_payload(
        status="ok",
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Report read successfully.",
        data={"report": text},
    )


def tool_run_scanpy_basic(project_id: str, input_path: str, species: str = "human") -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    project_dir = require_project_dir(project_id, create=True)
    warnings: list[str] = []
    try:
        candidate = validate_input_path(input_path)
    except Exception as exc:
        log_file = log_event("run_scanpy_basic", job_id, {"project_id": project_id, "status": "invalid_input", "error": str(exc)})
        return response_payload(
            status="error",
            project_id=safe_slug(project_id),
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message=str(exc),
            warnings=["All input paths must stay under mcp_server/workspace/inputs."],
        )
    if not candidate.exists():
        warnings.append("Input path does not exist under workspace/inputs.")
    if importlib.util.find_spec("scanpy") is None:
        warnings.append("Python package scanpy is not installed in the current environment.")
    log_file = log_event("run_scanpy_basic", job_id, {"project_id": project_id, "status": "not_run", "input": str(candidate), "species": species})
    return response_payload(
        status="error",
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="scanpy_basic is a guarded placeholder in this MCP version; install dependencies and place input under workspace/inputs before enabling execution.",
        warnings=warnings or ["Execution intentionally disabled in first MCP version; interface validated only."],
    )


def tool_run_seurat_basic(project_id: str, input_path: str, species: str = "human") -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    project_dir = require_project_dir(project_id, create=True)
    warnings: list[str] = []
    try:
        candidate = validate_input_path(input_path)
    except Exception as exc:
        log_file = log_event("run_seurat_basic", job_id, {"project_id": project_id, "status": "invalid_input", "error": str(exc)})
        return response_payload(
            status="error",
            project_id=safe_slug(project_id),
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message=str(exc),
            warnings=["All input paths must stay under mcp_server/workspace/inputs."],
        )
    if not candidate.exists():
        warnings.append("Input path does not exist under workspace/inputs.")
    if shutil.which("Rscript") is None:
        warnings.append("Rscript is not on PATH. Install R >= 4.3 and add Rscript to PATH.")
    warnings.append("Seurat execution is disabled in this first MCP version; use the main skill environment installer before enabling.")
    log_file = log_event("run_seurat_basic", job_id, {"project_id": project_id, "status": "not_run", "input": str(candidate), "species": species})
    return response_payload(
        status="error",
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="seurat_basic is a guarded placeholder. Install R/Seurat and place input under workspace/inputs before enabling execution.",
        warnings=warnings,
    )


TOOLS: dict[str, Callable[..., dict[str, Any]]] = {
    "ping": tool_ping,
    "list_available_pipelines": tool_list_available_pipelines,
    "create_project": tool_create_project,
    "run_demo_pipeline": tool_run_demo_pipeline,
    "list_result_files": tool_list_result_files,
    "read_report": tool_read_report,
    "run_scanpy_basic": tool_run_scanpy_basic,
    "run_seurat_basic": tool_run_seurat_basic,
}


TOOL_SCHEMAS = [
    {"name": "ping", "description": "Return service status, current time, and workspace paths.", "inputSchema": {"type": "object", "properties": {}}},
    {"name": "list_available_pipelines", "description": "List whitelisted local pipelines.", "inputSchema": {"type": "object", "properties": {}}},
    {"name": "create_project", "description": "Create a project under workspace/outputs.", "inputSchema": {"type": "object", "properties": {"project_name": {"type": "string"}}, "required": ["project_name"]}},
    {"name": "run_demo_pipeline", "description": "Generate lightweight synthetic result files without private data.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}}, "required": ["project_id"]}},
    {"name": "list_result_files", "description": "List files under one project output directory.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}}, "required": ["project_id"]}},
    {"name": "read_report", "description": "Read report_skeleton.md from one project.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}}, "required": ["project_id"]}},
    {"name": "run_scanpy_basic", "description": "Guarded placeholder for Scanpy execution.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}, "input_path": {"type": "string"}, "species": {"type": "string", "default": "human"}}, "required": ["project_id", "input_path"]}},
    {"name": "run_seurat_basic", "description": "Guarded placeholder for Seurat execution.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}, "input_path": {"type": "string"}, "species": {"type": "string", "default": "human"}}, "required": ["project_id", "input_path"]}},
]


def call_tool(name: str, arguments: dict[str, Any] | None = None) -> dict[str, Any]:
    if name not in TOOLS:
        return response_payload(
            status="error",
            project_id=(arguments or {}).get("project_id"),
            job_id=uuid.uuid4().hex,
            output_dir=None,
            log_file=None,
            message=f"Unknown or disallowed tool: {name}",
            warnings=["Only whitelisted tools are available. Arbitrary command execution is not implemented."],
        )
    try:
        return TOOLS[name](**(arguments or {}))
    except TypeError as exc:
        return response_payload(status="error", project_id=(arguments or {}).get("project_id"), job_id=uuid.uuid4().hex, message=f"Invalid arguments for {name}: {exc}", warnings=[])
    except Exception as exc:
        job_id = uuid.uuid4().hex
        log_file = log_event(name, job_id, {"status": "error", "error": str(exc)})
        return response_payload(status="error", project_id=(arguments or {}).get("project_id"), job_id=job_id, log_file=log_file, message=str(exc), warnings=[])


class LocalMCPHandler(BaseHTTPRequestHandler):
    server_version = "scMechanismLocalMCP/0.1"

    def _send_json(self, payload: Any, status: HTTPStatus = HTTPStatus.OK) -> None:
        body = json.dumps(payload, ensure_ascii=False, indent=2).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(body)

    def _read_json(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0"))
        if length <= 0:
            return {}
        raw = self.rfile.read(length).decode("utf-8")
        return json.loads(raw or "{}")

    def do_OPTIONS(self) -> None:
        self.send_response(HTTPStatus.NO_CONTENT)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "content-type, authorization")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.end_headers()

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/health":
            self._send_json(tool_ping())
            return
        if parsed.path == "/mcp":
            query = parse_qs(parsed.query)
            if "tool" in query:
                name = query["tool"][0]
                args = {k: v[0] for k, v in query.items() if k != "tool"}
                self._send_json(call_tool(name, args))
                return
            self._send_json({"status": "ok", "message": "MCP endpoint is available. Use POST /mcp for JSON-RPC tools/call.", "tools": TOOL_SCHEMAS})
            return
        if parsed.path == "/sse":
            body = (
                "event: health\n"
                f"data: {json.dumps({'status': 'ok', 'time_utc': utc_now()}, ensure_ascii=False)}\n\n"
                "event: tools\n"
                f"data: {json.dumps({'tools': [tool['name'] for tool in TOOL_SCHEMAS]}, ensure_ascii=False)}\n\n"
            ).encode("utf-8")
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "text/event-stream; charset=utf-8")
            self.send_header("Cache-Control", "no-cache")
            self.send_header("Access-Control-Allow-Origin", "*")
            self.end_headers()
            self.wfile.write(body)
            return
        self._send_json({"status": "error", "message": "Not found"}, HTTPStatus.NOT_FOUND)

    def do_POST(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path != "/mcp":
            self._send_json({"status": "error", "message": "Not found"}, HTTPStatus.NOT_FOUND)
            return
        try:
            payload = self._read_json()
        except Exception as exc:
            self._send_json({"status": "error", "message": f"Invalid JSON: {exc}"}, HTTPStatus.BAD_REQUEST)
            return

        if payload.get("jsonrpc") == "2.0":
            req_id = payload.get("id")
            method = payload.get("method")
            params = payload.get("params") or {}
            if method == "initialize":
                result = {
                    "protocolVersion": "2024-11-05",
                    "capabilities": {"tools": {}},
                    "serverInfo": {"name": "scmechanism-local-mcp", "version": "0.1.0"},
                }
                self._send_json({"jsonrpc": "2.0", "id": req_id, "result": result})
                return
            if method == "tools/list":
                self._send_json({"jsonrpc": "2.0", "id": req_id, "result": {"tools": TOOL_SCHEMAS}})
                return
            if method == "tools/call":
                tool_name = params.get("name")
                arguments = params.get("arguments") or {}
                result = call_tool(tool_name, arguments)
                self._send_json({"jsonrpc": "2.0", "id": req_id, "result": {"content": [{"type": "text", "text": json.dumps(result, ensure_ascii=False)}], "isError": result.get("status") == "error"}})
                return
            self._send_json({"jsonrpc": "2.0", "id": req_id, "error": {"code": -32601, "message": f"Unsupported method: {method}"}})
            return

        tool_name = payload.get("tool") or payload.get("name")
        arguments = payload.get("arguments") or payload.get("params") or {}
        self._send_json(call_tool(tool_name, arguments))

    def log_message(self, fmt: str, *args: Any) -> None:
        LOGS_DIR.mkdir(parents=True, exist_ok=True)
        with (LOGS_DIR / "http_access.log").open("a", encoding="utf-8") as handle:
            handle.write(f"{utc_now()} {self.client_address[0]} {fmt % args}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    ensure_workspace()
    server = ThreadingHTTPServer((args.host, args.port), LocalMCPHandler)
    print(f"scMechanism local MCP server listening on http://{args.host}:{args.port}")
    print("Health: /health  MCP: /mcp  SSE: /sse")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping server...")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


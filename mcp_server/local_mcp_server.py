#!/usr/bin/env python3
"""Local optional MCP-like HTTP backend for lightweight scRNA-seq tasks.

This server intentionally exposes only a small whitelist of tools. It never
executes user-provided shell commands and restricts user input files to
workspace/inputs and outputs to workspace/outputs.
"""

from __future__ import annotations

import argparse
import base64
import csv
import importlib.util
import ipaddress
import json
import mimetypes
import os
import re
import shutil
import socket
import subprocess
import sys
import time
import uuid
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Callable
from urllib.parse import parse_qs, urlparse
from urllib.request import Request, urlopen


SERVER_ROOT = Path(__file__).resolve().parent
WORKSPACE_ROOT = SERVER_ROOT / "workspace"
INPUTS_DIR = WORKSPACE_ROOT / "inputs"
OUTPUTS_DIR = WORKSPACE_ROOT / "outputs"
LOGS_DIR = WORKSPACE_ROOT / "logs"
CODE_DIR = WORKSPACE_ROOT / "code"
RUNTIME_DIR = WORKSPACE_ROOT / "runtime"

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

TEXT_FILE_SUFFIXES = {
    ".txt",
    ".tsv",
    ".csv",
    ".json",
    ".md",
    ".log",
    ".yaml",
    ".yml",
    ".py",
    ".r",
}

MAX_READ_BYTES = 1_000_000
MAX_WRITE_BYTES = 5_000_000
MAX_DOWNLOAD_BYTES = 100_000_000
MAX_RUN_TIMEOUT_SECONDS = 300

FORBIDDEN_CODE_PATTERNS = [
    r"\bsubprocess\b",
    r"\bos\.system\b",
    r"\bos\.popen\b",
    r"\bshutil\.rmtree\b",
    r"\bctypes\b",
    r"\bsocket\b",
    r"\brequests\b",
    r"\burllib\b",
    r"\bhttpx\b",
    r"\bftplib\b",
    r"\bparamiko\b",
    r"\bwinreg\b",
    r"\bwebbrowser\b",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def ensure_workspace() -> None:
    for path in [INPUTS_DIR, OUTPUTS_DIR, LOGS_DIR, CODE_DIR, RUNTIME_DIR]:
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


def workspace_relative(path: Path) -> str:
    return str(path.relative_to(WORKSPACE_ROOT)).replace("\\", "/")


def validate_workspace_path(relative_path: str, *, must_exist: bool = False, allow_dir: bool = True) -> Path:
    if not relative_path:
        raise ValueError("relative_path is required")
    lowered = relative_path.lower()
    if any(hint in lowered for hint in FORBIDDEN_PATH_HINTS):
        raise ValueError("relative_path appears to reference a forbidden private/system location")
    candidate = Path(relative_path)
    if candidate.is_absolute():
        raise ValueError("Only relative paths under mcp_server/workspace are allowed")
    resolved = (WORKSPACE_ROOT / candidate).resolve(strict=False)
    if not is_relative_to(resolved, WORKSPACE_ROOT):
        raise ValueError("Path must stay inside mcp_server/workspace")
    if must_exist and not resolved.exists():
        raise FileNotFoundError(f"Workspace path not found: {relative_path}")
    if resolved.exists() and resolved.is_dir() and not allow_dir:
        raise ValueError("Expected a file path, got a directory")
    return resolved


def validate_code_path(relative_path: str, *, must_exist: bool = False) -> Path:
    path = validate_workspace_path(relative_path, must_exist=must_exist, allow_dir=False)
    if not is_relative_to(path, CODE_DIR):
        raise ValueError("Executable scripts must be inside mcp_server/workspace/code")
    if path.suffix.lower() != ".py":
        raise ValueError("Only Python .py scripts are executable in this MCP version")
    return path


def validate_download_url(url: str) -> str:
    parsed = urlparse(url)
    if parsed.scheme not in {"http", "https"}:
        raise ValueError("Only http and https downloads are allowed")
    if not parsed.hostname:
        raise ValueError("Download URL must include a hostname")
    hostname = parsed.hostname.lower()
    if hostname in {"localhost"} or hostname.endswith(".local"):
        raise ValueError("Downloads from localhost/private hosts are not allowed")
    try:
        infos = socket.getaddrinfo(hostname, parsed.port or (443 if parsed.scheme == "https" else 80), type=socket.SOCK_STREAM)
    except socket.gaierror as exc:
        raise ValueError(f"Could not resolve download host: {exc}") from exc
    for info in infos:
        ip = ipaddress.ip_address(info[4][0])
        if ip.is_private or ip.is_loopback or ip.is_link_local or ip.is_multicast or ip.is_reserved:
            raise ValueError("Downloads from private, local, or reserved network addresses are not allowed")
    return url


def scan_code_for_sandbox_violations(script_path: Path) -> list[str]:
    text = script_path.read_text(encoding="utf-8", errors="replace")
    violations = []
    for pattern in FORBIDDEN_CODE_PATTERNS:
        if re.search(pattern, text):
            violations.append(pattern.replace("\\b", "").replace("\\", ""))
    lowered = text.lower()
    if any(hint in lowered for hint in FORBIDDEN_PATH_HINTS):
        violations.append("forbidden private path hint")
    if re.search(r"[A-Za-z]:\\", text) or re.search(r"(^|[^.])/(users|home|etc|var|root)\b", text, flags=re.IGNORECASE):
        violations.append("absolute filesystem path")
    return sorted(set(violations))


def write_python_sandbox_guard() -> Path:
    guard_dir = RUNTIME_DIR / "python_guard"
    guard_dir.mkdir(parents=True, exist_ok=True)
    guard_file = guard_dir / "sitecustomize.py"
    guard_file.write_text(
        f'''
import builtins
import os
import sys
from pathlib import Path

_WORKSPACE = Path(os.environ.get("SCMECHANISM_WORKSPACE", r"{WORKSPACE_ROOT}")).resolve()
_ALLOWED_ROOTS = [
    _WORKSPACE,
    Path(sys.base_prefix).resolve(),
    Path(sys.prefix).resolve(),
    Path(sys.exec_prefix).resolve(),
    Path(sys.executable).resolve().parent,
]
_orig_open = builtins.open
_orig_os_open = os.open
_orig_remove = os.remove
_orig_unlink = os.unlink
_orig_rename = os.rename
_orig_replace = os.replace
_orig_mkdir = os.mkdir
_orig_makedirs = os.makedirs
_orig_listdir = os.listdir
_orig_stat = os.stat
_orig_scandir = os.scandir

def _check_path(path):
    if path is None:
        return path
    try:
        resolved = Path(path).resolve()
    except Exception:
        resolved = (_WORKSPACE / Path(path)).resolve()
    for root in _ALLOWED_ROOTS:
        try:
            resolved.relative_to(root)
            return path
        except ValueError:
            pass
    raise PermissionError("Path escapes mcp_server/workspace sandbox: " + str(path))
    return path

def open(file, *args, **kwargs):
    _check_path(file)
    return _orig_open(file, *args, **kwargs)

builtins.open = open

def _wrap_one_arg(fn):
    def wrapped(path, *args, **kwargs):
        _check_path(path)
        return fn(path, *args, **kwargs)
    return wrapped

def _wrap_two_arg(fn):
    def wrapped(src, dst, *args, **kwargs):
        _check_path(src)
        _check_path(dst)
        return fn(src, dst, *args, **kwargs)
    return wrapped

os.open = _wrap_one_arg(_orig_os_open)
os.remove = _wrap_one_arg(_orig_remove)
os.unlink = _wrap_one_arg(_orig_unlink)
os.mkdir = _wrap_one_arg(_orig_mkdir)
os.makedirs = _wrap_one_arg(_orig_makedirs)
os.listdir = _wrap_one_arg(_orig_listdir)
os.stat = _wrap_one_arg(_orig_stat)
os.scandir = _wrap_one_arg(_orig_scandir)
os.rename = _wrap_two_arg(_orig_rename)
os.replace = _wrap_two_arg(_orig_replace)

_orig_path_open = Path.open
_orig_path_read_text = Path.read_text
_orig_path_read_bytes = Path.read_bytes
_orig_path_write_text = Path.write_text
_orig_path_write_bytes = Path.write_bytes
_orig_path_unlink = Path.unlink
_orig_path_mkdir = Path.mkdir
_orig_path_iterdir = Path.iterdir

def _path_checked(self):
    _check_path(self)
    return self

def _path_open(self, *args, **kwargs):
    return _orig_path_open(_path_checked(self), *args, **kwargs)

Path.open = _path_open
Path.read_text = lambda self, *a, **k: _orig_path_read_text(_path_checked(self), *a, **k)
Path.read_bytes = lambda self, *a, **k: _orig_path_read_bytes(_path_checked(self), *a, **k)
Path.write_text = lambda self, *a, **k: _orig_path_write_text(_path_checked(self), *a, **k)
Path.write_bytes = lambda self, *a, **k: _orig_path_write_bytes(_path_checked(self), *a, **k)
Path.unlink = lambda self, *a, **k: _orig_path_unlink(_path_checked(self), *a, **k)
Path.mkdir = lambda self, *a, **k: _orig_path_mkdir(_path_checked(self), *a, **k)
Path.iterdir = lambda self: _orig_path_iterdir(_path_checked(self))
'''.lstrip(),
        encoding="utf-8",
    )
    return guard_dir


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


def tool_list_workspace_files(relative_path: str = "", max_files: int = 200) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    root = WORKSPACE_ROOT if not relative_path else validate_workspace_path(relative_path, must_exist=True, allow_dir=True)
    if root.is_file():
        files = [{
            "path": workspace_relative(root),
            "size": root.stat().st_size,
            "is_dir": False,
            "modified_utc": datetime.fromtimestamp(root.stat().st_mtime, timezone.utc).isoformat(),
        }]
    else:
        entries = []
        for path in sorted(root.rglob("*")):
            if path.is_file():
                stat = path.stat()
                entries.append({
                    "path": workspace_relative(path),
                    "size": stat.st_size,
                    "is_dir": False,
                    "modified_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                })
            if len(entries) >= max_files:
                break
        files = entries
    log_file = log_event("list_workspace_files", job_id, {"relative_path": relative_path, "n_files": len(files)})
    return response_payload(
        status="ok",
        job_id=job_id,
        output_dir=WORKSPACE_ROOT,
        log_file=log_file,
        message="Workspace files listed.",
        warnings=["Results are limited by max_files."] if len(files) >= max_files else [],
        data={"files": files, "workspace": str(WORKSPACE_ROOT)},
    )


def tool_read_workspace_file(relative_path: str, max_bytes: int = MAX_READ_BYTES) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    path = validate_workspace_path(relative_path, must_exist=True, allow_dir=False)
    max_bytes = max(1, min(int(max_bytes), MAX_READ_BYTES))
    raw = path.read_bytes()
    truncated = len(raw) > max_bytes
    raw_slice = raw[:max_bytes]
    is_text = path.suffix.lower() in TEXT_FILE_SUFFIXES
    data: dict[str, Any] = {
        "path": workspace_relative(path),
        "size": len(raw),
        "truncated": truncated,
        "encoding": "utf-8" if is_text else "base64",
    }
    if is_text:
        data["content"] = raw_slice.decode("utf-8", errors="replace")
    else:
        data["content_base64"] = base64.b64encode(raw_slice).decode("ascii")
    log_file = log_event("read_workspace_file", job_id, {"relative_path": relative_path, "bytes": len(raw_slice), "truncated": truncated})
    return response_payload(
        status="ok",
        job_id=job_id,
        output_dir=WORKSPACE_ROOT,
        log_file=log_file,
        message="Workspace file read.",
        warnings=["File was truncated to max_bytes."] if truncated else [],
        data=data,
    )


def tool_write_workspace_file(
    relative_path: str,
    content: str = "",
    content_base64: str = "",
    overwrite: bool = False,
) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    path = validate_workspace_path(relative_path, must_exist=False, allow_dir=False)
    if path.exists() and not overwrite:
        raise FileExistsError("File exists. Set overwrite=true to replace it.")
    if content_base64:
        raw = base64.b64decode(content_base64)
    else:
        raw = content.encode("utf-8")
    if len(raw) > MAX_WRITE_BYTES:
        raise ValueError(f"File content exceeds MAX_WRITE_BYTES={MAX_WRITE_BYTES}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)
    log_file = log_event("write_workspace_file", job_id, {"relative_path": relative_path, "bytes": len(raw)})
    return response_payload(
        status="ok",
        job_id=job_id,
        output_dir=path.parent,
        log_file=log_file,
        message="Workspace file written.",
        data={"path": workspace_relative(path), "size": len(raw)},
    )


def tool_download_to_workspace(url: str, relative_path: str, overwrite: bool = False, max_bytes: int = MAX_DOWNLOAD_BYTES) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    safe_url = validate_download_url(url)
    path = validate_workspace_path(relative_path, must_exist=False, allow_dir=False)
    if path.exists() and not overwrite:
        raise FileExistsError("File exists. Set overwrite=true to replace it.")
    max_bytes = max(1, min(int(max_bytes), MAX_DOWNLOAD_BYTES))
    request = Request(safe_url, headers={"User-Agent": "scMechanism-local-mcp/0.1"})
    path.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    with urlopen(request, timeout=60) as response, path.open("wb") as handle:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            total += len(chunk)
            if total > max_bytes:
                handle.close()
                path.unlink(missing_ok=True)
                raise ValueError(f"Download exceeds max_bytes={max_bytes}")
            handle.write(chunk)
    log_file = log_event("download_to_workspace", job_id, {"url": safe_url, "relative_path": relative_path, "bytes": total})
    return response_payload(
        status="ok",
        job_id=job_id,
        output_dir=path.parent,
        log_file=log_file,
        message="Downloaded file into workspace.",
        data={"path": workspace_relative(path), "size": total},
    )


def tool_run_workspace_python(
    project_id: str,
    script_path: str,
    args: list[str] | None = None,
    timeout_seconds: int = 120,
) -> dict[str, Any]:
    project_dir = require_project_dir(project_id, create=True)
    job_id = uuid.uuid4().hex
    script = validate_code_path(script_path, must_exist=True)
    violations = scan_code_for_sandbox_violations(script)
    if violations:
        log_file = log_event("run_workspace_python", job_id, {"project_id": project_id, "status": "blocked", "violations": violations})
        return response_payload(
            status="error",
            project_id=safe_slug(project_id),
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message="Script blocked by sandbox safety scan.",
            warnings=violations,
        )
    safe_args = [str(arg) for arg in (args or [])]
    if any(Path(arg).is_absolute() or ".." in Path(arg).parts for arg in safe_args):
        raise ValueError("Script args must not contain absolute paths or parent-directory traversal")
    timeout_seconds = max(1, min(int(timeout_seconds), MAX_RUN_TIMEOUT_SECONDS))
    guard_dir = write_python_sandbox_guard()
    env = os.environ.copy()
    existing_pythonpath = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(guard_dir) if not existing_pythonpath else str(guard_dir) + os.pathsep + existing_pythonpath
    env["SCMECHANISM_WORKSPACE"] = str(WORKSPACE_ROOT)
    env["SCMECHANISM_INPUTS"] = str(INPUTS_DIR)
    env["SCMECHANISM_OUTPUTS"] = str(OUTPUTS_DIR)
    env["SCMECHANISM_PROJECT_DIR"] = str(project_dir)
    env["SCMECHANISM_LOGS"] = str(LOGS_DIR)
    env["PYTHONIOENCODING"] = "utf-8"
    completed = subprocess.run(
        [sys.executable, str(script), *safe_args],
        cwd=WORKSPACE_ROOT,
        env=env,
        text=True,
        capture_output=True,
        timeout=timeout_seconds,
        shell=False,
    )
    run_log = project_dir / "logs" / f"run_workspace_python-{job_id}.log"
    run_log.parent.mkdir(parents=True, exist_ok=True)
    run_log.write_text(
        json.dumps(
            {
                "time": utc_now(),
                "script": workspace_relative(script),
                "returncode": completed.returncode,
                "stdout": completed.stdout[-10000:],
                "stderr": completed.stderr[-10000:],
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    log_file = log_event("run_workspace_python", job_id, {"project_id": project_id, "status": "ok" if completed.returncode == 0 else "error", "returncode": completed.returncode})
    return response_payload(
        status="ok" if completed.returncode == 0 else "error",
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Workspace Python script completed." if completed.returncode == 0 else "Workspace Python script failed.",
        warnings=[] if completed.returncode == 0 else ["See stdout/stderr in data and project log."],
        data={
            "script": workspace_relative(script),
            "returncode": completed.returncode,
            "stdout": completed.stdout[-10000:],
            "stderr": completed.stderr[-10000:],
            "run_log": workspace_relative(run_log),
        },
    )


def workspace_resource_uri(path: Path) -> str:
    return "sandbox://workspace/" + workspace_relative(path)


def path_from_workspace_resource_uri(uri: str) -> Path:
    prefix = "sandbox://workspace/"
    if not uri.startswith(prefix):
        raise ValueError("Only sandbox://workspace/ resource URIs are supported")
    return validate_workspace_path(uri[len(prefix):], must_exist=True, allow_dir=False)


def list_workspace_resources(max_files: int = 200) -> list[dict[str, Any]]:
    resources = []
    for path in sorted(WORKSPACE_ROOT.rglob("*")):
        if not path.is_file():
            continue
        if is_relative_to(path, RUNTIME_DIR):
            continue
        stat = path.stat()
        mime_type = mimetypes.guess_type(path.name)[0] or ("text/plain" if path.suffix.lower() in TEXT_FILE_SUFFIXES else "application/octet-stream")
        resources.append({
            "uri": workspace_resource_uri(path),
            "name": workspace_relative(path),
            "mimeType": mime_type,
            "size": stat.st_size,
        })
        if len(resources) >= max_files:
            break
    return resources


def read_workspace_resource(uri: str, max_bytes: int = MAX_READ_BYTES) -> dict[str, Any]:
    path = path_from_workspace_resource_uri(uri)
    max_bytes = max(1, min(int(max_bytes), MAX_READ_BYTES))
    raw = path.read_bytes()
    truncated = len(raw) > max_bytes
    raw_slice = raw[:max_bytes]
    mime_type = mimetypes.guess_type(path.name)[0] or ("text/plain" if path.suffix.lower() in TEXT_FILE_SUFFIXES else "application/octet-stream")
    content: dict[str, Any] = {"uri": uri, "mimeType": mime_type}
    if path.suffix.lower() in TEXT_FILE_SUFFIXES:
        content["text"] = raw_slice.decode("utf-8", errors="replace")
    else:
        content["blob"] = base64.b64encode(raw_slice).decode("ascii")
    return {"contents": [content], "truncated": truncated, "size": len(raw)}


TOOLS: dict[str, Callable[..., dict[str, Any]]] = {
    "ping": tool_ping,
    "list_available_pipelines": tool_list_available_pipelines,
    "create_project": tool_create_project,
    "run_demo_pipeline": tool_run_demo_pipeline,
    "list_result_files": tool_list_result_files,
    "read_report": tool_read_report,
    "run_scanpy_basic": tool_run_scanpy_basic,
    "run_seurat_basic": tool_run_seurat_basic,
    "list_workspace_files": tool_list_workspace_files,
    "read_workspace_file": tool_read_workspace_file,
    "write_workspace_file": tool_write_workspace_file,
    "download_to_workspace": tool_download_to_workspace,
    "run_workspace_python": tool_run_workspace_python,
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
    {"name": "list_workspace_files", "description": "List files under the MCP workspace sandbox.", "inputSchema": {"type": "object", "properties": {"relative_path": {"type": "string", "default": ""}, "max_files": {"type": "integer", "default": 200}}}},
    {"name": "read_workspace_file", "description": "Read a text or base64-encoded file from the MCP workspace sandbox.", "inputSchema": {"type": "object", "properties": {"relative_path": {"type": "string"}, "max_bytes": {"type": "integer", "default": MAX_READ_BYTES}}, "required": ["relative_path"]}},
    {"name": "write_workspace_file", "description": "Write text or base64 content into the MCP workspace sandbox.", "inputSchema": {"type": "object", "properties": {"relative_path": {"type": "string"}, "content": {"type": "string", "default": ""}, "content_base64": {"type": "string", "default": ""}, "overwrite": {"type": "boolean", "default": False}}, "required": ["relative_path"]}},
    {"name": "download_to_workspace", "description": "Download an http/https file into the MCP workspace sandbox.", "inputSchema": {"type": "object", "properties": {"url": {"type": "string"}, "relative_path": {"type": "string"}, "overwrite": {"type": "boolean", "default": False}, "max_bytes": {"type": "integer", "default": MAX_DOWNLOAD_BYTES}}, "required": ["url", "relative_path"]}},
    {"name": "run_workspace_python", "description": "Run a Python script located under workspace/code with workspace-only file access guards.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}, "script_path": {"type": "string"}, "args": {"type": "array", "items": {"type": "string"}}, "timeout_seconds": {"type": "integer", "default": 120}}, "required": ["project_id", "script_path"]}},
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
    HEALTH_PATHS = {"/health", "/health/"}
    MCP_PATHS = {"/mcp", "/mcp/"}

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
        if parsed.path in self.HEALTH_PATHS:
            self._send_json(tool_ping())
            return
        if parsed.path in self.MCP_PATHS:
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
        if parsed.path not in self.MCP_PATHS:
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
                    "capabilities": {"tools": {}, "resources": {}},
                    "serverInfo": {"name": "scmechanism-local-mcp", "version": "0.1.0"},
                }
                self._send_json({"jsonrpc": "2.0", "id": req_id, "result": result})
                return
            if method == "tools/list":
                self._send_json({"jsonrpc": "2.0", "id": req_id, "result": {"tools": TOOL_SCHEMAS}})
                return
            if method == "resources/list":
                max_files = int(params.get("max_files", 200)) if isinstance(params, dict) else 200
                self._send_json({"jsonrpc": "2.0", "id": req_id, "result": {"resources": list_workspace_resources(max_files=max_files)}})
                return
            if method == "resources/read":
                uri = params.get("uri") if isinstance(params, dict) else None
                if not uri:
                    self._send_json({"jsonrpc": "2.0", "id": req_id, "error": {"code": -32602, "message": "resources/read requires params.uri"}})
                    return
                try:
                    result = read_workspace_resource(uri, int(params.get("max_bytes", MAX_READ_BYTES)))
                except Exception as exc:
                    self._send_json({"jsonrpc": "2.0", "id": req_id, "error": {"code": -32000, "message": str(exc)}})
                    return
                self._send_json({"jsonrpc": "2.0", "id": req_id, "result": result})
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

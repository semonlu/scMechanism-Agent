#!/usr/bin/env python3
"""Local optional MCP-like HTTP backend for sandboxed scRNA-seq tasks.

This server intentionally exposes only a small whitelist of tools. It never
executes user-provided shell commands and restricts user input files to
workspace/inputs and outputs to workspace/outputs.
"""

from __future__ import annotations

import argparse
import base64
import csv
import gzip
import html
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
import tarfile
import time
import uuid
import zipfile
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Callable
from urllib.parse import parse_qs, urljoin, urlparse, unquote
from urllib.request import Request, urlopen


SERVER_ROOT = Path(__file__).resolve().parent
SKILL_RUNTIME_DIR = SERVER_ROOT / "skill_runtime"
WORKSPACE_ROOT = SERVER_ROOT / "workspace"
INPUTS_DIR = WORKSPACE_ROOT / "inputs"
OUTPUTS_DIR = WORKSPACE_ROOT / "outputs"
LOGS_DIR = WORKSPACE_ROOT / "logs"
CODE_DIR = WORKSPACE_ROOT / "code"
RUNTIME_DIR = WORKSPACE_ROOT / "runtime"
DEFAULT_ANALYSIS_PYTHON = r"C:\ProgramData\miniconda3\envs\seuratv5-course-py\python.exe"

PIPELINES = [
    "scanpy_basic",
    "seurat_basic",
    "seurat_full_workflow",
    "singler_annotation",
    "marker_enrichment",
    "downstream_proposal",
    "cellchat",
    "monocle3",
    "report_skeleton",
    "geo_download",
]

SYNTHETIC_TEST_PIPELINES = [
    "demo_pipeline",
    "full_workflow_selftest",
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
MAX_DOWNLOAD_BYTES = 2_000_000_000
MAX_RUN_TIMEOUT_SECONDS = 7200
GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/"

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


def geo_series_bucket(gse_accession: str) -> str:
    match = re.fullmatch(r"GSE(\d+)", gse_accession.strip().upper())
    if not match:
        raise ValueError("gse_accession must look like GSE12345")
    digits = match.group(1)
    if len(digits) <= 3:
        bucket = "GSEnnn"
    else:
        bucket = f"GSE{digits[:-3]}nnn"
    return bucket


def geo_suppl_url(gse_accession: str) -> str:
    gse = gse_accession.strip().upper()
    return urljoin(GEO_FTP_BASE, f"{geo_series_bucket(gse)}/{gse}/suppl/")


def fetch_geo_supplementary_listing(gse_accession: str) -> list[dict[str, Any]]:
    base_url = geo_suppl_url(gse_accession)
    safe_url = validate_download_url(base_url)
    request = Request(safe_url, headers={"User-Agent": "scMechanism-local-mcp/0.1"})
    with urlopen(request, timeout=120) as response:
        text = response.read().decode("utf-8", errors="replace")
    files: list[dict[str, Any]] = []
    seen: set[str] = set()
    for href in re.findall(r'href=["\']([^"\']+)["\']', text, flags=re.IGNORECASE):
        href = html.unescape(href)
        if href.startswith("?") or href.startswith("/") or href in {"../", "./"}:
            continue
        full_url = urljoin(base_url, href)
        if not full_url.startswith(base_url):
            continue
        filename = Path(unquote(href.rstrip("/"))).name
        if not filename or filename in seen:
            continue
        seen.add(filename)
        files.append({"filename": filename, "url": full_url})
    return files


def download_url_to_workspace(url: str, relative_path: str, overwrite: bool = False, max_bytes: int = MAX_DOWNLOAD_BYTES) -> tuple[Path, int, int | None]:
    safe_url = validate_download_url(url)
    path = validate_workspace_path(relative_path, must_exist=False, allow_dir=False)
    if path.exists() and not overwrite:
        raise FileExistsError("File exists. Set overwrite=true to replace it.")
    max_bytes = max(1, min(int(max_bytes), MAX_DOWNLOAD_BYTES))
    path.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    expected_size: int | None = None
    attempts = 3
    last_error = ""
    for attempt in range(1, attempts + 1):
        total = 0
        request = Request(safe_url, headers={"User-Agent": "scMechanism-local-mcp/0.1"})
        try:
            with urlopen(request, timeout=120) as response, path.open("wb") as handle:
                content_length = response.headers.get("Content-Length")
                expected_size = int(content_length) if content_length and content_length.isdigit() else None
                if expected_size is not None and expected_size > max_bytes:
                    handle.close()
                    path.unlink(missing_ok=True)
                    raise ValueError(f"Download exceeds max_bytes={max_bytes}; remote size={expected_size}")
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
            if expected_size is None or total == expected_size:
                break
            last_error = f"incomplete download: got {total} bytes, expected {expected_size}"
            path.unlink(missing_ok=True)
        except Exception as exc:
            last_error = str(exc)
            path.unlink(missing_ok=True)
        if attempt < attempts:
            time.sleep(2 * attempt)
    if expected_size is not None and total != expected_size:
        raise ValueError(last_error or f"incomplete download: got {total} bytes, expected {expected_size}")
    if not path.exists():
        raise ValueError(last_error or "download failed")
    return path, total, expected_size


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


def validate_input_output_dir(relative_path: str) -> Path:
    if not relative_path:
        raise ValueError("output_dir is required")
    lowered = relative_path.lower()
    if any(hint in lowered for hint in FORBIDDEN_PATH_HINTS):
        raise ValueError("output_dir appears to reference a forbidden private/system location")
    candidate = Path(relative_path)
    if candidate.is_absolute():
        raise ValueError("Only relative output_dir values under workspace/inputs are allowed")
    resolved = (INPUTS_DIR / candidate).resolve(strict=False)
    if not is_relative_to(resolved, INPUTS_DIR):
        raise ValueError("output_dir must stay inside mcp_server/workspace/inputs")
    return resolved


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
    log_file = log_event("list_available_pipelines", job_id, {"pipelines": PIPELINES, "synthetic_test_pipelines": SYNTHETIC_TEST_PIPELINES})
    return response_payload(
        status="ok",
        project_id=None,
        job_id=job_id,
        output_dir=OUTPUTS_DIR,
        log_file=log_file,
        message="Real-data whitelisted pipelines returned. Synthetic demo/selftest tools are hidden from the normal pipeline list and require explicit confirmation.",
        data={"pipelines": PIPELINES, "synthetic_test_pipelines": SYNTHETIC_TEST_PIPELINES},
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


def tool_run_demo_pipeline(project_id: str, confirm_synthetic: bool = False) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    project_dir = require_project_dir(project_id, create=True)
    if not confirm_synthetic:
        log_file = log_event("run_demo_pipeline", job_id, {"project_id": project_id, "status": "blocked_requires_confirmation"})
        return response_payload(
            status="error",
            project_id=safe_slug(project_id),
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message="Synthetic demo generation is disabled by default. For real analysis, download/upload data into workspace/inputs and call run_seurat_basic or run_scanpy_basic. To intentionally test demo output, call run_demo_pipeline with confirm_synthetic=true.",
            warnings=["This tool must not be used as evidence that real GEO data were analyzed."],
        )
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
    report_candidates = [
        project_dir / "manuscript_report.md",
        project_dir / "result_quality_check.md",
        project_dir / "downstream_proposal.md",
        project_dir / "report_skeleton.md",
    ]
    report = next((path for path in report_candidates if path.exists()), report_candidates[0])
    if not report.exists():
        log_file = log_event("read_report", job_id, {"project_id": project_id, "status": "missing_report"})
        return response_payload(
            status="error",
            project_id=safe_slug(project_id),
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message="No readable report found for this project. Run a real workflow and validate/report the result first.",
            warnings=["Expected one of manuscript_report.md, result_quality_check.md, downstream_proposal.md, or report_skeleton.md."],
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


def tool_list_geo_supplementary_files(gse_accession: str, file_regex: str = "") -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    files = fetch_geo_supplementary_listing(gse_accession)
    if file_regex:
        pattern = re.compile(file_regex, flags=re.IGNORECASE)
        files = [item for item in files if pattern.search(item["filename"])]
    log_file = log_event("list_geo_supplementary_files", job_id, {"gse_accession": gse_accession, "n_files": len(files)})
    return response_payload(
        status="ok",
        project_id=None,
        job_id=job_id,
        output_dir=INPUTS_DIR,
        log_file=log_file,
        message="GEO supplementary files listed. No data were generated or analyzed.",
        warnings=["Choose real files to download before running Scanpy or Seurat. Listing a GEO directory is not an analysis."],
        data={"gse_accession": gse_accession.strip().upper(), "suppl_url": geo_suppl_url(gse_accession), "files": files},
    )


def tool_download_geo_supplementary(
    project_id: str,
    gse_accession: str,
    file_regex: str = "",
    max_files: int = 20,
    overwrite: bool = False,
    max_bytes_per_file: int = MAX_DOWNLOAD_BYTES,
) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    project_slug = safe_slug(project_id)
    require_project_dir(project_slug, create=True)
    files = fetch_geo_supplementary_listing(gse_accession)
    if file_regex:
        pattern = re.compile(file_regex, flags=re.IGNORECASE)
        files = [item for item in files if pattern.search(item["filename"])]
    files = files[: max(1, min(int(max_files), 100))]
    if not files:
        raise FileNotFoundError("No GEO supplementary files matched the requested accession/filter.")

    downloaded: list[dict[str, Any]] = []
    errors: list[str] = []
    target_dir = f"inputs/{project_slug}/{gse_accession.strip().upper()}"
    for item in files:
        filename = safe_slug(item["filename"], "geo_file")
        relative_path = f"{target_dir}/{filename}"
        try:
            path, size, expected_size = download_url_to_workspace(item["url"], relative_path, overwrite=overwrite, max_bytes=max_bytes_per_file)
            downloaded.append({"filename": item["filename"], "workspace_path": workspace_relative(path), "url": item["url"], "size": size, "expected_size": expected_size})
        except Exception as exc:
            errors.append(f"{item['filename']}: {exc}")

    log_file = log_event("download_geo_supplementary", job_id, {"project_id": project_slug, "gse_accession": gse_accession, "downloaded": len(downloaded), "errors": errors})
    status = "ok" if downloaded and not errors else ("warning" if downloaded else "error")
    return response_payload(
        status=status,
        project_id=project_slug,
        job_id=job_id,
        output_dir=INPUTS_DIR / project_slug / gse_accession.strip().upper(),
        log_file=log_file,
        message="Downloaded real GEO supplementary files into workspace/inputs." if downloaded else "No GEO supplementary files were downloaded.",
        warnings=errors + ["Archives must be extracted/inspected before choosing input_type if downloaded files are tar/tar.gz/tgz."],
        data={"downloaded": downloaded, "errors": errors, "input_dir": f"{project_slug}/{gse_accession.strip().upper()}"},
    )


def tool_extract_workspace_archive(archive_path: str, output_dir: str = "", overwrite: bool = False, max_members: int = 20000) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    archive = validate_input_path(archive_path)
    if not archive.exists() or not archive.is_file():
        raise FileNotFoundError("archive_path must point to an existing file under workspace/inputs")
    if not output_dir:
        name = archive.name
        for suffix in [".tar.gz", ".tgz", ".tar", ".zip", ".gz"]:
            if name.lower().endswith(suffix):
                name = name[: -len(suffix)]
                break
        output_dir = str(Path(archive_path).with_name(name))
    out_dir = validate_input_output_dir(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    extracted: list[str] = []

    def target_for(member_name: str) -> Path:
        normalized = Path(member_name)
        if normalized.is_absolute() or ".." in normalized.parts:
            raise ValueError(f"Unsafe archive member path: {member_name}")
        target = (out_dir / normalized).resolve(strict=False)
        if not is_relative_to(target, out_dir):
            raise ValueError(f"Archive member escapes output_dir: {member_name}")
        return target

    suffixes = "".join(archive.suffixes).lower()
    if suffixes.endswith((".tar", ".tar.gz", ".tgz")):
        mode = "r:gz" if suffixes.endswith((".tar.gz", ".tgz")) else "r:"
        with tarfile.open(archive, mode) as tf:
            members = tf.getmembers()
            if len(members) > max_members:
                raise ValueError(f"Archive has {len(members)} members, exceeding max_members={max_members}")
            for member in members:
                target = target_for(member.name)
                if member.isdir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                if not member.isfile():
                    continue
                if target.exists() and not overwrite:
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                src = tf.extractfile(member)
                if src is None:
                    continue
                with src, target.open("wb") as handle:
                    shutil.copyfileobj(src, handle)
                extracted.append(workspace_relative(target))
    elif suffixes.endswith(".zip"):
        with zipfile.ZipFile(archive) as zf:
            infos = zf.infolist()
            if len(infos) > max_members:
                raise ValueError(f"Archive has {len(infos)} members, exceeding max_members={max_members}")
            for info in infos:
                target = target_for(info.filename)
                if info.is_dir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                if target.exists() and not overwrite:
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                with zf.open(info) as src, target.open("wb") as handle:
                    shutil.copyfileobj(src, handle)
                extracted.append(workspace_relative(target))
    elif suffixes.endswith(".gz"):
        target = out_dir / archive.with_suffix("").name
        if target.exists() and not overwrite:
            raise FileExistsError("Extracted file exists. Set overwrite=true to replace it.")
        target.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(archive, "rb") as src, target.open("wb") as handle:
            shutil.copyfileobj(src, handle)
        extracted.append(workspace_relative(target))
    else:
        raise ValueError("Supported archives: .tar, .tar.gz, .tgz, .zip, or single-file .gz")

    log_file = log_event("extract_workspace_archive", job_id, {"archive_path": archive_path, "output_dir": output_dir, "n_files": len(extracted)})
    return response_payload(
        status="ok",
        project_id=None,
        job_id=job_id,
        output_dir=out_dir,
        log_file=log_file,
        message="Archive extracted inside workspace/inputs. Inspect extracted files before choosing the analysis input_type.",
        warnings=[] if extracted else ["No files were extracted; existing files may have been skipped because overwrite=false."],
        data={"archive": workspace_relative(archive), "output_dir": str(out_dir.relative_to(INPUTS_DIR)).replace("\\", "/"), "files": extracted[:500], "n_files": len(extracted)},
    )


def validate_optional_input_path(input_path: str) -> str:
    if not input_path:
        return ""
    return str(validate_input_path(input_path))


def infer_input_type(path: Path, input_type: str = "") -> str:
    if input_type:
        return input_type
    if path.is_dir():
        if is_nonstandard_10x_dir(path):
            return "10x_nonstandard"
        return "10x_mtx"
    suffixes = "".join(path.suffixes).lower()
    if suffixes.endswith(".h5ad"):
        return "h5ad"
    if suffixes.endswith(".h5") or suffixes.endswith(".hdf5"):
        return "10x_h5"
    if suffixes.endswith(".rds"):
        return "rds"
    if suffixes.endswith(".csv") or suffixes.endswith(".csv.gz"):
        return "csv"
    return "unknown"


def is_nonstandard_10x_dir(path: Path) -> bool:
    def has_triplet(directory: Path) -> bool:
        names = {child.name.lower() for child in directory.iterdir() if child.is_file()}
        has_matrix = bool({"count_matrix_sparse.mtx", "count_matrix_sparse.mtx.gz"} & names)
        has_barcodes = bool({"count_matrix_barcodes.tsv", "count_matrix_barcodes.tsv.gz"} & names)
        has_genes = bool({"count_matrix_genes.tsv", "count_matrix_genes.tsv.gz"} & names)
        return has_matrix and has_barcodes and has_genes

    if not path.exists() or not path.is_dir():
        return False
    if has_triplet(path):
        return True
    try:
        return any(child.is_dir() and has_triplet(child) for child in path.iterdir())
    except OSError:
        return False


def find_rscript() -> str | None:
    candidates = [
        os.environ.get("RSCRIPT"),
        shutil.which("Rscript"),
        r"E:\R-4.4.2\bin\Rscript.exe",
        r"C:\Program Files\R\R-4.4.2\bin\Rscript.exe",
        r"C:\Program Files\R\R-4.4.1\bin\Rscript.exe",
        r"C:\Program Files\R\R-4.3.3\bin\Rscript.exe",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return str(candidate)
    return None


def find_analysis_python() -> str:
    candidates = [
        os.environ.get("SCMECHANISM_PYTHON"),
        os.environ.get("ANALYSIS_PYTHON"),
        DEFAULT_ANALYSIS_PYTHON,
        r"C:\ProgramData\miniconda3\envs\scmechanism-agent\python.exe",
        r"C:\ProgramData\miniconda3\envs\singlecell-skill\python.exe",
        sys.executable,
    ]
    for candidate in candidates:
        if candidate and Path(candidate).exists():
            return str(candidate)
    return sys.executable


def check_python_packages(python_exe: str, packages: list[str], timeout_seconds: int = 120) -> tuple[dict[str, bool], str]:
    code = (
        "import importlib.util, json, sys; "
        f"mods={packages!r}; "
        "print(json.dumps({'executable': sys.executable, 'version': sys.version, 'packages': {m: importlib.util.find_spec(m) is not None for m in mods}}, ensure_ascii=False))"
    )
    completed = subprocess.run(
        [python_exe, "-c", code],
        cwd=SERVER_ROOT,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        timeout=timeout_seconds,
        shell=False,
    )
    if completed.returncode != 0:
        return {pkg: False for pkg in packages}, (completed.stderr or "").strip() or (completed.stdout or "").strip()
    payload = json.loads((completed.stdout or "").strip().splitlines()[-1])
    return {str(k): bool(v) for k, v in payload["packages"].items()}, payload["version"]


def skill_runtime_path(relative_path: str, *, must_exist: bool = True, allow_dir: bool = False) -> Path:
    candidate = Path(relative_path)
    if candidate.is_absolute() or ".." in candidate.parts:
        raise ValueError("Skill runtime paths must be relative and cannot use parent traversal")
    resolved = (SKILL_RUNTIME_DIR / candidate).resolve(strict=False)
    if not is_relative_to(resolved, SKILL_RUNTIME_DIR):
        raise ValueError("Skill runtime path escaped mcp_server/skill_runtime")
    if must_exist and not resolved.exists():
        raise FileNotFoundError(f"Skill runtime file not found: {relative_path}")
    if resolved.exists() and resolved.is_dir() and not allow_dir:
        raise ValueError("Expected a skill runtime file, got a directory")
    return resolved


def render_template_file(template: Path, out_file: Path, replacements: dict[str, str]) -> Path:
    text = template.read_text(encoding="utf-8")
    for key, value in replacements.items():
        text = text.replace("{{" + key + "}}", value)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(text, encoding="utf-8")
    return out_file


def r_path(path: Path | str) -> str:
    if not path:
        return ""
    return str(path).replace("\\", "/")


def r_string(value: str) -> str:
    """Escape a Python string for insertion inside an existing quoted R string."""
    return str(value or "").replace("\\", "\\\\").replace('"', '\\"')


def run_fixed_process(
    *,
    tool: str,
    job_id: str,
    project_dir: Path,
    command: list[str],
    timeout_seconds: int,
    env_extra: dict[str, str] | None = None,
) -> tuple[subprocess.CompletedProcess[str], Path]:
    safe_timeout = max(1, min(int(timeout_seconds), 24 * 60 * 60))
    run_log = project_dir / "logs" / f"{tool}-{job_id}.log"
    run_log.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    if env_extra:
        env.update(env_extra)
    completed = subprocess.run(
        command,
        cwd=SERVER_ROOT,
        env=env,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        timeout=safe_timeout,
        shell=False,
    )
    stdout = completed.stdout or ""
    stderr = completed.stderr or ""
    run_log.write_text(
        json.dumps(
            {
                "time": utc_now(),
                "command": command,
                "returncode": completed.returncode,
                "stdout": stdout[-20000:],
                "stderr": stderr[-20000:],
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    return completed, run_log


def collect_project_files(project_dir: Path, limit: int = 200) -> list[str]:
    return [
        str(path.relative_to(project_dir)).replace("\\", "/")
        for path in sorted(project_dir.rglob("*"))
        if path.is_file()
    ][:limit]


def safe_tool_status(tool_name: str, callback: Callable[..., dict[str, Any]], **kwargs: Any) -> dict[str, Any]:
    try:
        result = callback(**kwargs)
        if isinstance(result, dict):
            return result
        return {"status": "error", "message": f"{tool_name} returned a non-dict result", "data": {"repr": repr(result)}}
    except Exception as exc:
        return {"status": "error", "message": f"{tool_name} failed: {exc}", "data": {}}


def write_selftest_10x_input() -> dict[str, Any]:
    out = INPUTS_DIR / "mcp_selftest_10x_three_groups"
    out.mkdir(parents=True, exist_ok=True)
    for path in out.glob("*"):
        if path.is_file():
            path.unlink()
    genes = [
        "CD3D", "CD3E", "TRAC", "NKG7", "MS4A1", "CD79A", "CD79B", "CD74",
        "LYZ", "LST1", "AIF1", "S100A8", "COL1A1", "COL1A2", "DCN", "LUM",
        "PECAM1", "VWF", "MKI67", "TOP2A", "MALAT1", "ACTB", "GAPDH", "RPLP0",
        "CXCL12", "CXCR4", "CCL2", "CCR2", "IL6", "IL6R", "IL6ST", "VEGFA", "KDR",
        "SPP1", "CD44", "FN1", "ITGA5", "ITGB1", "MMP9", "TGFB1", "TGFBR1", "TGFBR2",
    ]
    groups = ["T_cell"] * 80 + ["B_cell"] * 80 + ["Myeloid"] * 80
    boosts = {
        "T_cell": {"CD3D", "CD3E", "TRAC", "NKG7", "CXCR4", "IL6R", "CD44", "ITGB1"},
        "B_cell": {"MS4A1", "CD79A", "CD79B", "CD74", "CXCR4", "IL6R", "ITGB1"},
        "Myeloid": {"LYZ", "LST1", "AIF1", "S100A8", "CCL2", "MMP9", "IL6", "TGFB1", "SPP1"},
    }
    shared_lr = {"CXCL12", "VEGFA", "FN1", "ITGA5", "TGFBR1", "TGFBR2"}
    housekeeping = {"MALAT1", "ACTB", "GAPDH", "RPLP0"}
    entries: list[tuple[int, int, int]] = []
    for cell_index, group in enumerate(groups, start=1):
        for gene_index, gene in enumerate(genes, start=1):
            value = 1
            if gene in boosts[group]:
                value += 10 + (cell_index % 3)
            if gene in shared_lr:
                value += 3 + (cell_index % 2)
            if gene in housekeeping:
                value += 4
            if value > 0:
                entries.append((gene_index, cell_index, value))
    with gzip.open(out / "matrix.mtx.gz", "wt", encoding="utf-8") as handle:
        handle.write("%%MatrixMarket matrix coordinate integer general\n")
        handle.write("% generated by scMechanism local MCP selftest\n")
        handle.write(f"{len(genes)} {len(groups)} {len(entries)}\n")
        for row, col, value in entries:
            handle.write(f"{row} {col} {value}\n")
    with gzip.open(out / "features.tsv.gz", "wt", encoding="utf-8") as handle:
        for gene in genes:
            handle.write(f"{gene}\t{gene}\tGene Expression\n")
    with gzip.open(out / "barcodes.tsv.gz", "wt", encoding="utf-8") as handle:
        for i in range(len(groups)):
            handle.write(f"SELFTEST_{i + 1:03d}-1\n")
    metadata = out / "metadata.csv"
    with metadata.open("w", encoding="utf-8") as handle:
        handle.write("barcode,known_group,sample_id,batch,condition,stage\n")
        for i, group in enumerate(groups):
            handle.write(f"SELFTEST_{i + 1:03d}-1,{group},selftest_1,batch1,test,{group}\n")
    summary = {
        "input_path": "mcp_selftest_10x_three_groups",
        "metadata_path": "mcp_selftest_10x_three_groups/metadata.csv",
        "cells": len(groups),
        "genes": len(genes),
        "groups": sorted(set(groups)),
    }
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    return summary


def render_seurat_workflow_scripts(
    project_dir: Path,
    input_path: Path,
    *,
    input_type: str,
    species: str,
    sample_id: str,
    metadata_path: str = "",
    batch_col: str = "",
    condition_col: str = "",
    reference_name: str = "",
    celltype_col: str = "singleR_label",
    subset_query: str = "",
    root_query: str = "",
) -> dict[str, Path]:
    scripts_dir = project_dir / "scripts"
    seurat_out = project_dir / "seurat_basic"
    annotation_out = project_dir / "singler_annotation"
    marker_out = project_dir / "marker_enrichment"
    cellchat_out = project_dir / "cellchat"
    monocle_out = project_dir / "monocle3"
    organism = species.lower()
    reference = reference_name or ("mouse_rnaseq" if organism == "mouse" else "hpca")
    rds_core = seurat_out / "objects" / "processed_seurat.rds"
    rds_annotated = annotation_out / "objects" / "seurat_singleR_annotated.rds"

    scripts = {
        "seurat_basic": render_template_file(
            skill_runtime_path("scripts/course_adapted/01_seurat_v5_core_pipeline.R"),
            scripts_dir / "01_seurat_v5_core_pipeline.R",
            {
                "INPUT_PATH": r_path(input_path),
                "METADATA_PATH": r_path(metadata_path),
                "OUTPUT_DIR": r_path(seurat_out),
                "ORGANISM": organism,
                "SAMPLE_ID": sample_id,
                "BATCH_COL": batch_col,
                "CONDITION_COL": condition_col,
                "INPUT_TYPE": input_type,
            },
        ),
        "singler_annotation": render_template_file(
            skill_runtime_path("scripts/course_adapted/05_singler_cell_annotation.R"),
            scripts_dir / "05_singler_cell_annotation.R",
            {
                "INPUT_RDS": r_path(rds_core),
                "OUTPUT_DIR": r_path(annotation_out),
                "ORGANISM": organism,
                "REFERENCE_NAME": reference,
                "CLUSTER_COL": "seurat_clusters",
            },
        ),
        "marker_enrichment": render_template_file(
            skill_runtime_path("scripts/course_adapted/02_marker_enrichment_from_seurat.R"),
            scripts_dir / "02_marker_enrichment_from_seurat.R",
            {
                "SEURAT_RDS": r_path(rds_annotated),
                "OUTPUT_DIR": r_path(marker_out),
                "CLUSTER_COL": celltype_col,
                "ORGANISM": organism,
                "TOP_N": "20",
            },
        ),
        "cellchat": render_template_file(
            skill_runtime_path("scripts/course_adapted/03_cellchat_from_seurat.R"),
            scripts_dir / "03_cellchat_from_seurat.R",
            {
                "SEURAT_RDS": r_path(rds_annotated),
                "OUTPUT_DIR": r_path(cellchat_out),
                "CELLTYPE_COL": celltype_col,
                "ORGANISM": organism,
                "MIN_CELLS": "10",
            },
        ),
        "monocle3": render_template_file(
            skill_runtime_path("scripts/course_adapted/04_monocle3_from_seurat.R"),
            scripts_dir / "04_monocle3_from_seurat.R",
            {
                "SEURAT_RDS": r_path(rds_annotated),
                "OUTPUT_DIR": r_path(monocle_out),
                "CELLTYPE_COL": celltype_col,
                "SUBSET_QUERY": r_string(subset_query),
                "ROOT_CELLS_FILE": "",
                "ROOT_QUERY": r_string(root_query),
                "MAX_CELLS": "8000",
                "BALANCE_COL": celltype_col,
                "USE_SEURAT_UMAP": "true",
                "RANDOM_SEED": "20260702",
            },
        ),
    }
    (scripts_dir / "README.md").write_text(
        "\n".join(
            [
                "# Rendered scMechanism MCP Workflow Scripts",
                "",
                "Run order:",
                "1. 01_seurat_v5_core_pipeline.R",
                "2. 05_singler_cell_annotation.R",
                "3. 02_marker_enrichment_from_seurat.R",
                "4. Generate downstream_proposal.md and ask the user for approval.",
                "5. Only after approval, run 03_cellchat_from_seurat.R or 04_monocle3_from_seurat.R.",
                "",
                "All inputs and outputs are restricted to mcp_server/workspace.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return scripts


def tool_list_skill_runtime_files(relative_path: str = "", max_files: int = 300) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    root = SKILL_RUNTIME_DIR if not relative_path else skill_runtime_path(relative_path, must_exist=True, allow_dir=True)
    files = []
    if root.is_file():
        paths = [root]
    else:
        paths = [path for path in sorted(root.rglob("*")) if path.is_file()]
    for path in paths[: max(1, min(int(max_files), 1000))]:
        files.append({
            "path": str(path.relative_to(SKILL_RUNTIME_DIR)).replace("\\", "/"),
            "size": path.stat().st_size,
        })
    log_file = log_event("list_skill_runtime_files", job_id, {"relative_path": relative_path, "n_files": len(files)})
    return response_payload(
        status="ok",
        job_id=job_id,
        output_dir=SKILL_RUNTIME_DIR,
        log_file=log_file,
        message="Read-only skill runtime files listed.",
        data={"files": files, "skill_runtime": str(SKILL_RUNTIME_DIR)},
    )


def tool_read_skill_runtime_file(relative_path: str, max_bytes: int = MAX_READ_BYTES) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    path = skill_runtime_path(relative_path, must_exist=True, allow_dir=False)
    max_bytes = max(1, min(int(max_bytes), MAX_READ_BYTES))
    raw = path.read_bytes()
    truncated = len(raw) > max_bytes
    raw_slice = raw[:max_bytes]
    data: dict[str, Any] = {
        "path": str(path.relative_to(SKILL_RUNTIME_DIR)).replace("\\", "/"),
        "size": len(raw),
        "truncated": truncated,
    }
    if path.suffix.lower() in TEXT_FILE_SUFFIXES:
        data["content"] = raw_slice.decode("utf-8", errors="replace")
        data["encoding"] = "utf-8"
    else:
        data["content_base64"] = base64.b64encode(raw_slice).decode("ascii")
        data["encoding"] = "base64"
    log_file = log_event("read_skill_runtime_file", job_id, {"relative_path": relative_path, "bytes": len(raw_slice)})
    return response_payload(
        status="ok",
        job_id=job_id,
        output_dir=SKILL_RUNTIME_DIR,
        log_file=log_file,
        message="Read-only skill runtime file returned.",
        warnings=["File was truncated to max_bytes."] if truncated else [],
        data=data,
    )


def tool_check_runtime_environment(project_id: str = "environment_check") -> dict[str, Any]:
    project_dir = require_project_dir(project_id, create=True)
    job_id = uuid.uuid4().hex
    packages = ["scanpy", "anndata", "pandas", "numpy", "scipy", "matplotlib", "seaborn", "celltypist", "gseapy", "bbknn", "scanorama", "harmonypy", "liana", "leidenalg", "igraph"]
    analysis_python = find_analysis_python()
    python_status, python_detail = check_python_packages(analysis_python, packages)
    rscript = find_rscript()
    r_status: dict[str, Any] = {"Rscript": rscript, "packages": {}}
    if rscript:
        check_expr = (
            "pkgs <- c('Seurat','SingleR','celldex','clusterProfiler','CellChat','monocle3','harmony'); "
            "cat(paste(pkgs, vapply(pkgs, requireNamespace, logical(1), quietly=TRUE), sep='=', collapse='\\n'))"
        )
        completed, run_log = run_fixed_process(
            tool="check_runtime_environment_r",
            job_id=job_id,
            project_dir=project_dir,
            command=[rscript, "-e", check_expr],
            timeout_seconds=120,
        )
        r_status["returncode"] = completed.returncode
        for line in (completed.stdout or "").splitlines():
            if "=" in line:
                key, value = line.split("=", 1)
                r_status["packages"][key] = value.strip().lower() == "true"
    else:
        run_log = project_dir / "logs" / f"check_runtime_environment-{job_id}.log"
        run_log.parent.mkdir(parents=True, exist_ok=True)
        run_log.write_text("Rscript not found\n", encoding="utf-8")
    report = {
        "time_utc": utc_now(),
        "server_python": sys.version,
        "analysis_python": analysis_python,
        "analysis_python_detail": python_detail,
        "python_packages": python_status,
        "r": r_status,
        "skill_runtime": str(SKILL_RUNTIME_DIR),
    }
    out_json = project_dir / "environment_status.json"
    out_json.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    warnings = []
    if not all(python_status.values()):
        warnings.append("Some Python analysis packages are missing in the analysis Python environment; scanpy_basic optional modules may be limited.")
    if not rscript:
        warnings.append("Rscript is not available. Add R >= 4.3 Rscript to PATH or set RSCRIPT.")
    elif not all(r_status.get("packages", {}).values()):
        warnings.append("Some R packages required by Seurat/SingleR/CellChat/Monocle3 are missing.")
    log_file = log_event("check_runtime_environment", job_id, {"project_id": project_id, "warnings": warnings})
    return response_payload(
        status="ok" if not warnings else "warning",
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Runtime environment checked.",
        warnings=warnings,
        data={"environment_status": str(out_json), **report},
    )


def tool_render_workflow_scripts(
    project_id: str,
    input_path: str,
    species: str = "human",
    input_type: str = "",
    metadata_path: str = "",
    sample_id: str = "",
    batch_col: str = "",
    condition_col: str = "",
    reference_name: str = "",
) -> dict[str, Any]:
    project_dir = require_project_dir(project_id, create=True)
    job_id = uuid.uuid4().hex
    input_candidate = validate_input_path(input_path)
    metadata_candidate = validate_optional_input_path(metadata_path)
    inferred_type = infer_input_type(input_candidate, input_type)
    if inferred_type not in {"10x_mtx", "10x_nonstandard", "10x_h5", "rds", "csv"}:
        raise ValueError("Seurat workflow supports 10x_mtx, 10x_nonstandard, 10x_h5, rds, and csv inputs.")
    scripts = render_seurat_workflow_scripts(
        project_dir,
        input_candidate,
        input_type=inferred_type,
        species=species,
        sample_id=sample_id or safe_slug(project_id),
        metadata_path=metadata_candidate,
        batch_col=batch_col,
        condition_col=condition_col,
        reference_name=reference_name,
    )
    log_file = log_event("render_workflow_scripts", job_id, {"project_id": project_id, "scripts": list(scripts)})
    return response_payload(
        status="ok",
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Workflow scripts rendered into the project folder. Downstream CellChat/Monocle3 scripts are approval-gated.",
        data={"scripts": {key: str(value) for key, value in scripts.items()}},
    )


def tool_run_scanpy_basic(
    project_id: str,
    input_path: str,
    species: str = "human",
    input_type: str = "",
    batch_col: str = "",
    metadata_path: str = "",
    sample_id: str = "",
    group_col: str = "",
    annotation_method: str = "marker_summary",
    timeout_seconds: int = 3600,
) -> dict[str, Any]:
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
        raise FileNotFoundError("Input path does not exist under workspace/inputs.")
    if input_type:
        inferred_type = input_type
    else:
        inferred_type = infer_input_type(candidate)
    if inferred_type not in {"h5ad", "10x_mtx", "10x_h5"}:
        raise ValueError("Scanpy workflow supports h5ad, 10x_mtx, and 10x_h5 inputs.")
    metadata_candidate = validate_optional_input_path(metadata_path)
    analysis_python = find_analysis_python()
    py_status, py_detail = check_python_packages(analysis_python, ["scanpy", "anndata", "pandas", "numpy", "scipy", "matplotlib", "seaborn", "leidenalg", "igraph"])
    if not py_status.get("scanpy"):
        log_file = log_event("run_scanpy_basic", job_id, {"project_id": project_id, "status": "missing_scanpy"})
        return response_payload(
            status="error",
            project_id=safe_slug(project_id),
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message="Python package scanpy is not installed in the analysis Python environment. Create or repair the environment before running scanpy_basic.",
            warnings=["Required environment: Python >= 3.10 with scanpy/anndata/pandas/numpy/scipy/matplotlib/seaborn/leidenalg/igraph.", f"analysis_python={analysis_python}", f"detail={py_detail}"],
        )
    output_dir = project_dir / "scanpy_basic"
    completed, run_log = run_fixed_process(
        tool="run_scanpy_basic",
        job_id=job_id,
        project_dir=project_dir,
        command=[
            analysis_python,
            str(SERVER_ROOT / "pipelines" / "run_scanpy_basic.py"),
            "--input-path",
            str(candidate),
            "--input-type",
            inferred_type,
            "--output-dir",
            str(output_dir),
            "--species",
            species,
            "--batch-col",
            batch_col,
            "--metadata-path",
            metadata_candidate,
            "--sample-id",
            sample_id or safe_slug(project_id),
            "--group-col",
            group_col,
            "--annotation-method",
            annotation_method,
        ],
        timeout_seconds=timeout_seconds,
        env_extra={"SCMECHANISM_WORKSPACE": str(WORKSPACE_ROOT), "SCMECHANISM_ANALYSIS_PYTHON": analysis_python},
    )
    status = "ok" if completed.returncode == 0 else "error"
    log_file = log_event("run_scanpy_basic", job_id, {"project_id": project_id, "status": status, "returncode": completed.returncode})
    return response_payload(
        status=status,
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Scanpy basic workflow completed." if status == "ok" else "Scanpy basic workflow failed. Read the run log for details.",
        warnings=warnings if status == "ok" else warnings + ["See run_log for stdout/stderr."],
        data={"analysis_python": analysis_python, "run_log": str(run_log), "files": collect_project_files(output_dir)},
    )


def tool_run_seurat_basic(
    project_id: str,
    input_path: str,
    species: str = "human",
    input_type: str = "",
    metadata_path: str = "",
    sample_id: str = "",
    batch_col: str = "",
    condition_col: str = "",
    reference_name: str = "",
    run_annotation: bool = True,
    run_marker_enrichment: bool = True,
    timeout_seconds: int = 7200,
) -> dict[str, Any]:
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
        raise FileNotFoundError("Input path does not exist under workspace/inputs.")
    rscript = find_rscript()
    if not rscript:
        log_file = log_event("run_seurat_basic", job_id, {"project_id": project_id, "status": "missing_rscript"})
        return response_payload(
            status="error",
            project_id=safe_slug(project_id),
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message="Rscript is not available. Install R >= 4.3 and add Rscript to PATH, or set RSCRIPT to Rscript.exe.",
            warnings=["Use the main skill environment scripts before running Seurat modules."],
        )
    metadata_candidate = validate_optional_input_path(metadata_path)
    inferred_type = infer_input_type(candidate, input_type)
    if inferred_type not in {"10x_mtx", "10x_nonstandard", "10x_h5", "rds", "csv"}:
        raise ValueError("Seurat workflow supports 10x_mtx, 10x_nonstandard, 10x_h5, rds, and csv inputs.")
    scripts = render_seurat_workflow_scripts(
        project_dir,
        candidate,
        input_type=inferred_type,
        species=species,
        sample_id=sample_id or safe_slug(project_id),
        metadata_path=metadata_candidate,
        batch_col=batch_col,
        condition_col=condition_col,
        reference_name=reference_name,
    )
    run_sequence = [("seurat_basic", scripts["seurat_basic"])]
    if run_annotation:
        run_sequence.append(("singler_annotation", scripts["singler_annotation"]))
    if run_marker_enrichment:
        run_sequence.append(("marker_enrichment", scripts["marker_enrichment"]))
    run_logs: list[str] = []
    for step_name, script in run_sequence:
        completed, run_log = run_fixed_process(
            tool=step_name,
            job_id=job_id,
            project_dir=project_dir,
            command=[rscript, str(script)],
            timeout_seconds=timeout_seconds,
        )
        run_logs.append(str(run_log))
        if completed.returncode != 0:
            log_file = log_event("run_seurat_basic", job_id, {"project_id": project_id, "status": "failed", "step": step_name, "returncode": completed.returncode})
            return response_payload(
                status="error",
                project_id=safe_slug(project_id),
                job_id=job_id,
                output_dir=project_dir,
                log_file=log_file,
                message=f"Seurat workflow failed at step: {step_name}.",
                warnings=["Read the step run log for stdout/stderr.", *warnings],
                data={"failed_step": step_name, "run_logs": run_logs, "files": collect_project_files(project_dir)},
            )
    proposal = safe_tool_status("propose_downstream_modules", tool_propose_downstream_modules, project_id=safe_slug(project_id))
    quality = safe_tool_status("validate_result_bundle", tool_validate_result_bundle, project_id=safe_slug(project_id))
    log_file = log_event("run_seurat_basic", job_id, {"project_id": project_id, "status": "ok", "steps": [name for name, _ in run_sequence]})
    return response_payload(
        status="ok",
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Seurat workflow completed through annotation/marker steps. Review downstream_proposal.md before approving CellChat or Monocle3.",
        warnings=warnings + ["CellChat and Monocle3 are not run automatically. Use approval-gated tools after reviewing the proposal."],
        data={
            "run_logs": run_logs,
            "downstream_proposal_status": proposal.get("status"),
            "quality_status": quality.get("status"),
            "files": collect_project_files(project_dir),
        },
    )


def tool_propose_downstream_modules(project_id: str) -> dict[str, Any]:
    project_dir = require_project_dir(project_id, create=False)
    job_id = uuid.uuid4().hex
    out_md = project_dir / "downstream_proposal.md"
    completed, run_log = run_fixed_process(
        tool="propose_downstream_modules",
        job_id=job_id,
        project_dir=project_dir,
        command=[
            sys.executable,
            str(skill_runtime_path("scripts/propose_downstream_modules.py")),
            "--result-dir",
            str(project_dir),
            "--out-md",
            str(out_md),
        ],
        timeout_seconds=300,
    )
    status = "ok" if completed.returncode == 0 else "error"
    log_file = log_event("propose_downstream_modules", job_id, {"project_id": project_id, "status": status})
    return response_payload(
        status=status,
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Downstream CellChat/Monocle3 proposal generated. User approval is required before running downstream modules." if status == "ok" else "Failed to generate downstream proposal.",
        warnings=[] if status == "ok" else ["See run_log for stdout/stderr."],
        data={"proposal_file": str(out_md), "run_log": str(run_log), "proposal": out_md.read_text(encoding="utf-8", errors="replace") if out_md.exists() else ""},
    )


def tool_validate_result_bundle(project_id: str) -> dict[str, Any]:
    project_dir = require_project_dir(project_id, create=False)
    job_id = uuid.uuid4().hex
    out_md = project_dir / "result_quality_check.md"
    completed, run_log = run_fixed_process(
        tool="validate_result_bundle",
        job_id=job_id,
        project_dir=project_dir,
        command=[
            sys.executable,
            str(skill_runtime_path("scripts/validate_result_bundle.py")),
            "--result-dir",
            str(project_dir),
            "--out-md",
            str(out_md),
        ],
        timeout_seconds=300,
    )
    status = "ok" if completed.returncode == 0 else "error"
    log_file = log_event("validate_result_bundle", job_id, {"project_id": project_id, "status": status})
    return response_payload(
        status=status,
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Result bundle validation report generated." if status == "ok" else "Result bundle validation failed.",
        warnings=[] if status == "ok" else ["See run_log for stdout/stderr."],
        data={"validation_file": str(out_md), "run_log": str(run_log), "validation": out_md.read_text(encoding="utf-8", errors="replace") if out_md.exists() else ""},
    )


def tool_run_cellchat(
    project_id: str,
    approval_token: str,
    celltype_col: str = "singleR_label",
    species: str = "human",
    min_cells: int = 10,
    timeout_seconds: int = 7200,
) -> dict[str, Any]:
    project_dir = require_project_dir(project_id, create=False)
    job_id = uuid.uuid4().hex
    if approval_token != "APPROVED_CELLCHAT":
        log_file = log_event("run_cellchat", job_id, {"project_id": project_id, "status": "approval_required"})
        return response_payload(
            status="approval_required",
            project_id=safe_slug(project_id),
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message="CellChat requires prior downstream proposal review and explicit user approval.",
            warnings=["Call propose_downstream_modules first, then call run_cellchat with approval_token='APPROVED_CELLCHAT'."],
        )
    rscript = find_rscript()
    if not rscript:
        raise FileNotFoundError("Rscript is not available. Add Rscript to PATH or set RSCRIPT.")
    seurat_rds = project_dir / "singler_annotation" / "objects" / "seurat_singleR_annotated.rds"
    if not seurat_rds.exists():
        seurat_rds = project_dir / "seurat_basic" / "objects" / "processed_seurat.rds"
    if not seurat_rds.exists():
        raise FileNotFoundError("No processed Seurat RDS found. Run run_seurat_basic first.")
    script = render_template_file(
        skill_runtime_path("scripts/course_adapted/03_cellchat_from_seurat.R"),
        project_dir / "scripts" / "03_cellchat_from_seurat.R",
        {
            "SEURAT_RDS": r_path(seurat_rds),
            "OUTPUT_DIR": r_path(project_dir / "cellchat"),
            "CELLTYPE_COL": celltype_col,
            "ORGANISM": species.lower(),
            "MIN_CELLS": str(min_cells),
        },
    )
    completed, run_log = run_fixed_process(
        tool="run_cellchat",
        job_id=job_id,
        project_dir=project_dir,
        command=[rscript, str(script)],
        timeout_seconds=timeout_seconds,
    )
    status = "ok" if completed.returncode == 0 else "error"
    log_file = log_event("run_cellchat", job_id, {"project_id": project_id, "status": status, "returncode": completed.returncode})
    return response_payload(
        status=status,
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="CellChat completed." if status == "ok" else "CellChat failed. Read the run log.",
        warnings=[] if status == "ok" else ["See run_log for stdout/stderr."],
        data={"run_log": str(run_log), "files": collect_project_files(project_dir / "cellchat")},
    )


def tool_run_monocle3(
    project_id: str,
    approval_token: str,
    celltype_col: str = "singleR_label",
    subset_query: str = "",
    root_query: str = "",
    timeout_seconds: int = 7200,
) -> dict[str, Any]:
    project_dir = require_project_dir(project_id, create=False)
    job_id = uuid.uuid4().hex
    if approval_token != "APPROVED_MONOCLE3":
        log_file = log_event("run_monocle3", job_id, {"project_id": project_id, "status": "approval_required"})
        return response_payload(
            status="approval_required",
            project_id=safe_slug(project_id),
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message="Monocle3 requires prior downstream proposal review and explicit user approval.",
            warnings=["Call propose_downstream_modules first, then call run_monocle3 with approval_token='APPROVED_MONOCLE3'."],
        )
    rscript = find_rscript()
    if not rscript:
        raise FileNotFoundError("Rscript is not available. Add Rscript to PATH or set RSCRIPT.")
    seurat_rds = project_dir / "singler_annotation" / "objects" / "seurat_singleR_annotated.rds"
    if not seurat_rds.exists():
        raise FileNotFoundError("No SingleR-annotated Seurat RDS found. Run run_seurat_basic with annotation first.")
    script = render_template_file(
        skill_runtime_path("scripts/course_adapted/04_monocle3_from_seurat.R"),
        project_dir / "scripts" / "04_monocle3_from_seurat.R",
        {
            "SEURAT_RDS": r_path(seurat_rds),
            "OUTPUT_DIR": r_path(project_dir / "monocle3"),
            "CELLTYPE_COL": celltype_col,
            "SUBSET_QUERY": r_string(subset_query),
            "ROOT_CELLS_FILE": "",
            "ROOT_QUERY": r_string(root_query),
            "MAX_CELLS": "8000",
            "BALANCE_COL": celltype_col,
            "USE_SEURAT_UMAP": "true",
            "RANDOM_SEED": "20260702",
        },
    )
    completed, run_log = run_fixed_process(
        tool="run_monocle3",
        job_id=job_id,
        project_dir=project_dir,
        command=[rscript, str(script)],
        timeout_seconds=timeout_seconds,
    )
    status = "ok" if completed.returncode == 0 else "error"
    log_file = log_event("run_monocle3", job_id, {"project_id": project_id, "status": status, "returncode": completed.returncode})
    return response_payload(
        status=status,
        project_id=safe_slug(project_id),
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Monocle3 completed." if status == "ok" else "Monocle3 failed. Read the run log.",
        warnings=[] if status == "ok" else ["See run_log for stdout/stderr."],
        data={"run_log": str(run_log), "files": collect_project_files(project_dir / "monocle3")},
    )


def tool_run_full_workflow_selftest(project_id: str = "mcp_selftest_full_modules", timeout_seconds: int = 7200, confirm_synthetic: bool = False) -> dict[str, Any]:
    job_id = uuid.uuid4().hex
    project_slug = safe_slug(project_id, "mcp_selftest_full_modules")
    project_dir = require_project_dir(project_slug, create=True)
    if not confirm_synthetic:
        log_file = log_event("run_full_workflow_selftest", job_id, {"project_id": project_slug, "status": "blocked_requires_confirmation"})
        return response_payload(
            status="error",
            project_id=project_slug,
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message="Synthetic full workflow selftest is disabled by default. For real testing, call list_geo_supplementary_files/download_geo_supplementary or download_to_workspace, then run_seurat_basic/run_scanpy_basic on the downloaded input. To intentionally run the synthetic selftest, pass confirm_synthetic=true.",
            warnings=["This selftest creates synthetic data and must not be treated as real GEO analysis."],
        )
    input_summary = write_selftest_10x_input()
    steps: list[dict[str, Any]] = []

    seurat = tool_run_seurat_basic(
        project_id=project_slug,
        input_path=input_summary["input_path"],
        species="human",
        input_type="10x_mtx",
        metadata_path=input_summary["metadata_path"],
        sample_id="mcp_selftest",
        batch_col="batch",
        condition_col="known_group",
        reference_name="hpca",
        run_annotation=True,
        run_marker_enrichment=True,
        timeout_seconds=timeout_seconds,
    )
    steps.append({"step": "seurat_singler_marker", "status": seurat.get("status"), "message": seurat.get("message")})
    if seurat.get("status") != "ok":
        log_file = log_event("run_full_workflow_selftest", job_id, {"project_id": project_slug, "status": "failed", "failed_step": "seurat_singler_marker"})
        return response_payload(
            status="error",
            project_id=project_slug,
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message="Full workflow selftest failed at Seurat/SingleR/marker step.",
            warnings=["Read nested run logs in the project output directory."],
            data={"input": input_summary, "steps": steps, "files": collect_project_files(project_dir)},
        )

    cellchat = tool_run_cellchat(
        project_id=project_slug,
        approval_token="APPROVED_CELLCHAT",
        celltype_col="known_group",
        species="human",
        min_cells=10,
        timeout_seconds=timeout_seconds,
    )
    steps.append({"step": "cellchat", "status": cellchat.get("status"), "message": cellchat.get("message")})
    if cellchat.get("status") != "ok":
        log_file = log_event("run_full_workflow_selftest", job_id, {"project_id": project_slug, "status": "failed", "failed_step": "cellchat"})
        return response_payload(
            status="error",
            project_id=project_slug,
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message="Full workflow selftest failed at CellChat step.",
            warnings=["Read nested run logs in the project output directory."],
            data={"input": input_summary, "steps": steps, "files": collect_project_files(project_dir)},
        )

    monocle3 = tool_run_monocle3(
        project_id=project_slug,
        approval_token="APPROVED_MONOCLE3",
        celltype_col="known_group",
        subset_query='known_group %in% c("T_cell","B_cell","Myeloid")',
        root_query='known_group == "T_cell"',
        timeout_seconds=timeout_seconds,
    )
    steps.append({"step": "monocle3", "status": monocle3.get("status"), "message": monocle3.get("message")})
    if monocle3.get("status") != "ok":
        log_file = log_event("run_full_workflow_selftest", job_id, {"project_id": project_slug, "status": "failed", "failed_step": "monocle3"})
        return response_payload(
            status="error",
            project_id=project_slug,
            job_id=job_id,
            output_dir=project_dir,
            log_file=log_file,
            message="Full workflow selftest failed at Monocle3 step.",
            warnings=["Read nested run logs in the project output directory."],
            data={"input": input_summary, "steps": steps, "files": collect_project_files(project_dir)},
        )

    quality = tool_validate_result_bundle(project_id=project_slug)
    steps.append({"step": "result_quality", "status": quality.get("status"), "message": quality.get("message")})

    report = project_dir / "full_workflow_selftest_report.md"
    report.write_text(
        "\n".join(
            [
                "# Local MCP Full Workflow Selftest",
                "",
                "## Input",
                f"- Cells: {input_summary['cells']}",
                f"- Genes: {input_summary['genes']}",
                f"- Groups: {', '.join(input_summary['groups'])}",
                "",
                "## Step Status",
                *[f"- {step['step']}: {step['status']} - {step['message']}" for step in steps],
                "",
                "## Key Outputs",
                "- seurat_basic/objects/processed_seurat.rds",
                "- singler_annotation/objects/seurat_singleR_annotated.rds",
                "- marker_enrichment/tables/cluster_markers_top.csv",
                "- marker_enrichment/tables/GO_BP_top_markers.csv",
                "- cellchat/tables/cellchat_ligand_receptor.csv",
                "- monocle3/tables/pseudotime.csv",
                "- result_quality_check.md",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    log_file = log_event("run_full_workflow_selftest", job_id, {"project_id": project_slug, "status": "ok"})
    return response_payload(
        status="ok",
        project_id=project_slug,
        job_id=job_id,
        output_dir=project_dir,
        log_file=log_file,
        message="Full local MCP workflow selftest completed.",
        warnings=["This is a synthetic selftest dataset and is not biological evidence."],
        data={"input": input_summary, "steps": steps, "report": str(report), "files": collect_project_files(project_dir)},
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
    path, total, expected_size = download_url_to_workspace(url, relative_path, overwrite=overwrite, max_bytes=max_bytes)
    log_file = log_event("download_to_workspace", job_id, {"url": safe_url, "relative_path": relative_path, "bytes": total, "expected_size": expected_size})
    return response_payload(
        status="ok",
        job_id=job_id,
        output_dir=path.parent,
        log_file=log_file,
        message="Downloaded file into workspace.",
        data={"path": workspace_relative(path), "size": total, "expected_size": expected_size},
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
    analysis_python = find_analysis_python()
    completed = subprocess.run(
        [analysis_python, str(script), *safe_args],
        cwd=WORKSPACE_ROOT,
        env=env,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        timeout=timeout_seconds,
        shell=False,
    )
    stdout = completed.stdout or ""
    stderr = completed.stderr or ""
    run_log = project_dir / "logs" / f"run_workspace_python-{job_id}.log"
    run_log.parent.mkdir(parents=True, exist_ok=True)
    run_log.write_text(
        json.dumps(
            {
                "time": utc_now(),
                "script": workspace_relative(script),
                "python": analysis_python,
                "returncode": completed.returncode,
                "stdout": stdout[-10000:],
                "stderr": stderr[-10000:],
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
            "python": analysis_python,
            "returncode": completed.returncode,
            "stdout": stdout[-10000:],
            "stderr": stderr[-10000:],
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
    "list_geo_supplementary_files": tool_list_geo_supplementary_files,
    "download_geo_supplementary": tool_download_geo_supplementary,
    "extract_workspace_archive": tool_extract_workspace_archive,
    "list_skill_runtime_files": tool_list_skill_runtime_files,
    "read_skill_runtime_file": tool_read_skill_runtime_file,
    "check_runtime_environment": tool_check_runtime_environment,
    "render_workflow_scripts": tool_render_workflow_scripts,
    "run_scanpy_basic": tool_run_scanpy_basic,
    "run_seurat_basic": tool_run_seurat_basic,
    "propose_downstream_modules": tool_propose_downstream_modules,
    "validate_result_bundle": tool_validate_result_bundle,
    "run_cellchat": tool_run_cellchat,
    "run_monocle3": tool_run_monocle3,
    "run_full_workflow_selftest": tool_run_full_workflow_selftest,
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
    {"name": "list_result_files", "description": "List files under one project output directory.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}}, "required": ["project_id"]}},
    {"name": "read_report", "description": "Read report_skeleton.md from one project.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}}, "required": ["project_id"]}},
    {"name": "list_geo_supplementary_files", "description": "List real GEO supplementary files for a GSE accession without downloading or generating data.", "inputSchema": {"type": "object", "properties": {"gse_accession": {"type": "string"}, "file_regex": {"type": "string", "default": ""}}, "required": ["gse_accession"]}},
    {"name": "download_geo_supplementary", "description": "Download selected real GEO supplementary files into workspace/inputs/<project>/<GSE>; does not generate synthetic data or run analysis.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}, "gse_accession": {"type": "string"}, "file_regex": {"type": "string", "default": ""}, "max_files": {"type": "integer", "default": 20}, "overwrite": {"type": "boolean", "default": False}, "max_bytes_per_file": {"type": "integer", "default": MAX_DOWNLOAD_BYTES}}, "required": ["project_id", "gse_accession"]}},
    {"name": "extract_workspace_archive", "description": "Extract a downloaded archive under workspace/inputs into workspace/inputs, with path traversal protection.", "inputSchema": {"type": "object", "properties": {"archive_path": {"type": "string"}, "output_dir": {"type": "string", "default": ""}, "overwrite": {"type": "boolean", "default": False}, "max_members": {"type": "integer", "default": 20000}}, "required": ["archive_path"]}},
    {"name": "list_skill_runtime_files", "description": "List read-only copied Skill assets: agents, references, workflows, templates, and course scripts.", "inputSchema": {"type": "object", "properties": {"relative_path": {"type": "string", "default": ""}, "max_files": {"type": "integer", "default": 300}}}},
    {"name": "read_skill_runtime_file", "description": "Read a read-only copied Skill runtime file.", "inputSchema": {"type": "object", "properties": {"relative_path": {"type": "string"}, "max_bytes": {"type": "integer", "default": MAX_READ_BYTES}}, "required": ["relative_path"]}},
    {"name": "check_runtime_environment", "description": "Check Python and R package readiness for Scanpy, Seurat, SingleR, CellChat, and Monocle3.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string", "default": "environment_check"}}}},
    {"name": "render_workflow_scripts", "description": "Render Seurat, SingleR, marker, CellChat, and Monocle3 scripts for a real input under workspace/inputs without running them.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}, "input_path": {"type": "string"}, "species": {"type": "string", "default": "human"}, "input_type": {"type": "string", "default": ""}, "metadata_path": {"type": "string", "default": ""}, "sample_id": {"type": "string", "default": ""}, "batch_col": {"type": "string", "default": ""}, "condition_col": {"type": "string", "default": ""}, "reference_name": {"type": "string", "default": ""}}, "required": ["project_id", "input_path"]}},
    {"name": "run_scanpy_basic", "description": "Run the whitelisted Scanpy QC/clustering/annotation-evidence workflow on a real h5ad/10x input under workspace/inputs.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}, "input_path": {"type": "string"}, "species": {"type": "string", "default": "human"}, "input_type": {"type": "string", "default": ""}, "batch_col": {"type": "string", "default": ""}, "metadata_path": {"type": "string", "default": ""}, "sample_id": {"type": "string", "default": ""}, "group_col": {"type": "string", "default": ""}, "annotation_method": {"type": "string", "default": "marker_summary"}, "timeout_seconds": {"type": "integer", "default": 3600}}, "required": ["project_id", "input_path"]}},
    {"name": "run_seurat_basic", "description": "Run whitelisted Seurat V5 workflow on a real input under workspace/inputs; supports 10x_mtx, 10x_nonstandard, 10x_h5, rds, and csv.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}, "input_path": {"type": "string"}, "species": {"type": "string", "default": "human"}, "input_type": {"type": "string", "default": ""}, "metadata_path": {"type": "string", "default": ""}, "sample_id": {"type": "string", "default": ""}, "batch_col": {"type": "string", "default": ""}, "condition_col": {"type": "string", "default": ""}, "reference_name": {"type": "string", "default": ""}, "run_annotation": {"type": "boolean", "default": True}, "run_marker_enrichment": {"type": "boolean", "default": True}, "timeout_seconds": {"type": "integer", "default": 7200}}, "required": ["project_id", "input_path"]}},
    {"name": "propose_downstream_modules", "description": "Generate downstream_proposal.md from upstream annotation and marker results before CellChat/Monocle3.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}}, "required": ["project_id"]}},
    {"name": "validate_result_bundle", "description": "Generate result_quality_check.md from project outputs.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}}, "required": ["project_id"]}},
    {"name": "run_cellchat", "description": "Run CellChat only after explicit approval_token='APPROVED_CELLCHAT'.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}, "approval_token": {"type": "string"}, "celltype_col": {"type": "string", "default": "singleR_label"}, "species": {"type": "string", "default": "human"}, "min_cells": {"type": "integer", "default": 10}, "timeout_seconds": {"type": "integer", "default": 7200}}, "required": ["project_id", "approval_token"]}},
    {"name": "run_monocle3", "description": "Run Monocle3 only after explicit approval_token='APPROVED_MONOCLE3'.", "inputSchema": {"type": "object", "properties": {"project_id": {"type": "string"}, "approval_token": {"type": "string"}, "celltype_col": {"type": "string", "default": "singleR_label"}, "subset_query": {"type": "string", "default": ""}, "root_query": {"type": "string", "default": ""}, "timeout_seconds": {"type": "integer", "default": 7200}}, "required": ["project_id", "approval_token"]}},
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

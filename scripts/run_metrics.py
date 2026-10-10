"""Run a child process while recording elapsed time, CPU time, and peak RSS."""
from __future__ import annotations

import subprocess
import time
from pathlib import Path
import hashlib
import importlib.metadata
import os
import platform
import sys


def _process_stats(pid: int) -> tuple[float | None, int | None]:
    """Return cumulative CPU seconds and resident bytes from `ps` where available."""
    try:
        result = subprocess.run(
            ["ps", "-o", "time=,rss=", "-p", str(pid)],
            check=False, capture_output=True, text=True, timeout=2,
        )
        fields = result.stdout.split()
        if result.returncode or len(fields) < 2:
            return None, None
        cpu_text, rss_text = fields[0], fields[1]
        parts = cpu_text.split(":")
        if len(parts) == 3:
            cpu = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        elif len(parts) == 2:
            cpu = int(parts[0]) * 60 + float(parts[1])
        else:
            cpu = float(cpu_text)
        return cpu, int(rss_text) * 1024
    except (OSError, ValueError, subprocess.TimeoutExpired):
        return None, None


def run_with_metrics(command: list[str], log_path: Path, timeout: float | None = None,
                     sample_interval: float = 0.05) -> dict:
    """Run command, writing merged stdout/stderr to log_path.

    Peak RSS is sampled from the direct child using `ps`; very short-lived
    processes may finish between samples, in which case RSS is null. CPU time
    is the child's accumulated process time as reported by `ps`.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    process = None
    cpu_seconds = None
    peak_rss_bytes = None
    timed_out = False
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                   text=True)
        deadline = None if timeout is None else started + timeout
        while process.poll() is None:
            cpu, rss = _process_stats(process.pid)
            if cpu is not None:
                cpu_seconds = cpu
            if rss is not None:
                peak_rss_bytes = max(peak_rss_bytes or 0, rss)
            if deadline is not None and time.monotonic() >= deadline:
                timed_out = True
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                break
            time.sleep(sample_interval)
        # Capture a final sample in case the last poll preceded process exit.
        cpu, rss = _process_stats(process.pid)
        if cpu is not None:
            cpu_seconds = cpu
        if rss is not None:
            peak_rss_bytes = max(peak_rss_bytes or 0, rss)
    return {
        "exit_code": process.returncode,
        "timed_out": timed_out,
        "wall_seconds": time.monotonic() - started,
        "cpu_seconds_sampled": cpu_seconds,
        "peak_rss_bytes_sampled": peak_rss_bytes,
        "resource_sampling": "ps 50ms; peak RSS is a sampled lower bound",
    }


def source_tree_sha256(root: Path) -> str:
    """Fingerprint tracked and non-ignored Rust/Python/dependency source files."""
    result = subprocess.run(["git", "ls-files", "-co", "--exclude-standard"],
                            cwd=root, capture_output=True, text=True, check=False)
    if result.returncode:
        return "unavailable"
    wanted = {".rs", ".py", ".toml", ".lock"}
    digest = hashlib.sha256()
    for relative in sorted(path for path in result.stdout.splitlines()
                           if Path(path).suffix in wanted):
        path = root / relative
        if path.is_file():
            digest.update(relative.encode("utf-8"))
            digest.update(bytes.fromhex(file_sha256(path)))
    return digest.hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_identity(root: Path, binary: Path) -> dict:
    root = root.resolve()
    try:
        revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=root,
                                  capture_output=True, text=True, check=True,
                                  timeout=10).stdout.strip()
        dirty = subprocess.run(["git", "diff", "--quiet", "HEAD"], cwd=root,
                               check=False, timeout=10).returncode != 0
    except (OSError, subprocess.SubprocessError):
        revision, dirty = "unknown", None
    versions = {}
    for package in ("lightkurve", "numpy", "astropy", "transitleastsquares"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            continue

    def version(command: list[str]) -> str:
        try:
            result = subprocess.run(command, capture_output=True, text=True,
                                    check=True, timeout=10)
            return (result.stdout or result.stderr).strip()
        except (OSError, subprocess.SubprocessError):
            return "unavailable"

    return {
        "source_revision": revision, "source_tree_sha256": source_tree_sha256(root),
        "working_tree_dirty": dirty, "binary_path": str(binary.resolve()),
        "binary_sha256": file_sha256(binary), "rustc": version(["rustc", "--version"]),
        "cargo": version(["cargo", "--version"]), "python": sys.version.split()[0],
        "platform": platform.platform(), "logical_cpus": os.cpu_count(),
        "python_packages": versions,
    }

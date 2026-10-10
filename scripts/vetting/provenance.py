"""Provenance for vetting outputs: code revision, input checksums, versions.

A result that cannot be tied to a source revision cannot be reproduced, so
every gate report records the git commit and whether the tree had
uncommitted changes, alongside SHA-256 digests of its inputs.
"""
from __future__ import annotations

import hashlib
import platform
import subprocess
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PACKAGES = ("numpy", "scipy", "astropy", "astroquery", "batman-package", "emcee",
            "lightkurve", "triceratops")


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def git_state(root: Path = ROOT) -> dict:
    """Commit hash and dirty flag; `None` fields when git or the repo is absent."""
    def run(*args):
        try:
            done = subprocess.run(["git", "-C", str(root), *args], capture_output=True,
                                  text=True, timeout=10)
        except (OSError, subprocess.TimeoutExpired):
            return None
        return done.stdout.strip() if done.returncode == 0 else None
    commit = run("rev-parse", "HEAD")
    status = run("status", "--porcelain", "--untracked-files=no")
    return {"commit": commit, "dirty": None if status is None else bool(status)}


def package_versions(names=PACKAGES) -> dict:
    versions = {}
    for name in names:
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            continue
    return versions


def record(inputs: dict[str, Path], argv: list[str] | None = None) -> dict:
    """Provenance block: code revision, input digests, interpreter, packages."""
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git": git_state(),
        "inputs": {name: {"path": str(path), "sha256": sha256_file(path)}
                   for name, path in inputs.items() if path and Path(path).exists()},
        "argv": list(sys.argv if argv is None else argv),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": package_versions(),
    }

"""Run-manifest helpers used by experiment drivers."""

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit(repository_root):
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository_root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    return result.stdout.strip() if result.returncode == 0 else "unknown"


def git_is_dirty(repository_root):
    result = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=repository_root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    return result.returncode != 0 or bool(result.stdout.strip())


def portable_path(path, repository_root):
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(Path(repository_root).resolve()))
    except ValueError:
        return str(resolved)


def write_manifest(path, experiment, config, repository_root, outputs=()):
    payload = {
        "schema_version": 1,
        "experiment": experiment,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(repository_root),
        "git_dirty": git_is_dirty(repository_root),
        "command": sys.argv,
        "python": sys.version,
        "platform": platform.platform(),
        "config": config,
        "outputs": [
            {
                "path": portable_path(output, repository_root),
                "sha256": file_sha256(output),
            }
            for output in outputs
            if Path(output).is_file()
        ],
    }
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return payload

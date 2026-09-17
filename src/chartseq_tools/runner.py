from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, TextIO


LOG = logging.getLogger("chartseq")


def executable(command: str) -> str:
    expanded = str(Path(command).expanduser()) if "/" in command else command
    found = shutil.which(expanded)
    if found:
        return found
    candidate = Path(expanded)
    if candidate.is_file():
        return str(candidate.resolve())
    raise FileNotFoundError(f"Required executable not found: {command}")


def format_command(command: Iterable[str]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def run(
    command: list[str],
    *,
    stdout: TextIO | None = None,
    cwd: Path | None = None,
    dry_run: bool = False,
) -> None:
    command = [str(part) for part in command]
    LOG.info("$ %s", format_command(command))
    if dry_run:
        return
    subprocess.run(command, check=True, stdout=stdout, cwd=cwd)


def run_to_file(command: list[str], output: Path, *, dry_run: bool = False) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if dry_run:
        run(command, dry_run=True)
        return
    with output.open("wb") as handle:
        run(command, stdout=handle)

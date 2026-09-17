from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any


CORE_TOOLS = ("umi_tools", "cutadapt", "star", "featurecounts", "samtools")


def check_environment(config: dict[str, Any]) -> list[tuple[str, str, bool]]:
    checks = []
    for key in CORE_TOOLS:
        value = str(config["tools"].get(key, key))
        candidate = shutil.which(value) or (str(Path(value).expanduser().resolve()) if Path(value).is_file() else "")
        checks.append((key, candidate or value, bool(candidate)))
    for key in ("star_index", "gtf"):
        value = config.get("references", {}).get(key, "")
        checks.append((f"references.{key}", str(value), bool(value and Path(value).expanduser().exists())))
    return checks


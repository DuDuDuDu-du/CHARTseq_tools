from __future__ import annotations

from copy import deepcopy
from importlib.resources import files
from pathlib import Path
from typing import Any

import yaml


DEFAULTS: dict[str, Any] = {
    "references": {},
    "tools": {
        "umi_tools": "umi_tools",
        "fastqc": "fastqc",
        "cutadapt": "cutadapt",
        "star": "STAR",
        "featurecounts": "featureCounts",
        "samtools": "samtools",
        "seqtk": "seqtk",
        "picard": "picard",
        "gene_body_coverage": "geneBody_coverage.py",
        "salmon": "salmon",
        "suppa": "suppa.py",
    },
    "library": {
        "bc_pattern": "CCCCCCCCXXXXXXXXXXXXXXXXXXX",
        "bc_pattern2": "CCCCCCCCNNNNNNXXXXXXXXXXXXXXXXXXX",
        "pattern_method": "string",
        "post_extract_trim_r1": 19,
        "post_extract_trim_r2": 19,
    },
    "pipeline": {
        "threads": 8,
        "minimum_length": 20,
        "quality_cutoff": "20,20",
        "adapters": [
            "AGATCGGAAGAG",
            "CTGTCTCTTATACACATCT",
            "A{7}", "T{7}", "C{7}", "G{7}",
        ],
        "run_fastqc": True,
        "keep_intermediates": False,
    },
}


def _merge(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _merge(base[key], value)
        else:
            base[key] = value
    return base


def load_config(path: str | Path) -> dict[str, Any]:
    config_path = Path(path).expanduser().resolve()
    with config_path.open(encoding="utf-8") as handle:
        user = yaml.safe_load(handle) or {}
    if not isinstance(user, dict):
        raise ValueError("The configuration root must be a YAML mapping.")
    return _merge(deepcopy(DEFAULTS), user)


def write_example(path: str | Path, force: bool = False) -> Path:
    destination = Path(path).expanduser().resolve()
    if destination.exists() and not force:
        raise FileExistsError(f"Refusing to overwrite existing file: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    source = files("chartseq_tools").joinpath("resources/config.example.yaml")
    destination.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    return destination


def require_reference(config: dict[str, Any], name: str) -> Path:
    value = config.get("references", {}).get(name)
    if not value:
        raise ValueError(f"Missing required references.{name} in the configuration.")
    path = Path(value).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Reference does not exist: {path}")
    return path

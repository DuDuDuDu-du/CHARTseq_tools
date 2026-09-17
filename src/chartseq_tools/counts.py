from __future__ import annotations

from pathlib import Path

import pandas as pd


def pivot_counts(input_file: str | Path, output_file: str | Path) -> Path:
    """Convert the long UMI-tools count table to a gene-by-cell matrix."""
    source, destination = Path(input_file), Path(output_file)
    frame = pd.read_csv(source, sep="\t")
    frame = frame.rename(columns={"gene": "gene", "cell": "cell_barcode"})
    required = {"gene", "cell_barcode", "count"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing columns in {source}: {', '.join(sorted(missing))}")
    # pivot_table safely combines duplicate gene/cell records if they occur.
    matrix = frame.pivot_table(
        index="gene", columns="cell_barcode", values="count", aggfunc="sum", fill_value=0
    )
    matrix.index.name = "gene"
    destination.parent.mkdir(parents=True, exist_ok=True)
    matrix.reset_index().to_csv(destination, sep="\t", index=False)
    return destination


def summarize_matrices(
    input_files: list[str | Path], output_file: str | Path, threshold: float = 0
) -> Path:
    """Count genes above *threshold* per cell in one or more count matrices."""
    results: dict[str, dict[str, int]] = {}
    for item in input_files:
        path = Path(item)
        frame = pd.read_csv(path, sep="\t")
        if frame.shape[1] < 2:
            raise ValueError(f"Count matrix has no cell columns: {path}")
        label = path.name
        for barcode in frame.columns[1:]:
            numeric = pd.to_numeric(frame[barcode], errors="coerce").fillna(0)
            results.setdefault(barcode, {})[label] = int((numeric > threshold).sum())
    summary = pd.DataFrame.from_dict(results, orient="index").fillna(0).astype(int)
    summary.index.name = "cell_barcode"
    destination = Path(output_file)
    destination.parent.mkdir(parents=True, exist_ok=True)
    summary.sort_index().reset_index().to_csv(destination, sep="\t", index=False)
    return destination


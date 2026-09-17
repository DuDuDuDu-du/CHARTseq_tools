from __future__ import annotations

from pathlib import Path

import pandas as pd


def plot_saturation(input_file: str | Path, output_file: str | Path, metric: str) -> Path:
    try:
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError("Plotting requires: pip install 'CHARTseq_tools[plot]'") from error

    frame = pd.read_csv(input_file, sep="\t", index_col=0)
    frame.index = pd.to_numeric(frame.index, errors="coerce")
    frame = frame.loc[frame.index.notna()].apply(pd.to_numeric, errors="coerce")
    if frame.empty:
        raise ValueError("No numeric saturation data found.")
    frame = frame.sort_index()
    mean, low, high = frame.mean(axis=1), frame.min(axis=1), frame.max(axis=1)
    label = "Detected genes" if metric == "genes" else "UMI count"
    destination = Path(output_file)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig, axis = plt.subplots(figsize=(8, 5))
    axis.fill_between(frame.index, low, high, color="#b2abd2", alpha=0.55, linewidth=0)
    axis.plot(frame.index, mean, "o-", color="#5e3c99", linewidth=1.8, markersize=4)
    axis.set(xlabel="Downsampled read pairs (millions)", ylabel=label, title=f"{label} saturation")
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", linestyle="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(destination, dpi=300)
    plt.close(fig)
    return destination


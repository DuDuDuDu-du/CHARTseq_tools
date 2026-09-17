from pathlib import Path

import pandas as pd

from chartseq_tools.counts import pivot_counts, summarize_matrices


def test_pivot_and_summary(tmp_path: Path):
    source = tmp_path / "counts.tsv"
    source.write_text("gene\tcell\tcount\nG1\tC1\t2\nG1\tC2\t0\nG2\tC1\t3\n", encoding="utf-8")
    matrix = pivot_counts(source, tmp_path / "matrix.tsv")
    frame = pd.read_csv(matrix, sep="\t")
    assert list(frame.columns) == ["gene", "C1", "C2"]
    assert frame.loc[frame["gene"] == "G1", "C1"].item() == 2
    summary = summarize_matrices([matrix], tmp_path / "summary.tsv", threshold=0)
    result = pd.read_csv(summary, sep="\t").set_index("cell_barcode")
    assert result.loc["C1"].iloc[0] == 2
    assert result.loc["C2"].iloc[0] == 0


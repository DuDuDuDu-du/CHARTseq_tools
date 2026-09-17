from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from .runner import executable, run_to_file


LOG = logging.getLogger("chartseq")


def downsample_pairs(
    read1: str | Path,
    read2: str | Path,
    output_dir: str | Path,
    depths_millions: list[float],
    *,
    seqtk: str = "seqtk",
    seed: int = 100,
    dry_run: bool = False,
) -> list[tuple[Path, Path]]:
    """Downsample paired reads with the same random seed to preserve pairing."""
    program = executable(seqtk)
    output = Path(output_dir).resolve()
    generated = []
    for depth in sorted(set(depths_millions)):
        reads = int(depth * 1_000_000)
        if reads <= 0:
            raise ValueError(f"Downsampling depth must be positive: {depth}")
        target = output / f"{depth:g}M_reads"
        r1_out, r2_out = target / "sampled_R1.fastq", target / "sampled_R2.fastq"
        run_to_file([program, "sample", "-s", str(seed), str(read1), str(reads)], r1_out, dry_run=dry_run)
        run_to_file([program, "sample", "-s", str(seed), str(read2), str(reads)], r2_out, dry_run=dry_run)
        generated.append((r1_out, r2_out))
    return generated


def summarize_downsampling(
    input_dir: str | Path,
    output_dir: str | Path,
    matrix_name: str = "sampled.gene_name.matrix.tsv",
) -> tuple[Path, Path]:
    root, output = Path(input_dir), Path(output_dir)
    rows_genes: dict[float, dict[str, int]] = {}
    rows_umis: dict[float, dict[str, float]] = {}
    for directory in root.glob("*M_reads"):
        try:
            level = float(directory.name.removesuffix("M_reads"))
        except ValueError:
            continue
        matrix_path = directory / matrix_name
        if not matrix_path.exists():
            LOG.warning("Skipping %s: %s not found", directory, matrix_name)
            continue
        frame = pd.read_csv(matrix_path, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
        rows_genes[level] = (frame > 0).sum().astype(int).to_dict()
        rows_umis[level] = frame.sum().to_dict()
    if not rows_genes:
        raise ValueError(f"No downsampling matrices named {matrix_name!r} were found under {root}")
    output.mkdir(parents=True, exist_ok=True)
    gene_path, umi_path = output / "detected_genes.tsv", output / "umi_counts.tsv"
    pd.DataFrame.from_dict(rows_genes, orient="index").sort_index().rename_axis("sampling_level").to_csv(gene_path, sep="\t")
    pd.DataFrame.from_dict(rows_umis, orient="index").sort_index().rename_axis("sampling_level").to_csv(umi_path, sep="\t")
    return gene_path, umi_path


def run_downsample_series(
    read1: str | Path,
    read2: str | Path,
    output_dir: str | Path,
    depths_millions: list[float],
    *,
    sample: str,
    cell_number: int,
    config: dict,
    seqtk: str = "seqtk",
    seed: int = 100,
    jobs: int = 1,
    force: bool = False,
) -> None:
    """Downsample paired FASTQs and run the CHART-seq workflow at every depth."""
    from .pipeline import run_pipeline

    pairs = downsample_pairs(read1, read2, output_dir, depths_millions, seqtk=seqtk, seed=seed)
    if jobs < 1:
        raise ValueError("jobs must be at least 1")
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {
            pool.submit(run_pipeline, r1, r2, sample, r1.parent, cell_number, config, force=force): r1.parent.name
            for r1, r2 in pairs
        }
        for future in as_completed(futures):
            level = futures[future]
            future.result()
            LOG.info("Completed downsampling analysis: %s", level)

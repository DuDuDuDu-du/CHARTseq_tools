from __future__ import annotations

from pathlib import Path

import pandas as pd


def normalize_featurecounts(
    counts_file: str | Path,
    summary_file: str | Path,
    output_file: str | Path,
    *,
    threshold: float = 1.0,
    paired: bool = False,
) -> Path:
    """Calculate RPKM/FPKM from featureCounts output and filter genes."""
    counts = pd.read_csv(counts_file, sep="\t", comment="#")
    sample_columns = [column for column in counts.columns if column.endswith("bam") or ".bam" in column]
    if not sample_columns or not {"Geneid", "Length"}.issubset(counts.columns):
        raise ValueError("featureCounts table must contain Geneid, Length and at least one BAM column")
    summary = pd.read_csv(summary_file, sep="\t", index_col=0)
    if "Assigned" not in summary.index:
        raise ValueError("featureCounts summary does not contain an Assigned row")
    assigned = pd.to_numeric(summary.loc["Assigned", sample_columns], errors="raise")
    values = counts[sample_columns].div(counts["Length"], axis=0).div(assigned, axis=1) * 1e9
    metric = "FPKM" if paired else "RPKM"
    keep = (values >= threshold).any(axis=1)
    result = counts.loc[keep, ["Geneid", "Length", *sample_columns]].copy()
    result = result.rename(columns={"Length": "Gene_Length"})
    for column in sample_columns:
        result[f"{metric}_{column}"] = values.loc[keep, column]
    destination = Path(output_file)
    destination.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(destination, sep="\t", index=False)
    return destination


def isoform_summary(
    quant_file: str | Path,
    filtered_genes_file: str | Path,
    gene_map_file: str | Path,
    transcript_lengths_file: str | Path,
    sample: str,
    output_dir: str | Path,
    *,
    tpm_cutoff: float = 0.1,
) -> Path:
    quant = pd.read_csv(quant_file, sep="\t")
    if not {"Name", "TPM"}.issubset(quant.columns):
        raise ValueError("Salmon quant file must contain Name and TPM")
    genes = set(pd.read_csv(filtered_genes_file, sep="\t")["Geneid"].astype(str))
    mapping = pd.read_csv(
        gene_map_file, sep="\t", header=None, names=["transcript_id", "gene_id", "gene_name"]
    ).drop_duplicates("transcript_id")
    detected = quant.merge(mapping, left_on="Name", right_on="transcript_id", how="inner")
    detected = detected[detected["gene_id"].astype(str).isin(genes) & (detected["TPM"] >= tpm_cutoff)].copy()
    if detected.empty:
        raise ValueError("No transcripts remain after gene and TPM filtering")
    lengths = pd.read_csv(transcript_lengths_file, sep="\t")
    detected = detected.merge(lengths, on="transcript_id", how="left")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    with (output / "iso_tpm_filtered.txt").open("w", encoding="utf-8") as handle:
        handle.write(sample + "\n")
        detected.sort_values("TPM", ascending=False)[["transcript_id", "TPM"]].to_csv(
            handle, sep="\t", index=False, header=False, float_format="%.6f"
        )
    per_gene = detected.groupby(["gene_id", "gene_name"], dropna=False)["transcript_id"].nunique().rename("isoform_count").reset_index()
    per_gene.sort_values("isoform_count", ascending=False).to_csv(output / "gene_isoform_count.tsv", sep="\t", index=False)
    detected.to_csv(output / "detected_transcripts.tsv", sep="\t", index=False)
    length = pd.to_numeric(detected.get("length"), errors="coerce").dropna()
    summary = pd.DataFrame([{
        "sample": sample,
        "transcript_filter": f"TPM>={tpm_cutoff:g}",
        "detected_isoforms": detected["transcript_id"].nunique(),
        "detected_genes": len(per_gene),
        "mean_isoforms_per_gene": per_gene["isoform_count"].mean(),
        "median_isoforms_per_gene": per_gene["isoform_count"].median(),
        "multi_isoform_ratio_ge2": (per_gene["isoform_count"] >= 2).mean() * 100,
        "multi_isoform_ratio_ge3": (per_gene["isoform_count"] >= 3).mean() * 100,
        "mean_transcript_length": length.mean(),
        "median_transcript_length": length.median(),
        "transcript_gt2kb": (length > 2000).mean() * 100,
        "transcript_gt5kb": (length > 5000).mean() * 100,
    }])
    destination = output / "benchmark_summary.tsv"
    summary.to_csv(destination, sep="\t", index=False)
    return destination


def summarize_splicing(sample: str, input_dir: str | Path, output_dir: str | Path) -> Path:
    source, output = Path(input_dir), Path(output_dir)
    events = ["SE", "RI", "A3SS", "A5SS", "MXE", "AF", "AL"]
    rows, detected_frames = [], []
    for event in events:
        path = source / f"{event}.psi"
        if not path.exists():
            continue
        frame = pd.read_csv(path, sep="\t")
        detected = frame[frame.iloc[:, -1].notna()].copy()
        detected["event_type"] = event
        detected_frames.append(detected)
        rows.append({"sample": sample, "event_type": event, "detected_events": len(detected)})
    if not rows:
        raise ValueError(f"No SUPPA PSI files found in {source}")
    output.mkdir(parents=True, exist_ok=True)
    long = pd.DataFrame(rows)
    wide = long.pivot(index="sample", columns="event_type", values="detected_events").fillna(0).astype(int)
    wide["total_AS_events"] = wide.sum(axis=1)
    long.to_csv(output / "AS_event_summary_long.tsv", sep="\t", index=False)
    pd.concat(detected_frames, ignore_index=True).to_csv(output / "detected_AS_events.tsv", sep="\t", index=False)
    destination = output / "AS_event_summary.tsv"
    wide.reset_index().to_csv(destination, sep="\t", index=False)
    return destination


from __future__ import annotations

import gzip
import logging
import shutil
from pathlib import Path
from typing import Any

from .config import require_reference
from .counts import pivot_counts, summarize_matrices
from .runner import executable, run


LOG = logging.getLogger("chartseq")


def _tool(config: dict[str, Any], key: str) -> str:
    return executable(config["tools"].get(key, key))


def _completed(target: Path, force: bool) -> bool:
    if target.exists() and target.stat().st_size > 0 and not force:
        LOG.info("Skipping completed output: %s", target)
        return True
    return False


def _move_featurecounts_bam(aligned: Path, destination: Path, dry_run: bool) -> None:
    generated = Path(str(aligned) + ".featureCounts.bam")
    if dry_run:
        LOG.info("Would move %s -> %s", generated, destination)
    else:
        if not generated.exists():
            raise FileNotFoundError(f"featureCounts did not create {generated}")
        generated.replace(destination)


def run_pipeline(
    read1: str | Path,
    read2: str | Path,
    sample: str,
    output_dir: str | Path,
    cell_number: int,
    config: dict[str, Any],
    *,
    force: bool = False,
    dry_run: bool = False,
) -> dict[str, Path]:
    """Run the paired-end CHART-seq preprocessing and UMI counting workflow."""
    if cell_number < 1:
        raise ValueError("cell_number must be at least 1")
    r1, r2 = Path(read1).resolve(), Path(read2).resolve()
    if not dry_run:
        for source in (r1, r2):
            if not source.is_file():
                raise FileNotFoundError(source)
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    threads = int(config["pipeline"]["threads"])
    library, settings = config["library"], config["pipeline"]
    star_index, gtf = require_reference(config, "star_index"), require_reference(config, "gtf")

    umi_tools = _tool(config, "umi_tools")
    cutadapt = _tool(config, "cutadapt")
    star = _tool(config, "star")
    featurecounts = _tool(config, "featurecounts")
    samtools = _tool(config, "samtools")

    whitelist = output / "whitelist.tsv"
    extracted1, extracted2 = output / f"{sample}.extracted_R1.fastq.gz", output / f"{sample}.extracted_R2.fastq.gz"
    trimmed1, trimmed2 = output / f"{sample}.trimmed_R1.fastq.gz", output / f"{sample}.trimmed_R2.fastq.gz"
    star_prefix = output / f"{sample}.star."
    aligned = output / f"{sample}.star.Aligned.sortedByCoord.out.bam"

    pattern_args = [
        f"--bc-pattern={library['bc_pattern']}",
        f"--bc-pattern2={library['bc_pattern2']}",
        f"--extract-method={library.get('pattern_method', 'string')}",
    ]
    if not _completed(whitelist, force):
        whitelist_command = [
            umi_tools, "whitelist", "--stdin", str(r1), "--read2-in", str(r2),
            *pattern_args, f"--set-cell-number={cell_number}", "--error-correct-threshold=0",
        ]
        if dry_run:
            run(whitelist_command, dry_run=True)
        else:
            with whitelist.open("w", encoding="utf-8") as handle:
                run(whitelist_command, stdout=handle)

    if not (_completed(extracted1, force) and _completed(extracted2, force)):
        run([
            umi_tools, "extract", "--stdin", str(r1), "--read2-in", str(r2),
            *pattern_args, f"--whitelist={whitelist}", f"--stdout={extracted1}", f"--read2-out={extracted2}",
        ], dry_run=dry_run)

    if not (_completed(trimmed1, force) and _completed(trimmed2, force)):
        command = [cutadapt]
        trim1, trim2 = int(library.get("post_extract_trim_r1", 0)), int(library.get("post_extract_trim_r2", 0))
        if trim1:
            command += ["-u", str(trim1)]
        if trim2:
            command += ["-U", str(trim2)]
        for adapter in settings["adapters"]:
            command += ["-a", adapter, "-A", adapter]
        command += [
            "--minimum-length", str(settings["minimum_length"]),
            "--quality-cutoff", str(settings["quality_cutoff"]),
            "-j", str(threads), "-o", str(trimmed1), "-p", str(trimmed2), str(extracted1), str(extracted2),
        ]
        run(command, dry_run=dry_run)

    if settings.get("run_fastqc"):
        fastqc = _tool(config, "fastqc")
        run([fastqc, "--threads", str(min(threads, 4)), "--outdir", str(output), str(trimmed1), str(trimmed2)], dry_run=dry_run)

    if not _completed(aligned, force):
        run([
            star, "--readFilesIn", str(trimmed1), str(trimmed2), "--readFilesCommand", "zcat",
            "--genomeDir", str(star_index), "--outFileNamePrefix", str(star_prefix),
            "--runThreadN", str(threads), "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outSAMattributes", "NH", "HI", "AS", "nM", "XS", "--quantMode", "TranscriptomeSAM", "GeneCounts",
        ], dry_run=dry_run)

    matrices: list[Path] = []
    outputs: dict[str, Path] = {"whitelist": whitelist, "aligned_bam": aligned}
    for label, gene_attribute in (("gene_name", "gene_name"), ("gene_id", "gene_id")):
        assigned = output / f"{sample}.{label}.assigned.bam"
        sorted_bam = output / f"{sample}.{label}.assigned.sorted.bam"
        count_gz = output / f"{sample}.{label}.counts.tsv.gz"
        matrix = output / f"{sample}.{label}.matrix.tsv"
        if not _completed(assigned, force):
            run([
                featurecounts, "-a", str(gtf), "-o", str(output / f"{sample}.{label}.featureCounts.txt"),
                "-R", "BAM", str(aligned), "-T", str(threads), "-p", "--countReadPairs", "-g", gene_attribute,
            ], dry_run=dry_run)
            _move_featurecounts_bam(aligned, assigned, dry_run)
        if not _completed(sorted_bam, force):
            run([samtools, "sort", "-@", str(threads), "-o", str(sorted_bam), str(assigned)], dry_run=dry_run)
            run([samtools, "index", "-@", str(threads), str(sorted_bam)], dry_run=dry_run)
        if not _completed(count_gz, force):
            run([
                umi_tools, "count", "--per-gene", "--gene-tag=XT", "--assigned-status-tag=XS", "--per-cell",
                "-I", str(sorted_bam), "-S", str(count_gz),
            ], dry_run=dry_run)
        if not dry_run and not _completed(matrix, force):
            with gzip.open(count_gz, "rt") as source, (output / f"{sample}.{label}.counts.tsv").open("w") as destination:
                shutil.copyfileobj(source, destination)
            pivot_counts(output / f"{sample}.{label}.counts.tsv", matrix)
        matrices.append(matrix)
        outputs[f"{label}_matrix"] = matrix

    summary = output / f"{sample}.summary.tsv"
    if not dry_run:
        summarize_matrices(matrices, summary)
    outputs["summary"] = summary

    if not settings.get("keep_intermediates") and not dry_run:
        for path in (extracted1, extracted2):
            path.unlink(missing_ok=True)
    return outputs

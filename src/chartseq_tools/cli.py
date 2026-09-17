from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from . import __version__
from .bam import split_bam
from .benchmark import isoform_summary, normalize_featurecounts, summarize_splicing
from .config import load_config, write_example
from .counts import pivot_counts, summarize_matrices
from .doctor import check_environment
from .downsample import downsample_pairs, run_downsample_series, summarize_downsampling
from .fastq import split_fastq, trim_prefix_if_matches
from .pipeline import run_pipeline
from .saturation import plot_saturation


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(prog="chartseq", description="Portable CHART-seq analysis toolkit")
    root.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    root.add_argument("-v", "--verbose", action="count", default=0)
    commands = root.add_subparsers(dest="command", required=True)

    item = commands.add_parser("init-config", help="write an editable example configuration")
    item.add_argument("output", nargs="?", default="chartseq.yaml")
    item.add_argument("--force", action="store_true")

    item = commands.add_parser("doctor", help="check core executables and references")
    item.add_argument("--config", required=True)

    item = commands.add_parser("run", help="run the paired-end CHART-seq workflow")
    item.add_argument("--read1", required=True)
    item.add_argument("--read2", required=True)
    item.add_argument("--sample", required=True)
    item.add_argument("--output", required=True)
    item.add_argument("--cells", required=True, type=int)
    item.add_argument("--config", required=True)
    item.add_argument("--force", action="store_true")
    item.add_argument("--dry-run", action="store_true")

    item = commands.add_parser("split-fastq", help="split paired FASTQs using barcodes in read names")
    item.add_argument("--read1", required=True)
    item.add_argument("--read2", required=True)
    item.add_argument("--whitelist", required=True)
    item.add_argument("--output", required=True)
    item.add_argument("--header-regex", help="regex with a named 'barcode' group or first capture group")
    item.add_argument("--header-field", type=int, default=-2)
    item.add_argument("--max-open-files", type=int, default=128)
    item.add_argument("--plain", action="store_true", help="write uncompressed FASTQ")

    item = commands.add_parser("trim-prefix", help="conditionally trim a prefix from FASTQ reads")
    item.add_argument("input")
    item.add_argument("output")
    item.add_argument("--prefix", required=True)
    item.add_argument("--length", required=True, type=int)

    item = commands.add_parser("split-bam", help="split BAM using barcodes in query names")
    item.add_argument("input_bam")
    item.add_argument("output")
    item.add_argument("--query-regex", help="regex with a named 'barcode' group or first capture group")
    item.add_argument("--query-field", type=int, default=-2)
    item.add_argument("--max-open-files", type=int, default=128)
    item.add_argument("--index", action="store_true")

    item = commands.add_parser("pivot-counts", help="convert long UMI counts to a gene-by-cell matrix")
    item.add_argument("input")
    item.add_argument("output")

    item = commands.add_parser("summarize-counts", help="summarize detected genes per cell")
    item.add_argument("inputs", nargs="+")
    item.add_argument("--output", required=True)
    item.add_argument("--threshold", type=float, default=0)

    item = commands.add_parser("downsample", help="create paired FASTQ subsets with seqtk")
    item.add_argument("--read1", required=True)
    item.add_argument("--read2", required=True)
    item.add_argument("--output", required=True)
    item.add_argument("--depths", nargs="+", type=float, required=True, metavar="MILLIONS")
    item.add_argument("--seqtk", default="seqtk")
    item.add_argument("--seed", type=int, default=100)
    item.add_argument("--dry-run", action="store_true")

    item = commands.add_parser("saturation-run", help="downsample and analyze every sequencing depth")
    item.add_argument("--read1", required=True)
    item.add_argument("--read2", required=True)
    item.add_argument("--output", required=True)
    item.add_argument("--depths", nargs="+", type=float, required=True, metavar="MILLIONS")
    item.add_argument("--sample", default="sampled")
    item.add_argument("--cells", required=True, type=int)
    item.add_argument("--config", required=True)
    item.add_argument("--seqtk", default="seqtk")
    item.add_argument("--seed", type=int, default=100)
    item.add_argument("--jobs", type=int, default=1)
    item.add_argument("--force", action="store_true")

    item = commands.add_parser("summarize-downsample", help="summarize gene and UMI saturation matrices")
    item.add_argument("input")
    item.add_argument("output")
    item.add_argument("--matrix-name", default="sampled.gene_name.matrix.tsv")

    item = commands.add_parser("plot-saturation", help="plot downsampling saturation")
    item.add_argument("input")
    item.add_argument("output")
    item.add_argument("--metric", choices=["genes", "umis"], required=True)

    item = commands.add_parser("normalize-featurecounts", help="calculate RPKM/FPKM and filter genes")
    item.add_argument("counts")
    item.add_argument("summary")
    item.add_argument("output")
    item.add_argument("--threshold", type=float, default=1)
    item.add_argument("--paired", action="store_true", help="label the metric FPKM instead of RPKM")

    item = commands.add_parser("isoform-summary", help="filter Salmon isoforms and create benchmark metrics")
    item.add_argument("--quant", required=True)
    item.add_argument("--genes", required=True)
    item.add_argument("--gene-map", required=True)
    item.add_argument("--transcript-lengths", required=True)
    item.add_argument("--sample", required=True)
    item.add_argument("--output", required=True)
    item.add_argument("--tpm-cutoff", type=float, default=0.1)

    item = commands.add_parser("summarize-splicing", help="summarize SUPPA PSI event files")
    item.add_argument("--sample", required=True)
    item.add_argument("--input", required=True)
    item.add_argument("--output", required=True)
    return root


def main(argv: list[str] | None = None) -> int:
    arguments = parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if arguments.verbose > 1 else logging.INFO if arguments.verbose else logging.WARNING,
        format="%(levelname)s: %(message)s",
    )
    try:
        command = arguments.command
        if command == "init-config":
            print(write_example(arguments.output, arguments.force))
        elif command == "doctor":
            checks = check_environment(load_config(arguments.config))
            for name, value, ok in checks:
                print(f"{'OK' if ok else 'MISSING'}\t{name}\t{value}")
            return 0 if all(row[2] for row in checks) else 1
        elif command == "run":
            outputs = run_pipeline(arguments.read1, arguments.read2, arguments.sample, arguments.output, arguments.cells,
                                   load_config(arguments.config), force=arguments.force, dry_run=arguments.dry_run)
            for name, path in outputs.items():
                print(f"{name}\t{path}")
        elif command == "split-fastq":
            counts = split_fastq(arguments.read1, arguments.read2, arguments.whitelist, arguments.output,
                                  header_regex=arguments.header_regex, header_field=arguments.header_field,
                                  max_open_files=arguments.max_open_files, gzip_output=not arguments.plain)
            for barcode, count in counts.most_common():
                print(f"{barcode}\t{count}")
        elif command == "split-bam":
            counts = split_bam(arguments.input_bam, arguments.output, query_regex=arguments.query_regex,
                               query_field=arguments.query_field, max_open_files=arguments.max_open_files,
                               index_outputs=arguments.index)
            for barcode, count in sorted(counts.items()):
                print(f"{barcode}\t{count}")
        elif command == "trim-prefix":
            total, trimmed = trim_prefix_if_matches(arguments.input, arguments.output, arguments.prefix, arguments.length)
            print(f"total_reads\t{total}\ntrimmed_reads\t{trimmed}")
        elif command == "pivot-counts":
            print(pivot_counts(arguments.input, arguments.output))
        elif command == "summarize-counts":
            print(summarize_matrices(arguments.inputs, arguments.output, arguments.threshold))
        elif command == "downsample":
            for pair in downsample_pairs(arguments.read1, arguments.read2, arguments.output, arguments.depths,
                                         seqtk=arguments.seqtk, seed=arguments.seed, dry_run=arguments.dry_run):
                print(*pair, sep="\t")
        elif command == "saturation-run":
            run_downsample_series(
                arguments.read1, arguments.read2, arguments.output, arguments.depths,
                sample=arguments.sample, cell_number=arguments.cells, config=load_config(arguments.config),
                seqtk=arguments.seqtk, seed=arguments.seed, jobs=arguments.jobs, force=arguments.force,
            )
        elif command == "summarize-downsample":
            print(*summarize_downsampling(arguments.input, arguments.output, arguments.matrix_name), sep="\n")
        elif command == "plot-saturation":
            print(plot_saturation(arguments.input, arguments.output, arguments.metric))
        elif command == "normalize-featurecounts":
            print(normalize_featurecounts(arguments.counts, arguments.summary, arguments.output,
                                          threshold=arguments.threshold, paired=arguments.paired))
        elif command == "isoform-summary":
            print(isoform_summary(arguments.quant, arguments.genes, arguments.gene_map, arguments.transcript_lengths,
                                  arguments.sample, arguments.output, tpm_cutoff=arguments.tpm_cutoff))
        elif command == "summarize-splicing":
            print(summarize_splicing(arguments.sample, arguments.input, arguments.output))
        return 0
    except (FileNotFoundError, ValueError, RuntimeError) as error:
        print(f"chartseq: error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

import gzip
from pathlib import Path

from chartseq_tools.fastq import barcode_from_header, split_fastq, trim_prefix_if_matches


def test_barcode_from_header():
    header = "@instrument_read_ACGTACGT_UMI 1:N:0:1\n"
    assert barcode_from_header(header) == "ACGTACGT"
    assert barcode_from_header(header, r"_(?P<barcode>[ACGT]{8})_") == "ACGTACGT"


def test_split_fastq(tmp_path: Path):
    r1, r2 = tmp_path / "R1.fastq", tmp_path / "R2.fastq"
    record1 = "@read_ACGT_UMI 1:N:0:1\nAAAA\n+\nIIII\n"
    record2 = "@read_ACGT_UMI 2:N:0:1\nTTTT\n+\nIIII\n"
    r1.write_text(record1, encoding="utf-8")
    r2.write_text(record2, encoding="utf-8")
    whitelist = tmp_path / "whitelist.tsv"
    whitelist.write_text("ACGT\n", encoding="utf-8")
    counts = split_fastq(r1, r2, whitelist, tmp_path / "out")
    assert counts["ACGT"] == 1
    with gzip.open(tmp_path / "out/ACGT_R1.fastq.gz", "rt") as handle:
        assert handle.read() == record1


def test_trim_prefix(tmp_path: Path):
    source, output = tmp_path / "input.fastq", tmp_path / "output.fastq"
    source.write_text("@r\nAACCGG\n+\nIIIIII\n", encoding="utf-8")
    total, trimmed = trim_prefix_if_matches(source, output, "AACC", 4)
    assert (total, trimmed) == (1, 1)
    assert output.read_text(encoding="utf-8") == "@r\nGG\n+\nII\n"

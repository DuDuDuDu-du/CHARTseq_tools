from __future__ import annotations

import re
import tempfile
from collections import OrderedDict
from pathlib import Path


def split_bam(
    input_bam: str | Path,
    output_dir: str | Path,
    *,
    query_regex: str | None = None,
    query_field: int = -2,
    max_open_files: int = 128,
    index_outputs: bool = False,
) -> dict[str, int]:
    try:
        import pysam
    except ImportError as error:
        raise RuntimeError("split-bam requires: pip install 'CHARTseq_tools[bam]'") from error

    source, output = Path(input_bam), Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    handles: OrderedDict[str, object] = OrderedDict()
    counts: dict[str, int] = {}

    # Spool SAM records to temporary text files through an LRU cache, then
    # convert each file to BAM. This remains valid when thousands of barcodes
    # exceed the process file-descriptor limit.
    with tempfile.TemporaryDirectory(prefix="chartseq_split_bam_") as temp_name, \
            pysam.AlignmentFile(source, "rb") as incoming:
        temp_dir = Path(temp_name)

        def writer(barcode: str):
            if barcode in handles:
                handles.move_to_end(barcode)
                return handles[barcode]
            if len(handles) >= max_open_files:
                _, old = handles.popitem(last=False)
                old.close()
            handle = (temp_dir / f"{barcode}.sam.txt").open("a", encoding="utf-8")
            handles[barcode] = handle
            return handle

        try:
            for read in incoming.fetch(until_eof=True):
                if query_regex:
                    match = re.search(query_regex, read.query_name)
                    barcode = match.group("barcode") if match and "barcode" in match.groupdict() else (
                        match.group(1) if match else None
                    )
                else:
                    fields = read.query_name.split("_")
                    try:
                        barcode = fields[query_field]
                    except IndexError:
                        barcode = None
                if barcode:
                    if not re.fullmatch(r"[A-Za-z0-9._-]+", barcode):
                        raise ValueError(f"Unsafe barcode for an output filename: {barcode!r}")
                    writer(barcode).write(read.to_string() + "\n")
                    counts[barcode] = counts.get(barcode, 0) + 1
        finally:
            for handle in handles.values():
                handle.close()

        for barcode in counts:
            target = output / f"{barcode}.bam"
            with pysam.AlignmentFile(target, "wb", template=incoming) as outgoing, \
                    (temp_dir / f"{barcode}.sam.txt").open(encoding="utf-8") as records:
                for line in records:
                    outgoing.write(pysam.AlignedSegment.fromstring(line.rstrip("\n"), incoming.header))

    if index_outputs:
        for barcode in counts:
            pysam.index(str(output / f"{barcode}.bam"))
    return counts

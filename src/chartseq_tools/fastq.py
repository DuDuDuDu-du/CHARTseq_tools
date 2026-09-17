from __future__ import annotations

import gzip
import re
from collections import Counter, OrderedDict
from contextlib import ExitStack
from pathlib import Path
from typing import Iterator, TextIO


def open_text(path: str | Path, mode: str = "rt") -> TextIO:
    path = Path(path)
    return gzip.open(path, mode) if path.suffix == ".gz" else path.open(mode)


def records(handle: TextIO) -> Iterator[tuple[str, str, str, str]]:
    while True:
        record = tuple(handle.readline() for _ in range(4))
        if not record[0]:
            return
        if any(line == "" for line in record):
            raise ValueError("Truncated FASTQ record encountered.")
        yield record  # type: ignore[misc]


def barcode_from_header(header: str, regex: str | None = None, field: int = -2) -> str | None:
    token = header.rstrip().split(maxsplit=1)[0]
    if regex:
        match = re.search(regex, token)
        if not match:
            return None
        if "barcode" in match.groupdict():
            return match.group("barcode")
        return match.group(1)
    pieces = token.split("_")
    try:
        return pieces[field]
    except IndexError:
        return None


def split_fastq(
    read1: str | Path,
    read2: str | Path,
    whitelist_path: str | Path,
    output_dir: str | Path,
    *,
    header_regex: str | None = None,
    header_field: int = -2,
    max_open_files: int = 128,
    gzip_output: bool = True,
) -> Counter[str]:
    """Split paired FASTQs by a barcode encoded in the read name.

    An LRU file-handle cache prevents the original implementation's
    ``too many open files`` failure for large cell counts.
    """
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    with Path(whitelist_path).open(encoding="utf-8") as handle:
        whitelist = {line.split("\t", 1)[0].strip() for line in handle if line.strip()}
    if not whitelist:
        raise ValueError("Barcode whitelist is empty.")

    suffix = ".fastq.gz" if gzip_output else ".fastq"
    handles: OrderedDict[str, tuple[TextIO, TextIO]] = OrderedDict()
    created: set[str] = set()
    counts: Counter[str] = Counter()

    def get_handles(barcode: str) -> tuple[TextIO, TextIO]:
        if barcode in handles:
            handles.move_to_end(barcode)
            return handles[barcode]
        if len(handles) >= max_open_files:
            _, old = handles.popitem(last=False)
            old[0].close()
            old[1].close()
        mode = "at" if barcode in created else "xt"
        pair = (
            open_text(output / f"{barcode}_R1{suffix}", mode),
            open_text(output / f"{barcode}_R2{suffix}", mode),
        )
        created.add(barcode)
        handles[barcode] = pair
        return pair

    try:
        with ExitStack() as stack:
            r1 = stack.enter_context(open_text(read1))
            r2 = stack.enter_context(open_text(read2))
            iterator1, iterator2 = records(r1), records(r2)
            while True:
                try:
                    item1 = next(iterator1)
                except StopIteration:
                    try:
                        next(iterator2)
                    except StopIteration:
                        break
                    raise ValueError("R2 contains more records than R1.")
                try:
                    item2 = next(iterator2)
                except StopIteration as error:
                    raise ValueError("R1 contains more records than R2.") from error
                barcode = barcode_from_header(item1[0], header_regex, header_field)
                if barcode in whitelist:
                    if not re.fullmatch(r"[A-Za-z0-9._-]+", barcode):
                        raise ValueError(f"Unsafe barcode for an output filename: {barcode!r}")
                    out1, out2 = get_handles(barcode)
                    out1.writelines(item1)
                    out2.writelines(item2)
                    counts[barcode] += 1
    finally:
        for pair in handles.values():
            pair[0].close()
            pair[1].close()
    return counts


def trim_prefix_if_matches(
    input_fastq: str | Path,
    output_fastq: str | Path,
    prefix: str,
    trim_length: int,
) -> tuple[int, int]:
    """Trim a configurable number of bases only when a read starts with *prefix*."""
    if trim_length < len(prefix):
        raise ValueError("trim_length cannot be shorter than the matching prefix")
    total = trimmed = 0
    destination = Path(output_fastq)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with open_text(input_fastq) as incoming, open_text(destination, "xt") as outgoing:
        for header, sequence, plus, quality in records(incoming):
            total += 1
            sequence_text, quality_text = sequence.rstrip("\r\n"), quality.rstrip("\r\n")
            if sequence_text.startswith(prefix):
                sequence_text, quality_text = sequence_text[trim_length:], quality_text[trim_length:]
                trimmed += 1
            outgoing.write(f"{header}{sequence_text}\n{plus}{quality_text}\n")
    return total, trimmed

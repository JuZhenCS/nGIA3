#!/usr/bin/env python3
"""Create an nGIA3 packed database from a FASTA file.

The implementation is intentionally independent from the C++ code in
src/makedb.cpp. It uses explicit little-endian, fixed-width fields so the
database layout is reproducible on every 64-bit platform supported by nGIA3.
NumPy is used automatically when available to accelerate MinHash generation;
the standard-library implementation is the correctness reference and fallback.
"""

from __future__ import annotations

import argparse
import hashlib
import os
from dataclasses import dataclass
from pathlib import Path
import struct
import sys
import tempfile
from typing import Iterable, Literal, Sequence as TypingSequence


UINT32_MASK = 0xFFFFFFFF
KMER_MASK = (1 << 25) - 1
SIGNATURE_COUNT = 128
MAX_SEQUENCE_LENGTH = 0xFFFF  # The C++ implementation uses length < 0xffff.
MAX_SEQUENCE_COUNT = 0x7FFFFFFF
HASH_SEEDS = tuple(
    (index * 0x9E3779B1 + 0x85EBCA6B) & UINT32_MASK
    for index in range(SIGNATURE_COUNT)
)

_U32 = struct.Struct("<I")
_U64 = struct.Struct("<Q")


class FastaFormatError(ValueError):
    """Raised when the input cannot be interpreted unambiguously as FASTA."""


@dataclass(frozen=True, slots=True)
class FastaRecord:
    """A normalized FASTA record stored as bytes for lossless output."""

    name: bytes
    sequence: bytes


def _build_translation_table() -> bytes:
    table = bytearray(256)
    groups = {
        1: b"aA",
        2: b"cC",
        3: b"dD",
        4: b"eEqQ",
        5: b"fFwWyY",
        6: b"gG",
        7: b"hH",
        8: b"iIvV",
        9: b"kK",
        10: b"lLmM",
        11: b"nN",
        12: b"pP",
        13: b"rR",
        14: b"sS",
        15: b"tT",
    }
    for code, residues in groups.items():
        for residue in residues:
            table[residue] = code
    return bytes(table)


TRANSLATION_TABLE = _build_translation_table()


def read_fasta(path: os.PathLike[str] | str) -> list[FastaRecord]:
    """Read and normalize FASTA records using the nGIA3 input contract.

    Sequence lines are concatenated without changing case or residue bytes.
    Records of 65535 residues or longer are omitted, matching the C++ code.
    Equal-length records retain their input order.
    """

    fasta_path = Path(path)
    records: list[FastaRecord] = []
    name: bytes | None = None
    sequence_parts: list[bytes] = []
    sequence_count = 0
    reached_limit = False

    def finish_record() -> None:
        if name is None:
            return
        sequence = b"".join(sequence_parts)
        if len(sequence) < MAX_SEQUENCE_LENGTH:
            records.append(FastaRecord(name, sequence))

    with fasta_path.open("rb") as source:
        for line_number, raw_line in enumerate(source, start=1):
            if not raw_line.endswith(b"\n"):
                raise FastaFormatError("FASTA file must end with a newline (LF)")
            line = raw_line[:-1]
            if line.endswith(b"\r"):
                line = line[:-1]
            if b"\r" in line:
                raise FastaFormatError("FASTA file contains a bare carriage return")
            if line.startswith(b">"):
                finish_record()
                if sequence_count >= MAX_SEQUENCE_COUNT:
                    reached_limit = True
                    break
                sequence_count += 1
                name = line
                sequence_parts = []
            else:
                if name is None:
                    raise FastaFormatError(
                        f"line {line_number}: sequence data appears before a FASTA header"
                    )
                sequence_parts.append(line)

    if not reached_limit:
        finish_record()

    records.sort(key=lambda record: len(record.sequence), reverse=True)
    return records


def _iter_kmers(sequence: bytes) -> Iterable[int]:
    """Yield the exact rolling k-mers signed by the C++ implementation."""

    kmer = 0
    length = len(sequence)
    for index, residue in enumerate(sequence):
        code = TRANSLATION_TABLE[residue]
        kmer = ((kmer << 4) + code) & KMER_MASK
        # For short sequences the C++ code signs only the final partial k-mer.
        if index < 6 and index < length - 1:
            continue
        yield kmer


def _wang_hash(value: int, seed: int) -> int:
    """Return the C++ uint32_t Wang-style hash, including every overflow."""

    result = (value ^ seed) & UINT32_MASK
    result = ((result ^ 61) ^ (result >> 16)) & UINT32_MASK
    result = (result + ((result << 3) & UINT32_MASK)) & UINT32_MASK
    result = (result ^ (result >> 4)) & UINT32_MASK
    result = (result * 0x27D4EB2D) & UINT32_MASK
    result = (result ^ (result >> 15)) & UINT32_MASK
    result = ((result ^ (result >> 8)) * 0x9E3779B1) & UINT32_MASK
    return (result ^ (result >> 14)) & UINT32_MASK


def signatures_python(sequence: bytes) -> tuple[int, ...]:
    """Compute 128 MinHash signatures using only the Python standard library."""

    minimum = [UINT32_MASK] * SIGNATURE_COUNT
    for kmer in _iter_kmers(sequence):
        for index, seed in enumerate(HASH_SEEDS):
            value = _wang_hash(kmer, seed)
            if value < minimum[index]:
                minimum[index] = value
    return tuple(minimum)


def signatures_numpy(sequence: bytes) -> tuple[int, ...]:
    """Compute signatures with NumPy's explicitly typed uint32 operations."""

    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise RuntimeError("the NumPy signature engine requires numpy") from exc

    kmers = np.fromiter(_iter_kmers(sequence), dtype=np.uint32)
    if kmers.size == 0:
        return (UINT32_MASK,) * SIGNATURE_COUNT

    with np.errstate(over="ignore"):
        indices = np.arange(SIGNATURE_COUNT, dtype=np.uint32)
        seeds = indices * np.uint32(0x9E3779B1) + np.uint32(0x85EBCA6B)
        minimum = np.full(SIGNATURE_COUNT, UINT32_MASK, dtype=np.uint32)
        # Bound peak memory for exceptionally long accepted records.
        for start in range(0, kmers.size, 4096):
            values = kmers[start : start + 4096, None] ^ seeds[None, :]
            values = (values ^ np.uint32(61)) ^ (values >> np.uint32(16))
            values = values + (values << np.uint32(3))
            values = values ^ (values >> np.uint32(4))
            values = values * np.uint32(0x27D4EB2D)
            values = values ^ (values >> np.uint32(15))
            values = (values ^ (values >> np.uint32(8))) * np.uint32(0x9E3779B1)
            values = values ^ (values >> np.uint32(14))
            minimum = np.minimum(minimum, values.min(axis=0))
    return tuple(int(value) for value in minimum)


SignatureEngine = Literal["auto", "python", "numpy"]


def make_signatures(sequence: bytes, engine: SignatureEngine = "auto") -> tuple[int, ...]:
    """Compute signatures with the selected deterministic implementation."""

    if engine == "python":
        return signatures_python(sequence)
    if engine == "numpy":
        return signatures_numpy(sequence)
    if engine != "auto":
        raise ValueError(f"unknown signature engine: {engine}")
    try:
        return signatures_numpy(sequence)
    except RuntimeError:
        return signatures_python(sequence)


def pack_sequence(sequence: bytes) -> tuple[int, ...]:
    """Encode residues into nGIA3's four 32-residue bit planes."""

    words = [0] * (1 + ((len(sequence) + 31) // 32) * 4)
    words[0] = len(sequence)
    for index, residue in enumerate(sequence):
        code = TRANSLATION_TABLE[residue]
        bit = 1 << (index % 32)
        base = 1 + (index // 32) * 4
        for plane in range(4):
            if code & (1 << plane):
                words[base + plane] |= bit
    return tuple(words)


def _write_u32_values(file_object, values: TypingSequence[int]) -> None:
    if values:
        file_object.write(struct.pack(f"<{len(values)}I", *values))


def _write_u64_values(file_object, values: TypingSequence[int]) -> None:
    if values:
        file_object.write(struct.pack(f"<{len(values)}Q", *values))


def _calculate_offsets(
    records: TypingSequence[FastaRecord],
) -> tuple[list[int], list[int], int]:
    count = len(records)
    offset = _U32.size + count * _U32.size * 2
    offset += count * _U64.size * 2
    offset += count * SIGNATURE_COUNT * _U32.size
    packed_offsets: list[int] = []
    for record in records:
        packed_offsets.append(offset)
        offset += (1 + ((len(record.sequence) + 31) // 32) * 4) * _U32.size

    fasta_offsets: list[int] = []
    for record in records:
        fasta_offsets.append(offset)
        offset += len(record.name) + len(record.sequence) + 2
    return packed_offsets, fasta_offsets, offset


def pack_database(
    fasta_path: os.PathLike[str] | str,
    packed_path: os.PathLike[str] | str,
    *,
    engine: SignatureEngine = "auto",
) -> int:
    """Create *packed_path* atomically and return the retained record count."""

    source = Path(fasta_path)
    destination = Path(packed_path)
    if destination.exists() and os.path.samefile(source, destination):
        raise FastaFormatError("input FASTA and packed output must be different files")

    records = read_fasta(source)
    packed_offsets, fasta_offsets, final_size = _calculate_offsets(records)
    destination.parent.mkdir(parents=True, exist_ok=True)

    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w+b",
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        ) as output:
            temporary_name = output.name
            output.write(_U32.pack(len(records)))
            _write_u32_values(output, [len(record.name) for record in records])
            _write_u32_values(output, [len(record.sequence) for record in records])
            _write_u64_values(output, packed_offsets)
            _write_u64_values(output, fasta_offsets)

            for record in records:
                _write_u32_values(output, make_signatures(record.sequence, engine))

            for expected_offset, record in zip(packed_offsets, records):
                if output.tell() != expected_offset:
                    raise AssertionError("internal packed-offset calculation mismatch")
                _write_u32_values(output, pack_sequence(record.sequence))

            for expected_offset, record in zip(fasta_offsets, records):
                if output.tell() != expected_offset:
                    raise AssertionError("internal FASTA-offset calculation mismatch")
                output.write(record.name + b"\n" + record.sequence + b"\n")

            if output.tell() != final_size:
                raise AssertionError("internal final-size calculation mismatch")
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary_name, destination)
        temporary_name = None
    finally:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass
    return len(records)


def md5sum(path: os.PathLike[str] | str) -> str:
    """Return a file's hexadecimal MD5 digest without loading it all at once."""

    digest = hashlib.md5()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_args(argv: TypingSequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--fasta", required=True, type=Path, help="input FASTA file")
    parser.add_argument("-p", "--packed", required=True, type=Path, help="output packed file")
    parser.add_argument(
        "--engine",
        choices=("auto", "python", "numpy"),
        default="auto",
        help="signature implementation (default: NumPy when installed, otherwise Python)",
    )
    return parser.parse_args(argv)


def main(argv: TypingSequence[str] | None = None) -> int:
    args = _parse_args(argv)
    try:
        count = pack_database(args.fasta, args.packed, engine=args.engine)
    except (OSError, FastaFormatError, RuntimeError) as exc:
        print(f"makedb: error: {exc}", file=sys.stderr)
        return 1
    print(f"packed {count} sequences into {args.packed}")
    print(f"md5: {md5sum(args.packed)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

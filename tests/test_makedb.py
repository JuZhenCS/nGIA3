from __future__ import annotations

from pathlib import Path
import subprocess
import struct
import sys
import tempfile
import unittest

from tools.makedb import (
    FastaFormatError,
    SIGNATURE_COUNT,
    make_signatures,
    pack_database,
    pack_sequence,
    read_fasta,
    signatures_python,
)


class MakeDbTests(unittest.TestCase):
    def test_multiline_input_is_stably_sorted_and_normalized(self) -> None:
        fasta = b">short\nAC\nDX\n>long\nACDEFGH\n>same\nacdefgh\n"
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "input.fas")
            source.write_bytes(fasta)
            records = read_fasta(source)

        self.assertEqual(
            [record.name for record in records],
            [b">long", b">same", b">short"],
        )
        self.assertEqual(
            [record.sequence for record in records],
            [b"ACDEFGH", b"acdefgh", b"ACDX"],
        )

    def test_pack_sequence_uses_four_bit_planes_and_zero_for_unknowns(self) -> None:
        # A=1, C=2, D=3, X=0. For positions 0..3 plane 0 is 0101b and
        # plane 1 is 0110b; the upper planes are clear.
        self.assertEqual(pack_sequence(b"ACDX"), (4, 0b0101, 0b0110, 0, 0))

    def test_database_offsets_and_payload_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "input.fas")
            output = Path(directory, "output.packed")
            source.write_bytes(b">a\nACDX\n>b\nACDEFGH\n")
            self.assertEqual(pack_database(source, output, engine="python"), 2)
            data = output.read_bytes()

        count = struct.unpack_from("<I", data, 0)[0]
        self.assertEqual(count, 2)
        packed_offsets = struct.unpack_from("<2Q", data, 4)
        fasta_offsets = struct.unpack_from("<2Q", data, 20)
        signature_end = 4 + count * 16 + count * SIGNATURE_COUNT * 4
        self.assertEqual(packed_offsets[0], signature_end)
        self.assertEqual(struct.unpack_from("<I", data, packed_offsets[0])[0], 7)
        self.assertEqual(struct.unpack_from("<I", data, packed_offsets[1])[0], 4)
        self.assertEqual(data[fasta_offsets[0] :], b">b\nACDEFGH\n>a\nACDX\n")

    def test_command_line_entry_point_creates_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "input.fas")
            output = Path(directory, "output.packed")
            source.write_bytes(b">a\nACDX\n")
            script = Path(__file__).parents[1] / "tools" / "makedb.py"
            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "-f",
                    str(source),
                    "-p",
                    str(output),
                    "--engine",
                    "python",
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(output.is_file())
            self.assertIn("md5:", result.stdout)

    def test_numpy_and_reference_signatures_match(self) -> None:
        try:
            actual = make_signatures(b"ACDEFGHIKLMNPQRSTVWYXX", "numpy")
        except RuntimeError:
            self.skipTest("NumPy is not installed")
        self.assertEqual(
            actual,
            signatures_python(b"ACDEFGHIKLMNPQRSTVWYXX"),
        )

    def test_empty_sequence_has_all_maximum_signatures(self) -> None:
        self.assertEqual(
            signatures_python(b""),
            (0xFFFFFFFF,) * SIGNATURE_COUNT,
        )

    def test_accepts_and_normalizes_crlf(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "input.fas")
            source.write_bytes(b">a\r\nAC\r\n")
            records = read_fasta(source)
        self.assertEqual(records[0].name, b">a")
        self.assertEqual(records[0].sequence, b"AC")

    def test_rejects_input_output_alias_without_data_loss(self) -> None:
        fixture = b">record\nACDEFGH\n"
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "input.fas")
            source.write_bytes(fixture)
            with self.assertRaises(FastaFormatError):
                pack_database(source, source, engine="python")
            self.assertEqual(source.read_bytes(), fixture)

    def test_rejects_ambiguous_input(self) -> None:
        invalid_inputs = (
            b">a\rAC\n",
            b">a\nAC",
            b"AC\n>a\nAC\n",
        )
        for content in invalid_inputs:
            with (
                self.subTest(content=content),
                tempfile.TemporaryDirectory() as directory,
            ):
                source = Path(directory, "input.fas")
                source.write_bytes(content)
                with self.assertRaises(FastaFormatError):
                    read_fasta(source)


if __name__ == "__main__":
    unittest.main()

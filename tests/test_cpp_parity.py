from __future__ import annotations

import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from tools.makedb import pack_database


class CppParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        binary_value = os.environ.get("NGIA3_BINARY")
        if not binary_value:
            raise unittest.SkipTest("NGIA3_BINARY is not set")
        cls.binary = Path(binary_value).resolve()
        if not cls.binary.is_file():
            raise unittest.SkipTest(f"nGIA3 binary does not exist: {cls.binary}")

    def run_cpp(
        self,
        source: Path,
        output: Path,
        *,
        threads: int | None = None,
    ) -> subprocess.CompletedProcess[str]:
        environment = os.environ.copy()
        if threads is not None:
            environment["OMP_NUM_THREADS"] = str(threads)
        return subprocess.run(
            [
                str(self.binary),
                "makedb",
                "-f",
                str(source),
                "-p",
                str(output),
            ],
            capture_output=True,
            text=True,
            check=False,
            env=environment,
        )

    def test_cpp_and_python_outputs_are_byte_identical(self) -> None:
        fixture = (
            b">long\r\nACDEFGHIK\r\n"
            b">same\r\nacdefghik\r\n"
            b">short\r\nACDX\r\n"
            b">empty\r\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "input.fas")
            cpp_output = Path(directory, "cpp.packed")
            python_output = Path(directory, "python.packed")
            source.write_bytes(fixture)

            result = self.run_cpp(source, cpp_output)
            self.assertEqual(result.returncode, 0, result.stderr)
            pack_database(source, python_output, engine="python")
            self.assertEqual(cpp_output.read_bytes(), python_output.read_bytes())

    def test_cpp_output_is_deterministic_across_thread_counts(self) -> None:
        fixture = b">a\nACDEFGHIKLMNPQRSTVWY\n>b\nacdx\n>empty\n"
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "input.fas")
            one_thread = Path(directory, "one.packed")
            four_threads = Path(directory, "four.packed")
            source.write_bytes(fixture)

            one_result = self.run_cpp(source, one_thread, threads=1)
            four_result = self.run_cpp(source, four_threads, threads=4)
            self.assertEqual(one_result.returncode, 0, one_result.stderr)
            self.assertEqual(four_result.returncode, 0, four_result.stderr)
            self.assertEqual(one_thread.read_bytes(), four_threads.read_bytes())

    def test_cpp_atomically_replaces_existing_output(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "input.fas")
            output = Path(directory, "output.packed")
            expected = Path(directory, "expected.packed")
            source.write_bytes(b">record\nACDEFGH\n")
            output.write_bytes(b"previous database")
            pack_database(source, expected, engine="python")

            result = self.run_cpp(source, output)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(output.read_bytes(), expected.read_bytes())
            self.assertEqual(list(Path(directory).glob("output.packed.tmp.*")), [])

    def test_cpp_cleans_temporary_output_if_commit_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "input.fas")
            output = Path(directory, "database")
            marker = Path(output, "keep.txt")
            source.write_bytes(b">record\nACDEFGH\n")
            output.mkdir()
            marker.write_text("keep", encoding="utf-8")

            result = self.run_cpp(source, output)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(marker.read_text(encoding="utf-8"), "keep")
            self.assertEqual(list(Path(directory).glob("database.tmp.*")), [])

    def test_cpp_preserves_existing_output_for_missing_final_newline(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "invalid.fas")
            output = Path(directory, "output.packed")
            source.write_bytes(b">record\nACDEFGH")
            output.write_bytes(b"previous database")

            result = self.run_cpp(source, output)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(output.read_bytes(), b"previous database")
            self.assertEqual(list(Path(directory).glob("output.packed.tmp.*")), [])

    def test_cpp_rejects_input_output_alias_without_data_loss(self) -> None:
        fixture = b">record\nACDEFGH\n"
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "input.fas")
            source.write_bytes(fixture)

            result = self.run_cpp(source, source)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(source.read_bytes(), fixture)

    def test_cpp_returns_nonzero_for_invalid_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory, "invalid.fas")
            output = Path(directory, "output.packed")
            source.write_bytes(b">record\rAC\n")

            result = self.run_cpp(source, output)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()

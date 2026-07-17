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

    def run_cpp(self, source: Path, output: Path) -> subprocess.CompletedProcess[str]:
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

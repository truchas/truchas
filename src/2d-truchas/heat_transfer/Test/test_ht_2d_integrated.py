#!/usr/bin/env python3

"""Small end-to-end tests for the JSON-driven ht_2d application."""

import pathlib
import re
import struct
import subprocess
import sys
import tempfile


def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} HT_2D JSON_INPUT CASE", file=sys.stderr)
        return 2

    executable = pathlib.Path(sys.argv[1]).resolve()
    input_file = pathlib.Path(sys.argv[2]).resolve()
    case = sys.argv[3]

    output_dir = pathlib.Path(tempfile.mkdtemp(prefix=f"ht_2d_{case}_"))
    try:
        result = subprocess.run(
            [str(executable), str(input_file)],
            cwd=output_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        (output_dir / "stdout.txt").write_text(result.stdout)

        if result.returncode != 0:
            print(result.stdout, end="")
            print(f"FAIL: ht_2d returned {result.returncode}")
            return 1

        log_file = output_dir / "run.log"
        if not log_file.exists():
            print("FAIL: ht_2d did not produce run.log")
            return 1
        log = log_file.read_text()

        if "unrecoverable integration failure" in log:
            print("FAIL: unrecoverable integration failure")
            return 1
        if re.search(r"\b(?:nan|inf)\b", log, re.IGNORECASE):
            print("FAIL: non-finite value reported in run.log")
            return 1

        completed = re.findall(r"Completed integration to T =\s*([0-9.E+-]+)", log)
        if not completed:
            print("FAIL: completion record not found in run.log")
            return 1
        final_time = float(completed[-1])
        if abs(final_time - 1.0e-2) > 1.0e-12:
            print(f"FAIL: final time {final_time:g} differs from 0.01")
            return 1

        xmf = (output_dir / "out.xmf").read_text()
        fields = re.findall(
            r'<Attribute Name="T".*?<DataItem[^>]*Dimensions="(\d+)"'
            r'[^>]*Seek="(\d+)"', xmf, re.DOTALL)
        if len(fields) < 2:
            print("FAIL: temperature fields not found in XDMF output")
            return 1
        count, offset = map(int, fields[-1])
        with (output_dir / "out.bin").open("rb") as binary:
            binary.seek(offset)
            temperature = struct.unpack(f"<{count}d", binary.read(8 * count))

        if case == "uniform":
            error = max(abs(value - 2.0) for value in temperature)
            if error > 1.0e-10:
                print(f"FAIL: uniform temperature error {error:g}")
                return 1
        elif case == "linear":
            expected = sorted(1.0 - (i + 0.5) / 8.0 for i in range(8) for _ in range(8))
            error = max(abs(a - b) for a, b in zip(sorted(temperature), expected))
            if error > 1.0e-10:
                print(f"FAIL: steady linear profile error {error:g}")
                return 1

        print(f"PASS: {case} reached final time {final_time:g}")
        print(f"      artifacts: {output_dir}")
        return 0

    finally:
        # Keep the run directory available for post-test diagnosis.
        pass


if __name__ == "__main__":
    sys.exit(main())

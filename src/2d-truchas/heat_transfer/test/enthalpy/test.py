#!/usr/bin/env python3

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from ht_2d_test_util import check_close, execute, finish, run_case


def test():
    executable = sys.argv[1]
    case_dir = Path(__file__).parent
    specific_heat = run_case(
        executable, case_dir / "input_specific_heat.json", expected_final_time=0.5
    )
    specific_enthalpy = run_case(
        executable, case_dir / "input_specific_enthalpy.json", expected_final_time=0.5
    )

    check_close(
        specific_heat.data.field(specific_heat.final_step, "T"),
        specific_enthalpy.data.field(specific_enthalpy.final_step, "T"),
        1.0e-12,
        "equivalent temperature",
    )
    check_close(
        specific_heat.data.field(specific_heat.final_step, "H"),
        specific_enthalpy.data.field(specific_enthalpy.final_step, "H"),
        1.0e-12,
        "equivalent enthalpy",
    )
    finish("single-phase enthalpy equivalence")


if __name__ == "__main__":
    sys.exit(execute(test))

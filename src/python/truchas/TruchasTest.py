#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import scipy.linalg as spla


def compare_max(field1, field2, tol, name, time):
    """Calculates the maximum absolute difference between two fields,
    and compares against an input tolerance. Returns 0 if error is below
    tolerance and 1 otherwise. Prints result of comparison to terminal.
    Fields may be a scalar value."""
    err = abs(field1 - field2).max()
    report_test(name, time, err, tol, "max")
    return 1 if not err <= tol else 0 # checks if err is nan


def compare_l2(field1, field2, tol, name, time):
    """Calculates the L2 norm between two fields,
    and compares against an input tolerance. Returns 0 if error is below
    tolerance and 1 otherwise. Prints result of comparison to terminal.
    Fields may be a scalar value."""
    err = spla.norm(abs(field1 - field2))
    report_test(name, time, err, tol, "l2")
    return 1 if not err <= tol else 0 # checks if err is nan


def compare_max_rel(field1, field2, tol, name, time):
    """Calculates the maximum relative absolute difference between two fields,
    and compares against an input tolerance. Returns 0 if error is below
    tolerance and 1 otherwise. Prints result of comparison to terminal.
    Fields may be a scalar value."""
    err = abs((field1 - field2) / field2).max()
    report_test(name, time, err, tol, "max rel")
    return 1 if not err <= tol else 0 # checks if err is nan


def report_test(var, time, err, tol, compare_type):
    status = "PASS" if err <= tol else "FAIL"
    print("{:s}: {:s} at t={:8.2e} {:s} error={:8.2e} (tol={:8.2e})"
          .format(status, var, time, compare_type, err, tol))


def report_summary(nfail):
    if nfail == 0:
        print("All tests PASS")
    else:
        print(nfail, " tests FAIL")

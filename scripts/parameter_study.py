#!/usr/bin/env python3

#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import numpy as np

import truchas


# read the maximum temperature in block 1
def max_temperature_str(input_parameters, output):
    sid = output.num_series()
    temperature = output.field(sid, "Z_TEMP")
    block1 = output.region(1)
    max_temperature = max(temperature[block1])
    print(f"Max temperature = {max_temperature:.2f} K.")
    return f"{max_temperature:.2f}"


# compute the average temperature in block 1
def avg_temperature_str(input_parameters, output):
    sid = output.num_series()
    temperature = output.field(sid, "Z_TEMP")
    cell_volume = output.volumes()
    block1 = output.region(1)
    avg_temp = np.average(temperature[block1], weights=cell_volume[block1])
    print(f"Avg temperature = {avg_temp:.2f} K.")
    return f"{avg_tmp:.2f}"


if __name__ == "__main__":
    # Store the results from each run in a database. If we wish to run this
    # parameter study again with additional points, each run will be cached so
    # we don't need to re-run data we've already calculated.
    nproc = 8
    tenv = truchas.TruchasEnvironment.default(overwrite_output=True)
    tdb = truchas.TruchasDatabase("parameter_study")
    tstudy = truchas.TruchasStudy(tenv, tdb, nproc)

    # Suppose the template file has multiple values that need to be inserted,
    # but we're only doing a parameter study on "heat_flux2".
    input_parameters = {"heat_flux1": 4e4,
                        "heat_flux2": 1e4,
                        }
    heat_flux2_points = [1e4, 2e4, 4e4]

    # The output will be a text file with the following columns: "heat_flux2",
    # "max_temperature", and "average_temperature". The latter two are custom
    # output metrics, tied to functions defined above.
    output_metrics = {"max_temperature": max_temperature_str,
                      "average_temperature": avg_temperature_str,
                      }

    tstudy.do_1d_parameter_study("template-heatup.inp",
                                 input_parameters, "heat_flux2", heat_flux2_points,
                                 output_metrics, "vary_heat_flux.txt")

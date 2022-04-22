#!/usr/bin/env python3

#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import scipy.optimize as spopt

import truchas

def final_top_temperature(heat_flux, tenv, nprocs):
    # Generate input deck
    parameters = {"heat_flux": heat_flux}
    tenv.generate_input_deck(parameters, "template-heatup.inp", "heatup.inp")

    # Run truchas with the new input
    stdout, output = tenv.truchas(nprocs, "block-heatup.inp")

    # Read the temperature probe
    temperature = output.probe_data("top-temperature.txt")

    # It is useful to print both the input and output for every iteration,
    # since each iteration is a full simulation and can be slow.
    print(f"heat_flux = {heat_flux:.4e}, temperature = {temperature[-1,1]:.4e}\n")

    # Return the last element (final time) of the second column (temperature column)
    return temperature[-1,1]


if __name__ == "__main__":
    nprocs = 8

    # We need to allow overwriting the output directory, since each iteration of
    # the optimizer will re-run Truchas.
    tenv = truchas.TruchasEnvironment.default(overwrite_output=True)

    # We'll wrap our function final_top_temperature in a lambda, so we can pass
    # in the TruchasEnvironment and nprocs variables. We'll also offset by our
    # target value, since we're using a root finder.
    fx = lambda x: final_top_temperature(x, tenv, nprocs) - 500

    # The Scipy.Optimize.brentq root finder is used. Here we're passing our
    # lambda, an interval for our heat flux of [2.5e4, 5e4], a tolerance of 50
    # K, and absolute and relative tolerances on the heat flux.
    heat_flux = spopt.brentq(fx, 2.5e4, 5e4, xtol=50, rtol=5e-3)
    print(heat_flux)

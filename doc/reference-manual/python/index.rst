``truchas`` Python Module
=========================

Truchas provides a Python package for scripting. This is useful for generating
custom restarts, performing automated optimization, and performing automated
parameter studies, or computing custom metrics.

.. contents:: Contents
   :local:

Examples
--------

Custom Restart
^^^^^^^^^^^^^^

Frequently it is necessary to generate custom mapped restarts. For example, if
we need to map data to a new mesh, where there is a new block which needs to be
initialized. Generally, we will need to initialize this block with a material
and a temperature. The below script shows how to read in a Truchas h5 file, map
it to a new mesh, and modify some field data before generating a restart file.

.. code-block:: python

    import numpy as np

    import truchas

    # Load truchas data from the preheat run, map it to a new mesh, and scale that mesh.
    tdata = truchas.TruchasMappedData("preheat_output/preheat.h5", "mesh2.gen", 0.001)

    # Select the last cycle in the preheat output. Further calls will read data from
    # this cycle, overwrite data at this cycle, and generate a restart with our
    # modifications.
    sid = tdata.num_series()

    # Set the VOF to [0, 1] (purely the final material in the input) in block 2
    tdata.assign_value_block(sid, "Z_TEMP", 1, f)
    tdata.assign_value_block(sid, "Z_TEMP", 3, f)

    # Compute the function f = 1 + 2*x at cell centers across the mesh.
    xc = tdata.centroids()
    f = 1 + 2 * xc[0]

    # Reassign the temperature to the function f, but only in blocks 1 and 3.
    tdata.assign_value_block(sid, "Z_TEMP", 1, xc[:,0])

    # Alternatively, we could read the whole field and modify it before reassigning.
    temp = tdata.field(sid, "Z_TEMP")
    graphite_region = tdata.region(1, 3)
    temp[graphite_region] = f[graphite_region]
    tdata.assign_field(sid, "Z_TEMP", temp)

    # Alternatively, we could use Numpy functions
    temp = tdata.field(sid, "Z_TEMP")
    temp = np.where(tdata.blockid() == 1 or tdata.blockid() == 3, f, temp)
    tdata.assign_field(sid, "Z_TEMP", temp)

    # Write the custom restart file
    tdata.write_restart("preheat.restart", sid)

:superscript:`Example Truchas Python script. This script reads a Truchas output file, generates a mapped restart onto a new mesh, and overwrites data for some variables on a specific block on that new mesh.`


Optimization
^^^^^^^^^^^^

Python scripting can be used to execute Truchas in a loop, for example to modify
some heat flux boundary condition until the temperature at some point in the
mesh fits an expected value. We can wrap this entire process into a Python
function, then use an optimization function from the popular Scipy package to
solve :math:`T(F) = 500` for the temperature :math:`T` at :math:`x = (0, 0, 0)`.

.. code-block::

    ...

    &PROBE
      coord = 0, 0, 0
      data = 'temperature'
      data_file = 'top-temperature.txt'
    /

    ...

    &THERMAL_BC
      name = 'heat-coil'
      face_set_ids = 14
      type = 'flux'
      flux = -{heat_flux:.4e} ! [W / m^2]
    /

    ...

:superscript:`Part of a Truchas template input, here named template-heatup.inp. The text {heat_flux:.4e} will be replaced by the Python-generated input deck below.`

.. code-block:: python

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


Parameter Study
^^^^^^^^^^^^^^^

Suppose we wish to perform the following parameter study: Vary the heat flux on
one boundary condition, and measure the maximum and average temperature in block
1. The following script will define the custom output metrics (max and avg
temperature on block 1), the parameter dictionary, and perform the parameter
study. This will produce a text file `vary_heat_flux.txt` which contains a table
with the input and output values.

This script will also generate a database which caches Truchas outputs. This way,
if we re-run the parameter study with more heat flux values later, we don't need
to re-run entire simulations which have already been computed.

.. code-block:: python

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


Components
----------

TruchasEnvironment
^^^^^^^^^^^^^^^^^^

.. automodule:: truchas.TruchasEnvironment
   :members:

TruchasData
^^^^^^^^^^^^^^^^^^

.. automodule:: truchas.TruchasData
   :members:

TruchasMappedData
^^^^^^^^^^^^^^^^^^

.. automodule:: truchas.TruchasMappedData
   :members:
   :show-inheritance:

TruchasStudy
^^^^^^^^^^^^^^^^^^

.. automodule:: truchas.TruchasStudy
   :members:

TruchasDatabase
^^^^^^^^^^^^^^^^^^

.. automodule:: truchas.TruchasDatabase
   :members:

TruchasTest
^^^^^^^^^^^

.. automodule:: truchas.TruchasTest
   :members:

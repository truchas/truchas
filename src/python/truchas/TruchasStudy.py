#!/usr/bin/env python3

#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import time
import os

from .TruchasDatabase import TruchasDatabase

class TruchasStudy:
    """
    A class for performing 1D Truchas parameter studies. This class expects a
    Truchas environment which will perform the tests (or generate input decks)
    and a Truchas Database which will cache the results.

    One will generally want to execute either do_1d_parameter_study, or
    generate_1d_parameter_study_inputs and register_1d_parameter_study_outputs
    and do_1d_parameter_study in that order.

    In the former case, the parameter study will be performed by sequentially
    running Truchas on the supplied metric. In the latter case, a collection
    of input files will be generated, which might be executed on a cluster
    via a SLURM array job. The results of these studies are then registered
    to the supplied database via the register_1d_parameter_study_outputs method.
    At this point, running do_1d_parameter_study will read outputs from the
    database of Truchas outputs instead of executing simulations.
    """

    def __init__(self, tenv, tdb, nprocs):
        self._tenv = tenv
        self._tdb = tdb
        self._nprocs = nprocs
        self._working_dir = "parameter_study_inputs"


    def do_1d_parameter_study(self, initial_parameters, variable, points,
                              output_metrics, output_filename):
        """Do a 1D parameter study of the given variable, across the given
        points. Output metrics are recorded to a text file. Input arguments are:

        **initial_parameters**: A dictionary of input parameters to supply to
        TruchasEnvironment.generate_input_deck. Keys are the python format
        variables present in the template input deck. Values are either values
        to replace those format strings, or functions. If a function, it is
        evaluated with a single input: the dict of parameters in the current
        iteration.

        **variable**: The name of a variable in the initial_parameters
        dictionary which is to be varied for the parameter study.

        **points**: The values of the variable.

        **output_metrics**: A dictionary of keys (names) and functions for
        output metrics to be printed to the output file. Functions are expected
        to accept two arguments: the input parameter dictionary and a
        TruchasData output set.

        **output_filename**: A string with a filename, where a table of inputs
        and outputs will be written.

        TODO: Some current hardwired shortcomings...
            - Templates are hardcoded to be template.inp and template.json.
            - Only a single json toolpath file is permitted, or none. The
              parameter key "toolpath_file" is disregarded when generating an
              input deck's unique identifier for the TruchasDatabase. This is
              because the toolpath file ought to be named with the same hash as
              the input file.
        """
        print(f"Performing a parameter study on {variable}...")
        parameters = initial_parameters.copy()

        # write the header from the input parameter names and output names
        line = "\t".join(parameters) + "\t"
        line += "\t".join(output_metrics)
        with open(output_filename, "w") as fh: fh.write("# " + line + "\n")

        for i, p in enumerate(points):
            parameters[variable] = p
            print(f"Running with {variable} = {parameters[variable]:.2e}")

            # Parameters may be functions of other non-function parameters. If
            # a function which depends on other data is needed, it should be
            # constructed with that data hardcoded, or via a lambda.
            for key, val in initial_parameters.items():
                if key != variable and callable(val):
                    parameters[key] = val(parameters)
            toolpath_file = "{}.json".format(TruchasDatabase.identifier(parameters))
            parameters["toolpath_file"] = toolpath_file

            if not self._tdb.exists(parameters):
                print("Generating input deck ... ", end="")
                input_file = "{}-{}.inp".format(os.path.splitext(output_filename)[0], i)
                self._tenv.generate_input_deck(parameters, "template.inp", input_file)
                if os.path.isfile("template.json"):
                    self._tenv.generate_input_deck(parameters, "template.json", toolpath_file)
                print("done.")

                print("Running Truchas ... ")
                elapsed = time.time()
                stdout, output = self._tenv.truchas(self._nprocs, input_file)
                elapsed = time.time() - elapsed
                print("done.")
                print("Elapsed {:.0f} seconds.".format(elapsed))

                self._tdb.donate(parameters, input_file, output.directory)

            print("Grabbing from database ... ", end="")
            output = self._tdb.truchas_output(parameters)
            print("done.")

            # construct the output line and write it
            line = "\t".join([f"{v:.3e}" for k,v in parameters.items()
                              if not k == "toolpath_file"]) + "\t" # dump the inputs
            line += "\t".join([f(parameters, output) for f in output_metrics.values()]) # outputs
            with open(output_filename, "a") as fh: fh.write(line + "\n")

            # print("Fetching max temperature history ... ", end="")
            # temp_fname = "{}-temperature-{}.txt".format(os.path.splitext(output_filename)[0], i)
            # with open(temp_fname, "w") as fh:
            #     fh.write("time\tmax_temperature\n")
            #     for sid in range(1,output.num_series()+1):
            #         t = output.time(sid)
            #         max_temperature = max(output.field(sid, "Z_TEMP"))
            #         fh.write(f"{t:.2f}\t{max_temperature:.2f}\n")
            # print("done.")
            print()


    def generate_1d_parameter_study_inputs(self, initial_parameters, variable, points):
        """Generate input decks across a given span of data points. This is
        intended for generating inputs to then be simulated on a cluster.
        """
        parameters = initial_parameters.copy()
        os.makedirs(self._working_dir, exist_ok=True)

        for p in points:
            parameters[variable] = p
            # Parameters may be functions of other non-function parameters. If
            # a function which depends on other data is needed, it should be
            # constructed with that data hardcoded, or via a lambda.
            for key, val in initial_parameters.items():
                if key != variable and callable(val):
                    parameters[key] = val(parameters)
            identifier = TruchasDatabase.identifier(parameters)
            parameters["toolpath_file"] = f"{identifier}.json"

            input_file = os.path.join(self._working_dir, f"{identifier}.inp")
            toolpath_file = os.path.join(self._working_dir, f"{identifier}.json")
            print(identifier, os.path.isfile(input_file), parameters)

            if not os.path.isfile(input_file):
                self._tenv.generate_input_deck(parameters, "template.inp", input_file)
                if os.path.isfile("template.json"):
                    self._tenv.generate_input_deck(parameters, "template.json", toolpath_file)
        print()


    def register_1d_parameter_study_outputs(self, initial_parameters, variable, points):
        """Register Truchas simulation outputs and input decks to a database."""
        parameters = initial_parameters.copy()

        for p in points:
            parameters[variable] = p
            # Parameters may be functions of other non-function parameters. If
            # a function which depends on other data is needed, it should be
            # constructed with that data hardcoded, or via a lambda.
            for key, val in initial_parameters.items():
                if key != variable and callable(val):
                    parameters[key] = val(parameters)
            identifier = TruchasDatabase.identifier(parameters)
            parameters["toolpath_file"] = f"{identifier}.json"

            input_file = os.path.join(self._working_dir, f"{identifier}.inp")
            toolpath_file = os.path.join(self._working_dir, f"{identifier}.json")
            output_directory = os.path.join(self._working_dir, f"{identifier}_output")

            if os.path.isfile(input_file):
                self._tdb.donate(parameters, input_file, output_directory)

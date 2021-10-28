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
                              output_metrics, output_filename,
                              replacement_parameters=None,
                              extra_outputs=None):
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

        **replacement_parameters** (optional): A function of one argument, the
        parameter dictionary at a given instance of the parameter study. If
        provided, it can be used to generate a completely new set of template
        replacements from the given dictionary of parameters.

        **extra_outputs** (optional): A function of four arguments, the
        parameter dictionary, the Truchas output, the study filename (intended
        for generating new filenames), and the parameters study iteration index.
        This is for generating new files out of each run, for instance recording
        the entire maximum temperature evolution of each run.

        NOTES:
            - Toolpath filename should be specified in the template.inp file
              with the parameter toolpath_file. This parameter is skipped when
              generating the instance unique identifier, so both the Truchas
              input file and toolpath file are renamed to a unique identifier
              independent of the toolpath filename.

        TODO: Some current hardwired shortcomings...
            - Templates are hardcoded to be template.inp and template.json.
            - Only a single json toolpath file is permitted, or none. The
              parameter key "toolpath_file" is disregarded when generating an
              input deck's unique identifier for the TruchasDatabase. This is
              because the toolpath file ought to be named with the same hash as
              the input file.

        """
        print(f"Performing a parameter study on {variable}...")

        # write the header from the input parameter names and output names
        line = "\t".join(initial_parameters) + "\t"
        line += "\t".join(output_metrics)
        with open(output_filename, "w") as fh: fh.write("# " + line + "\n")

        for i, p in enumerate(points):
            # Parameters may be functions of other non-function parameters. If
            # a function which depends on other data is needed, it should be
            # constructed with that data hardcoded, or via a lambda.
            parameters = initial_parameters.copy()
            parameters[variable] = p
            print(f"Running with {variable} = {parameters[variable]:.2e}")
            parameters = self._eval_parameter_functions(parameters)
            replacements = self._replacements(parameters, replacement_parameters)

            if not self._tdb.exists(replacements):
                print("Generating input deck ... ", end="", flush=True)
                identifier = TruchasDatabase.identifier(replacements)
                input_file = f"{identifier}.inp"
                toolpath_file = f"{identifier}.json"
                self._tenv.generate_input_deck(replacements, "template.inp", input_file)
                if os.path.isfile("template.json"):
                    self._tenv.generate_input_deck(replacements, "template.json", toolpath_file)
                print("done.")

                print("Running Truchas ... ")
                elapsed = time.time()
                stdout, output = self._tenv.truchas(self._nprocs, input_file)
                elapsed = time.time() - elapsed
                print("done.")
                print("Elapsed {:.0f} seconds.".format(elapsed))

                self._tdb.donate(replacements, input_file, output.directory)

            print("Grabbing from database ... ", end="", flush=True)
            output = self._tdb.truchas_output(replacements)
            print("done.")

            # construct the output line and write it
            line = "\t".join([f"{v:.3e}" for k,v in parameters.items()
                              if not k == "toolpath_file"]) + "\t" # dump the inputs
            line += "\t".join([f(parameters, output) for f in output_metrics.values()]) # outputs
            with open(output_filename, "a") as fh: fh.write(line + "\n")

            if callable(extra_outputs): extra_outputs(parameters, output, output_filename, i)
            print()


    def generate_1d_parameter_study_inputs(self, initial_parameters, variable, points,
                                           replacement_parameters=None):
        """Generate input decks across a given span of data points. This is
        intended for generating inputs to then be simulated on a cluster.
        """
        os.makedirs(self._working_dir, exist_ok=True)

        for p in points:
            # Parameters may be functions of other non-function parameters. If
            # a function which depends on other data is needed, it should be
            # constructed with that data hardcoded, or via a lambda.
            parameters = initial_parameters.copy()
            parameters[variable] = p
            parameters = self._eval_parameter_functions(parameters)
            replacements = self._replacements(parameters, replacement_parameters)

            identifier = TruchasDatabase.identifier(replacements)
            input_file = os.path.join(self._working_dir, f"{identifier}.inp")
            toolpath_file = os.path.join(self._working_dir, f"{identifier}.json")
            print(identifier, os.path.isfile(input_file), replacements)

            if not os.path.isfile(input_file):
                self._tenv.generate_input_deck(replacements, "template.inp", input_file)
                if os.path.isfile("template.json"):
                    self._tenv.generate_input_deck(replacements, "template.json", toolpath_file)
        print()


    def register_1d_parameter_study_outputs(self, initial_parameters, variable, points,
                                            replacement_parameters=None):
        """Register Truchas simulation outputs and input decks to a database."""

        for p in points:
            # Parameters may be functions of other non-function parameters. If
            # a function which depends on other data is needed, it should be
            # constructed with that data hardcoded, or via a lambda.
            parameters = initial_parameters.copy()
            parameters[variable] = p
            parameters = self._eval_parameter_functions(parameters)
            replacements = self._replacements(parameters, replacement_parameters)

            identifier = TruchasDatabase.identifier(replacements)
            input_file = os.path.join(self._working_dir, f"{identifier}.inp")
            toolpath_file = os.path.join(self._working_dir, f"{identifier}.json")
            output_directory = os.path.join(self._working_dir, f"{identifier}_output")

            if os.path.isfile(input_file):
                self._tdb.donate(replacements, input_file, output_directory)


    @staticmethod
    def _eval_parameter_functions(input_parameters):
        """Convert from a dictionary of input parameters, and the optional extra
        conversion function, to the dictionary of template replacements.
        """
        parameters = input_parameters.copy()
        for key, val in input_parameters.items():
            if callable(val):
                parameters[key] = val(input_parameters)
        return parameters


    @staticmethod
    def _replacements(input_parameters, replacement_parameters):
        """Convert from a dictionary of input parameters, and the optional extra
        conversion function, to the dictionary of template replacements.
        """
        replacements = input_parameters.copy()
        if callable(replacement_parameters):
            replacements = replacement_parameters(replacements)
        identifier = TruchasDatabase.identifier(replacements)
        replacements["toolpath_file"] = f"{identifier}.json"
        return replacements

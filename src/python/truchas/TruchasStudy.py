#!/usr/bin/env python3

#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import time
import os
import multiprocessing
import copy
import itertools

from .TruchasDatabase import TruchasDatabase


def _write_basic_restart(output, restart_filename):
    sid = output.num_series()
    output.write_restart(restart_filename, sid)


class TruchasStudy:
    """A class for performing Truchas parameter studies which vary in one
    parameter (1D). This class expects a Truchas environment which will perform
    the tests (or generate input decks) and a Truchas Database which will cache
    the results.

    One will generally want to execute either :func:`do_1d_parameter_study`, or
    :func:`generate_1d_parameter_study_inputs` and
    :func:`register_1d_parameter_study_outputs` and
    :func:`do_1d_parameter_study` in that order.

    In the former case, the parameter study will be performed by sequentially
    running Truchas on the supplied metric. In the latter case, a collection of
    input files will be generated, which might be executed on a cluster via a
    SLURM array job. The results of these studies are then registered to the
    supplied database via the :func:`register_1d_parameter_study_outputs`
    method. At this point, running :func:`do_1d_parameter_study` will read
    outputs from the database of Truchas outputs instead of executing
    simulations.

    To initialize:

    :param tenv: A :class:`TruchasEnvironment` for running the parameter study.
    :type tenv: :class:`TruchasEnvironment`

    :param tdb: A :class:`TruchasDatabase` where results will be cached.
    :type tdb: :class:`TruchasDatabase`

    :param nproc: Number of MPI ranks to use. When used with
        :func:`generate_1d_parameter_study_inputs`, this parameter is ignored.

    :type nproc: int

    :param njobs: Number of Truchas simulations to run in simultaneously.

    :type njobs: int, optional (default 1)

    :param restart_file: Filename of Truchas restart file. This restart will be
        used for all runs.

    :type restart_file: str, optional (default None)

    :param generate_only: Only generate inputs. Intended for subsequent runs on
        HPC.

    :type generate_only: bool, optional (default False)
    """

    def __init__(self, tenv, tdb, nprocs, njobs=1, restart_file=None, generate_only=False):
        self._tenv = tenv
        self._tdb = tdb
        self._nprocs = nprocs
        self._njobs = njobs
        self._restart_file = restart_file
        self._working_dir = "parameter_study_inputs"
        self._generate_only = generate_only


    def do_1d_parameter_study(self, template_input_file, initial_parameters, variable, points,
                              output_metrics, output_filename,
                              extra_outputs=None):
        """Do a 1D parameter study of the given variable, across the given
        points. Output metrics are recorded to a text file.

        :param initial_parameters: A dictionary of input parameters to supply to
            :func:`TruchasEnvironment.generate_input_deck`. Keys are the python
            format variables present in the template input deck. Values are
            either values to replace those format strings, or functions. If a
            function, it is evaluated with a single input: the dict of
            parameters in the current iteration. This allows values to be
            generated based on other values at this iteration. E.g., the "h2"
            key could hold a lambda which always returns double the value of the
            "h1" key.

        :type initial_parameters: dict with str keys

        :param variable: The name of a variable in the initial_parameters
            dictionary which is to be varied for the parameter study.

        :type variable: str

        :param points: The values of the variable to iterate over for the
            parameter study.

        :type points: list of values

        :param output_metrics: A dictionary of output metrics to be printed to
            the output file. Functions are expected to accept two arguments: the
            input parameter dictionary and a :class:`TruchasData` output set.
            These functions are to be written by the user to generate custom
            outputs. The result of each function must be a string.

        :type output_metrics: dict with str keys and function values.

        :param output_filename: A filename where the table of inputs and outputs
            will be written.

        :type output_filename: str

        :param extra_outputs: A function of four arguments, the parameter
            dictionary, the Truchas output, the study filename (intended for
            generating new filenames), and the parameter study iteration index.
            This is for custom operations to be performed at each run. For
            instance, generating unique new files each containing the maximum
            temperature evolution for a different run.

        :type extra_outputs: function, optional

        .. note::

            - Toolpath filename should be specified in the template.inp file
              with the parameter ``{toolpath_file}``. This parameter is skipped
              when generating the instance unique identifier, so both the
              Truchas input file and toolpath file are renamed to a unique
              identifier independent of the toolpath filename.

        .. warning::
            TODO Some current hardwired shortcomings...

            - Toolpath templates are hardcoded to be "template.json".

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

        # Generate list of replacements and new input files (those not already in the database)
        replacements = [self._instance_replacements(initial_parameters, variable, p)
                        for p in points]
        self.run_inputs(template_input_file, replacements)

        # Mine database for requested data
        print("Reading database ... ", end="", flush=True)
        for i, r in enumerate(replacements):
            output = self._tdb.truchas_output(r)

            # construct the output line and write it
            line = "\t".join([f"{v:.3e}" if isinstance(v, (float, int)) else f"{v}"
                              for k,v in r.items()
                              if not k == "toolpath_file"]) + "\t" # dump the inputs
            line += "\t".join([f(r, output) for f in output_metrics.values()]) # outputs
            with open(output_filename, "a") as fh: fh.write(line + "\n")

            if callable(extra_outputs): extra_outputs(r, output, output_filename, i)
        print("done.")


    def write_restarts(self, replacements, restart_filename_generator="{}.restart",
                       restart_writer=_write_basic_restart):
        """Generate restart files from the database.

        :param restart_writer: An optional function to write the restart. Can pass in a function
        that handles mapping or running from some sid other than the last.

        :type restart_writer: function of TruchasData and filename, optional (default
        _write_basic_restart). Must be global (not lambda) so Multiprocessing can see it.

        """
        assert all(self._tdb.exists(r) for r in replacements)

        restart_filenames = (restart_filename_generator.format(self._tdb.identifier(r))
                             for r in replacements)
        args = [(self, r,
                 restart_filename,
                 restart_writer)
                for r, restart_filename in zip(replacements, restart_filenames)
                if self._fresh_restart_needed(r, restart_filename)]

        print(f"Generating {len(args)} restart files ... ")
        with multiprocessing.Pool(self._njobs * self._nprocs) as pool:
            pool.starmap(_write_restart, args)
        print("done generating restart files.\n")


    def _fresh_restart_needed(self, r, restart_filename):
        """Check if a fresh restart file is needed for the given replacements and restart filename.
        If the given restart filename is not an existing file on the system, then a new restart is
        needed. If there already is such a file on the system, then check if it is older than the H5
        file corresponding to the given replacements r.

        """
        if not os.path.isfile(restart_filename):
            return True

        # TODO: Construct filename without opening HDF5 file.
        output_filename = self._tdb.truchas_output(r).filename
        output_write_time = os.path.getmtime(output_filename)
        restart_write_time = os.path.getmtime(restart_filename)
        return restart_write_time < output_write_time


    def run_inputs(self, template_input_file, replacements,
                   restart_filename_generator=None):
        # Generate list of replacements and new input files (those not already in the database)
        print("Generating input decks ... ", end="", flush=True)
        new_replacements = [r for r in replacements if not self._tdb.exists(r)]
        new_input_files = self.generate_inputs(template_input_file, new_replacements)
        print(f"done. {len(new_input_files)} new input files.\n")

        if not self._generate_only:
            # Run Truchas on all new inputs
            if restart_filename_generator is None: restart_filename_generator = self._restart_file
            args = [(self, i, r,
                     (restart_filename_generator.format(self._tdb.identifier(r)))
                     if restart_filename_generator is not None else None)
                    for i, r in zip(new_input_files, new_replacements)]
            with multiprocessing.Pool(self._njobs) as pool:
                pool.starmap(_run, args)
            print("Completed all simulations.\n")


    def run_if(self, template_input, replacements, restart=None):
        if self._tdb.exists(replacements):
            output = self._tdb.truchas_output(replacements)
            return output.stdout, output

        input_file = self.generate_input_deck(template_input, replacements)
        stdout, output = self._tenv.truchas(self._nprocs, input_file,
                                            restart_file=restart,
                                            overwrite_output=True)
        self._tdb.donate(replacements, input_file, output.directory)

        output = self._tdb.truchas_output(replacements)
        return output.stdout, output


    @staticmethod
    def replacements_from_variation(params_default, params_variation):
        """Given a default parameter set, and the a dictionary of variations for some/all of those
        parameters, return a single list of all unique replacements.

        """
        param_set = {}
        for v, x in params_variation.items():
            pvi = copy.deepcopy(params_default)
            pvi[v] = x
            psi = [inputs for inputs in _dict_product(pvi)]
            for p in psi:
                # frozenset appears to distinguish int from float,
                # while identifier doesn't, preventing collisions.
                identifier = TruchasDatabase.identifier(p) #frozenset(p.items())
                param_set[identifier] = copy.deepcopy(p)

        return param_set.values()


    def generate_inputs(self, template_input_file, replacements, use_working_dir=False):
        """Generate Truchas input decks from a given template file(s) and a list
        of dictionaries of replacements.

        .. tip::
            If you need to delete the input decks generated by this input without
            deleting the templates, carefully use a command like the following:

            .. code-block:: bash

                ls -1 *inp | grep -v template | xargs rm

        """
        return [self.generate_input_deck(template_input_file, r, use_working_dir)
                for r in replacements]


    def register_outputs(self, replacements):
        for r in replacements: self.register_study_output(r)


    def _instance_replacements(self, initial_parameters, variable, p):
        # Parameters may be functions of other non-function parameters. If
        # a function which depends on other data is needed, it should be
        # constructed with that data hardcoded, or via a lambda.
        parameters = initial_parameters.copy()
        parameters[variable] = p
        return self._replacements(parameters)


    def generate_input_deck(self, template_input_file, replacements, use_working_dir=False):
        """Generate a usable Truchas input deck from a given template file(s)
        and a dictionary of replacements.
        """
        working_dir = self._working_dir if use_working_dir else ""
        if use_working_dir: os.makedirs(self._working_dir, exist_ok=True)
        identifier = TruchasDatabase.identifier(replacements)
        input_file = os.path.join(working_dir, f"{identifier}.inp")
        toolpath_file = os.path.join(working_dir, f"{identifier}.json")
        self._tenv.generate_input_deck(replacements, template_input_file, input_file)
        if os.path.isfile("template.json"):
            self._tenv.generate_input_deck(replacements, "template.json", toolpath_file)
        return input_file


    def register_study_output(self, replacements, use_working_dir=False):
        working_dir = self._working_dir if use_working_dir else ""
        identifier = TruchasDatabase.identifier(replacements)
        input_file = os.path.join(working_dir, f"{identifier}.inp")
        toolpath_file = os.path.join(working_dir, f"{identifier}.json")
        output_directory = os.path.join(working_dir, f"{identifier}_output")
        if os.path.isfile(input_file) and os.path.isdir(output_directory) \
           and not self._tdb.exists(replacements):
            self._tdb.donate(replacements, input_file, output_directory)


    def _replacements_1d(self, template_input_file,
                         initial_parameters, variable, points):
        replacements = []
        for p in points:
            # Parameters may be functions of other non-function parameters. If
            # a function which depends on other data is needed, it should be
            # constructed with that data hardcoded, or via a lambda.
            parameters = initial_parameters.copy()
            parameters[variable] = p
            replacements.append(self._replacements(parameters))
        return replacements


    def generate_1d_parameter_study_inputs(self, template_input_file,
                                           initial_parameters, variable, points):
        """Generate input decks across a given span of data points. This is
        intended for generating inputs to then be simulated on a cluster via
        SLURM, not in a Python environment. These output files will be stored
        in the directory ``parameter_study_inputs``. With this approach, it is
        expected to:

        1. Use :func:`generate_1d_parameter_study_inputs` to generate Truchas
           input files in the ``parameter_study_inputs`` directory.

        2. Run Truchas independent of Python, e.g. via a SLURM script.

        3. Use :func:`register_1d_parameter_study_inputs` to identify the
           Truchas outputs and move them into a :class:`TruchasDatabase`.

        4. Use :func:`do_1d_parameter_study` to "run" the parameter study. This
           time, all of the outputs will be recognized as already-generated in
           the database, so only the postprocessing will be performed.

        See :func:`do_1d_parameter_study` for argument documentation.

        An example SLURM script is shown below for running each of these
        generated input files in a single array job.

        .. code-block:: bash

            #SBATCH -a 1-16  # Run 16 duplicates of this job in parallel

            cd <path-to-parameter-study-inputs>

            # exit if this run doesn't correspond to a file
            nfiles=$(ls -1 *.inp | wc -l)
            (( $SLURM_ARRAY_TASK_ID > $nfiles )) && exit

            # select a unique input file from an alphabatized list and my job array ID
            inputfile=$(ls -1 *.inp | sed -n "${SLURM_ARRAY_TASK_ID}p")
            ml truchas
            srun -c 1 truchas $inputfile

        """
        replacements = self._replacements_1d(template_input_file, initial_parameters, variable,
                                             points)
        self.generate_inputs(template_input_file, replacements)
        #print(identifier, os.path.isfile(input_file), replacements)
        #print()


    def register_1d_parameter_study_outputs(self, initial_parameters, variable, points):
        """Register Truchas simulation outputs and input decks to a database. If
        :func:`generate_1d_parameter_study_inputs` is used to generate input
        files in the directory ``parameter_study_inputs``, then Truchas is used
        externally to generate outputs for each of these input files, this
        function is used to register those outputs to a
        :class:`TruchasDatabase`.

        See :func:`do_1d_parameter_study` for argument documentation. See
        :func:`generate_1d_parameter_study_inputs` for the workflow description.
        """
        replacements = self._replacements_1d(template_input_file, initial_parameters, variable,
                                             points)
        self.register_outputs(replacements)


    @staticmethod
    def _replacements(input_parameters):
        """Convert from a dictionary of input parameters, to the dictionary of
        template replacements.
        """
        replacements = input_parameters.copy()
        for key, val in input_parameters.items():
            if callable(val):
                replacements[key] = val(input_parameters)
        identifier = TruchasDatabase.identifier(replacements)
        replacements["toolpath_file"] = f"{identifier}.json"
        return replacements


def _dict_product(d):
    """
    Yield one dict per combination of the values in d.
    """
    keys = d.keys()
    for combo in itertools.product(*d.values()):
        yield dict(zip(keys, combo))


### Multiprocessing variables & routines need to be global ###
_lock = multiprocessing.Lock()
def _run(self, input_file, replacements, restart_file):
    elapsed = time.time()
    stdout, output = self._tenv.truchas(self._nprocs, input_file, restart_file=restart_file)
    elapsed = time.time() - elapsed
    print(f"Finished {input_file}. Elapsed {elapsed:.0f} seconds.")
    _lock.acquire()
    try:
        self._tdb.donate(replacements, input_file, output.directory)
    finally:
        _lock.release()


def _write_restart(self, replacements, restart_filename, restart_writer):
    elapsed = time.time()
    output = self._tdb.truchas_output(replacements)
    restart_writer(output, restart_filename)
    elapsed = time.time() - elapsed
    print(f"    Finished {restart_filename}. Elapsed {elapsed:.0f} seconds.")

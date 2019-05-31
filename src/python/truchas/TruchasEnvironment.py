#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import sys
import os
import argparse
import json
import subprocess

from .TruchasData import TruchasData
try:
    from .TruchasConfigBuild import TruchasConfig
except ImportError:
    from .TruchasConfigInstall import TruchasConfig


class TruchasEnvironment:
    """
    Contains a configuration to run Truchas.
    """

    def __init__(self, mpiexec, truchas_executable, input_dir, write_restart_executable, \
                 working_dir="."):
        # read configuration from json file generated by cmake

        self._mpiexec = mpiexec
        self._truchas_executable = truchas_executable
        self._input_dir = input_dir
        self._write_restart_executable = write_restart_executable
        self._working_dir = working_dir

        # make the working directory if it doesn't exist
        if not os.path.isdir(self._working_dir):
            os.makedirs(os.path.abspath(self._working_dir))

        assert os.path.isfile(self._truchas_executable)
        assert os.path.isdir(self._input_dir)


    @classmethod
    def default(cls, input_dir=None):
        """Initialize using the CMake-generated default configuration"""
        if input_dir is None:
            # The input file is expected to be in the same folder as the test python script.
            input_dir = os.path.dirname(sys.argv[0])

        # Working directory is in the build path if the script is in the build path.
        # This is used for testing. Otherwise, the working directory is the current
        # directory.
        try:
            test_dir = os.path.relpath(input_dir, start=TruchasConfig.test_source_dir)
            working_dir = os.path.join(TruchasConfig.test_build_dir, test_dir)
        except AttributeError:
            working_dir = "."

        return cls(TruchasConfig.mpiexec, TruchasConfig.truchas_executable, input_dir, \
                   TruchasConfig.write_restart_executable, working_dir)


    @classmethod
    def from_config_file(cls, config_file, input_dir, working_dir="."):
        """Initialize from a configuration file, input directory,
        and an optional working directory."""
        assert os.path.isfile(config_file)

        # Configuration file is json.
        with open(config_file, 'r') as f:
            config = json.load(f)
            mpiexec = config["mpiexec"]
            truchas_executable = config["truchas-executable"]
            write_restart_executable = config["write-restart-executable"]

        return cls(mpiexec, truchas_executable, input_dir, working_dir)


    @classmethod
    def from_argv(cls):
        """Initialize from command line arguments."""
        # The input file is expected to be in the same folder
        # as the test python script.
        input_dir = os.path.dirname(sys.argv[0])

        parser = argparse.ArgumentParser()
        parser.add_argument('-c', "--config", type=str, required=True,
                            help='CMake-generated configuration file')
        parser.add_argument('-o', "--output", type=str, default=".",
                            help='Output directory')
        args = parser.parse_args()

        return cls.from_config_file(args.config, input_dir, args.output)


    def truchas(self, nprocs, input_file, restart_file=None, output_dir=None):
        """Runs truchas with specified number of MPI ranks, and returns the
        terminal output and the output data. Input name should not include
        the .inp extension."""

        # find the absolute path to the input file
        input_name = os.path.splitext(input_file)[0]
        input_file_abs = os.path.join(self._input_dir, input_file)
        assert os.path.isfile(input_file_abs)

        # build the command
        command = "{:s} -n {:d} {:s}" \
            .format(self._mpiexec, nprocs, self._truchas_executable)

        # The user may be running this python script from
        # the same directory as the Truchas working directory.
        # If so, no need to specify output directory to Truchas.
        if output_dir is None:
            output_dir_abs = os.path.join(self._working_dir, input_name + "_output")
            if not os.path.samefile(self._working_dir,"."):
                command += " -o:" + output_dir_abs
        else:
            output_dir_abs = os.path.join(self._working_dir, output_dir)
            command += " -o:" + output_dir_abs


        if restart_file is not None:
            restart_file_abs = os.path.join(self._input_dir, restart_file)
            assert os.path.isfile(restart_file_abs)
            command += " -r:" + restart_file_abs

        command += " " + input_file_abs

        # run truchas
        print(command)
        process = subprocess.run(command, shell=True, universal_newlines=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # WARN: Some cases, like false input file, cause Truchas to exit with
        #       an error message and a *zero* exit code (success). We need to
        #       either parse this or change Truchas's behavior.
        try:
            process.check_returncode()
        except:
            print("ERROR: Truchas returned a nonzero exit code. Printing stdout, stderr.")
            print(process.stdout)
            print(process.stderr)
            raise

        # read the output and return
        output = TruchasData(os.path.join(output_dir_abs, input_name + ".h5"))
        return process.stdout + process.stderr, output

    def write_restart(self, h5file, cycle_number, restart_file):
        """Write a restart file from the given H5 dump and cycle number to
        the given output file. Files expected to be relative to working directory."""

        command = "{:s} -n {:d} -o '{:s}' '{:s}'" \
            .format(self._write_restart_executable, cycle_number, restart_file, \
                    h5file)

        print(command)
        process = subprocess.run(command, shell=True, universal_newlines=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        try:
            process.check_returncode()
        except:
            print("ERROR: write-restart returned a nonzero exit code. Printing stdout, stderr.")
            print(process.stdout)
            print(process.stderr)
            raise

    def open_data(self, h5file):
        """Returns an output data object from an h5file, expected to be
        relative to the working directory."""
        return TruchasData(os.path.join(self._working_dir, h5file))

    def output(self, output_file):
        """Return an object for the output indicated by output_file, relative
        to the input directory. Primarily used for reading golden output."""
        return TruchasData(os.path.join(self._input_dir, output_file))
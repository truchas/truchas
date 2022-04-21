#!/usr/bin/env python3

#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

import os
import shutil
import json
import collections
import hashlib

import truchas

_TruchasRun = collections.namedtuple("_TruchasRun", ("parameters",
                                                   "input_file",
                                                   "dump_dir",
                                                   "dump_file",
                                                   ))


class TruchasDatabase:
    """Records Truchas outputs into a database for easy lookup. This class is
    used by :class:`TruchasStudy`, and for most user code will not need to be
    used directly except for initialization.

    :param directory: A directory where the database will be stored. Should be
        either an empty directory, nonexistent directory, or the path of an
        already-existing database.
    :type directory: str
    """

    def __init__(self, directory):
        self._directory = os.path.abspath(os.path.normpath(directory))
        self._database_filename = os.path.join(self._directory, "database.json")
        if not os.path.isdir(self._directory):
            os.makedirs(self._directory)


    def _read_database(self):
        """Read the database file if it exists, generate an empty database if
        it doesn't."""
        if not os.path.isfile(self._database_filename):
            database = {}
        else:
            with open(self._database_filename, "r") as fh:
                database = json.load(fh)
            database = {k: _TruchasRun(*v) for k, v in database.items()}
        return database


    def _write_database(self, database):
        with open(self._database_filename, "w") as fh:
            json.dump(database, fh)


    @staticmethod
    def identifier(parameters):
        """Generates a unique identifier for the given parameter dictionary.

        :param parameters: A dictionary of template replacements. Unlike
            :func:`TruchasStudy.do_1d_parameter_study`, this cannot contain
            functions (rather, it expects the functions have already been
            evaluated). This is so that a unique identifier can be produced.

        :type parameters: a dict with str keys and non-function values

        :return: Unique identifier hash
        :rtype: str
        """

        # Hashing is done by formatting the dictionary as a sorted JSON string,
        # and hashing the result by SHA1. A direct hash() call won't work, since
        # this seems to be seeded differently every time Python is started.

        # The toolpath_file parameter is special, and not used to create
        # the unique identifier. This is because we typically want to
        # name the toolpath file with the same identifier as the input
        # file, and that filename also has to be written to the input
        # file.
        parameters_ = parameters.copy()
        if "toolpath_file" in parameters_: del parameters_["toolpath_file"]

        json_str = json.dumps(parameters_, sort_keys=True)
        identifier = hashlib.sha1(json_str.encode()).hexdigest()
        return identifier


    def exists(self, parameters):
        """Checks if the given parameter combination is present in the database.

        :param parameters: A dictionary of template replacements. Unlike
            :func:`TruchasStudy.do_1d_parameter_study`, this cannot contain
            functions (rather, it expects the functions have already been
            evaluated). This is so that a unique identifier can be produced.

        :type parameters: a dict with str keys and non-function values

        :return: True if the Truchas results associated with ``parameters`` is
            already present in the database, False otherwise.
        :rtype: bool
        """
        database = self._read_database()
        identifier = self.identifier(parameters)
        return identifier in database


    def truchas_output(self, parameters):
        """Return the Truchas output corresponding to the given parameters. If
        such an output does not exist, return None.

        :param parameters: A dictionary of template replacements. Unlike
            :func:`TruchasStudy.do_1d_parameter_study`, this cannot contain
            functions (rather, it expects the functions have already been
            evaluated). This is so that a unique identifier can be produced.

        :type parameters: a dict with str keys and non-function values

        :return: Truchas output associated with ``parameters``
        :rtype: :class:`TruchasData`
        """
        database = self._read_database()
        identifier = self.identifier(parameters)
        if identifier in database:
            dump_file = os.path.join(self._directory, database[identifier].dump_file)
            output = truchas.TruchasData(dump_file)
        else:
            output = None
        return output


    def donate(self, parameters, infile, dump_dir):
        """Add the input file and Truchas output directory associated with a set
        of parameters to the database. Note this will overwrite any existing
        item associated with the parameter set in the database.

        :param parameters: A dictionary of template replacements. Unlike
            :func:`TruchasStudy.do_1d_parameter_study`, this cannot contain
            functions (rather, it expects the functions have already been
            evaluated). This is so that a unique identifier can be produced.

        :type parameters: a dict with str keys and non-function values

        :param infile: A Truchas input file
        :type infile: str

        :parameter dump_dir: A Truchas output dir associated with infile
        :type dump_dir: str
        """
        database = self._read_database()
        identifier = self.identifier(parameters)

        # Get the unique in-database names.
        # Note the h5 dump is always named after original input.
        my_infile = f"{identifier}.inp"
        my_dump_dir = f"{identifier}_output"
        my_dump_file = os.path.splitext(os.path.basename(infile))[0] + ".h5"
        my_dump_file = os.path.join(my_dump_dir, my_dump_file)

        # Move the original input file and output directory into the database.
        # Here we take ownership of the files, and move them to a new place on
        # disk. This also replaces any existing entry in the database which
        # already has these parameters.
        abs_infile = os.path.join(self._directory, my_infile)
        abs_dump_dir = os.path.join(self._directory, my_dump_dir)
        if os.path.isdir(abs_dump_dir): shutil.rmtree(abs_dump_dir)
        os.rename(dump_dir, abs_dump_dir)
        os.rename(infile, abs_infile)

        database[identifier] = _TruchasRun(parameters, my_infile, my_dump_dir, my_dump_file)
        self._write_database(database)


    def delete(self, parameters):
        """Delete the data associated with the given parameter set from the
        database. Note this will delete all associated data from disk.

        :param parameters: A dictionary of template replacements. Unlike
            :func:`TruchasStudy.do_1d_parameter_study`, this cannot contain
            functions (rather, it expects the functions have already been
            evaluated). This is so that a unique identifier can be produced.

        :type parameters: a dict with str keys and non-function values
        """
        database = self._read_database()
        identifier = self.identifier(parameters)
        if identifier in database:
            if os.path.isdir(database[identifier].dump_dir):
                shutil.rmtree(database[identifier].dump_dir)
            if os.path.isfile(database[identifier].input_file):
                os.remove(database[identifier].input_file)
            del database[identifier]
            self._write_database(database)


    def append_new_parameters(self, parameters):
        """Append a new set of parameters to the database. This requires
        rehashing all identifiers.

        :param parameters: A dictionary of template replacements, including any
            new keys that were not present in previous database entries. For
            accuracy, any new keys must contain the value that was used in old
            runs (where presumably they were hardcoded in the input file).

        :type parameters: a dict with str keys and non-function values

        """
        old_database = self._read_database()
        new_database = {}
        for old_id, old_values in old_database.items():
            # Any newly-listed items in parameters take the given default.
            new_parameters = parameters.copy()
            for p, v in old_values.parameters.items():
                new_parameters[p] = v

            new_id = self.identifier(new_parameters)
            new_input_file = os.path.join(self._directory, f"{new_id}.inp")
            new_dump_dir = os.path.join(self._directory, f"{new_id}_output")
            new_dump_file = os.path.join(new_dump_dir, os.path.basename(old_values.dump_file))
            new_database[new_id] = _TruchasRun(new_parameters, new_input_file,
                                               new_dump_dir, new_dump_file)

            old_abs_input_file = os.path.join(self._directory, old_values.input_file)
            new_abs_input_file = os.path.join(self._directory, new_input_file)
            old_abs_dump_dir = os.path.join(self._directory, old_values.dump_dir)
            new_abs_dump_dir = os.path.join(self._directory, new_dump_dir)
            os.rename(old_abs_input_file, new_abs_input_file)
            os.rename(old_abs_dump_dir, new_abs_dump_dir)

        self._write_database(new_database)


    def backup_database(self):
        """Create a backup file ``database-backup.json`` in the database
        directory. This backs up the database structure, not the Truchas inputs
        nor output data.
        """
        backup_filename = os.path.join(self._directory, "database-backup.json")
        shutil.copy2(self._database_filename, backup_filename)

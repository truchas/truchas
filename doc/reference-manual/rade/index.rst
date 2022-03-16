Radiation Enclosure Tools
=========================

Introduction
------------

Truchas includes several tools for working with radiation enclosure datasets.

- ``genre``: Used to compute view factors and generate radiation enclosure files.
- ``vizre``: Used to generate visualization files for a given enclosure.
- ``cmpre``: Used to compare two radiation enclosure files.


Invoking GENRE
--------------

Genre generates a radiation enclosure dataset.

:Serial: ``genre [-f] input_file encl_file``
:Parallel: ``mpirun -np N genre [-f] input_file encl_file``

``encl_file`` is the name of the output view factor radiation enclosure file to be written. If a file with this name already exists, it will not be overwritten unless the ``-f`` flag is provided. The input file is a :ref:`Fortran namelist file <introduction/index:Input File Format>` containing the following namelists:

.. toctree::
   :maxdepth: 1

   ENCLOSURE Namelist <ENCLOSURE_Namelist>
   CHAPARRAL Namelist <CHAPARRAL_Namelist>
   PATCHES Namelist <PATCHES_Namelist>

.. note::
   - The order of the namelists does not matter, and there may be other stuff in the input file; behavior is similar to Truchas input files.
   - There must be a single ENCLOSURE namelist, at most one CHAPARRAL namelist, and at most one PATCHES namelist. If there is a CHAPARRAL namelist, the specified enclosure surface is generated and written to the enclosure file along with the calculated view factors. If there is *no* CHAPARRAL namelist, just the enclosure surface is written; this is useful for examining the surface for correctness prior to performing the expensive view factor calculation.
   - Displacing surfaces is optional. Any displaced surfaces must be totally disconnected from the remaining enclosure surfaces (this is checked). The constant displacement must be compatibile with the specified symmetries (this is not checked). This is useful when the mesh consists of disconnected parts that will be shifted relative to one another. Instead of generating multiple 3D meshes, a single 3D mesh plus use of this displacement option is all that is needed.


Invoking VIZRE
--------------

Writes a GMV-format visualization file for the specified enclosure. If the enclosure includes view factor data, the ambient view factor and view factor matrix row sums are written as face variables.

:Serial: ``vizre [options] enclosure_file gmv_file``
:Parallel: ``mpirun -np N vizre [options] enclosure_file gmv_file`` (only useful if the matrix is too large to be held on a single process.)

Vizre Command Line Options:

-r list     Write the specified *rows* of the view factor matrix as face variables. List is a comma-separated list of ranges; a range is index, first:last, or first:last:stride.
-c list     Write the specified *columns* of the view factor matrix as face variables. List is a comma-separated list of ranges; a range is index, first:last, or first:last:stride.
-s          Write the fully-developed enclosure surface defined by the enclosure's symmetries. The default is to write just the generating surface.
-h, --help  Display this help and exit.

.. note::
   - Only surface faces are written.
   - Faces are oriented and GMV can optionally darken the "bottom" side (option "test normals")
   - Surface side set info is written as GMV materials.
   - When using the ``-s`` option, info about the copies of the generating surface are written as GMV flags.
   - May need a parallel version simply because the VF matrix is too large to be held on a single process.


Invoking CMPRE
--------------

If you've generated multiple datasets for the same enclosure using different sets of Chaparral parameters, it is useful to compare the calculated view factor matrices. You can do this with cmpre:

:Serial: ``cmpre encl_file1 encl_file2``
:Parallel: ``mpirun -np N cmpre encl_file1 encl_file2`` (only useful if the matrix is too large to be held on a single process.)

The program simply prints out several VF operator norms of the difference of the two view factor matrices to the screen. It also prints the largest single difference in matrix and its row and column location. Naturally both files must contain VF data, and the enclosure surfaces must be the same.

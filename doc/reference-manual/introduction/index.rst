.. _um_introduction:

.. toctree::
   :maxdepth: 1

Introduction
=============

Invoking Truchas
-----------------
Truchas is executed in serial using a command of the form

   truchas [-h] [-v:`n`] [-o:`outdir`] [-r:`rstfile`] [-g:`m`] `infile`

assuming "truchas" is the name of the executable. The brackets denote optional
arguments that are described in :numref:`Table %s <truchas_command_line_option>`.
The only required argument is the path to the input file `infile`. This file
name must end with the extension ".inp". The general format of the input file
is described in the next section, and the following chapters describe the
various Fortran namelists that go into the input file to describe the problem
to be simulated.

All of the output files are written to directory whose name is derived from the
base name of the input file. For example, if the input file is "myprob.inp",
the output directory will be named "myprob_output". The name of the output
directory can be overridden using the -o option. The directory will be created
if necessary.

The precise manner of executing Truchas in parallel depends on the MPI
implementation being used. This may be as simple as prefixing the serial
invocation above with "mpirun -np n", where n is the number of processes.
But this varies widely and providing specific instructions is beyond the
scope of this document. There is no difference in the Truchas arguments
between serial and parallel, however.

.. _truchas_command_line_option:
.. csv-table:: Truchas Command Line Options
   :class: tight-table
   :widths: 1 5

   "Option", "Description"
   -g:`m`, "Sets the number `m` of MPI ranks used for parallel output to the HDF5 output file. The default is 1, and must be no greater than the number of MPI ranks."
   "-h", "Print a usage summary of the command line options and exit."
   "-f", "Force overwrite of output directory contents."
   "-m", "Turn on memory diagnostics. Writes to a .mem file in the output directory."
   "-o:`outdir`", "Causes all output files to be written to the directory `outdir` instead of the default directory. The directory is created if necessary."
   "-r:`rstfile`", "Executes in restart mode, restarting from the data in the file `rstfile`. This file is generated from the output of a previous Truchas simulation using post-processing utilities."
   "-v:`n`", "Sets the verbosity level of terminal output to `n`. The default level is 1. Level 0 produces no output, and level 2 produces somewhat more output than default."

Stopping Truchas
-----------------
There are occasions where one would like to gracefully terminate a running
Truchas simulation before it has reached the final simulation time given in
the input file. This is easily done by sending the running process the
``SIGUSR2`` signal.

::

    kill -s SIGUSR2 pid

where ``pid`` is the process id. When Truchas receives this signal, it
continues until it reaches the end of the current time step, where it writes
the final solution and then exits normally.

Input File Format
------------------
The Truchas input file is composed of a sequence of Fortran namelist inputs.
Each namelist input has the form

::

   &namelist-group-name
      namelist-input-body
   /

The input begins with a line where the first nonblank is the "&" character
immediately followed by the name of the namelist group. The input continues
until the "/" character. The body of the namelist input consists of a
sequence of "`name = value`" pairs, separated by commas or newlines, that
assign values to namelist variables. Namelist input is a feature of Fortran
and a complete description of the syntax can be found in any Fortran reference,
however the basic syntax is very intuitive and a few examples like those in
the following input sample should suffice to explain its essentials.

.. code-block::

   This is a comment.  Anything outside a namelist input is ignored.

   &MESH
      ! Within a namelist input "!" introduces a comment.
      mesh_file = "my-big-mesh.exo" ! character string value
   /

   Another comment.

   &PHYSICS
      heat_conduction = .true.     ! logical values are .true./.false.
      body_force = 0.0, 0.0, -9.8  ! assigning values to an array
      !This would be an equivalent method ...
      !body_force(1) =  0.0
      !body_force(2) =  0.0
      !body_force(3) = -9.8
   /

   Newlines in a namelist are optional.

   &PHYSICAL_CONSTANTS
      stefan_boltzmann=0.1, absolute_zero=0.0
   /

The namelists and the variables they contain are described in the following
chapters. A particular namelist will be required or optional, and it may
appear only once or multiple times. Not all variables of a namelist need
to be specified. Some may not be relevant in a given context, and others
may have acceptable default values; only those that are used and need to
be assigned a value need to be specified.

The order of the namelist inputs in the input file is not significant; they
may appear in any order. Any text outside of a namelist input is ignored
and can be regarded as a comment as illustrated in the example.

Fortran is case-insensitive when interpreting the namelist group names and
variable names; they may be written in any mixture of upper and lower case.
However character string *values*, which are interpreted by Truchas, are
case-sensitive unless documented otherwise. In the event of a namelist input
syntax error, Truchas will report that it was unable to read the namelist
and include what limited information it gets from the namelist parser.
Common syntax errors are misspelt variable names, variables that don't
belong to the namelist, a blank written in place of an underscore, etc.


.. _physcial_units:

Physical Units
---------------

Truchas does not require the use any particular system of physical units, nor
does it provide a means for the user to specify the dimension of a numerical
value. The user is simply expected to ensure that all input values are given
using a consistent system of units. To assist in this end, the dimension for
all dimensional quantities is documented using the following abstract units:
mass :math:`M`, length :math:`L`, time :math:`T`, thermodynamic temperature
:math:`\Theta`, and electric current :math:`I`. Thus mass density, for
example, will be documented as having dimension :math:`M/L^3`. The following
derived abstract units are also used: force :math:`F (= M L/T^2)` and energy
:math:`E (= M L^2/T^2)`.

There are a few physical constants, like the Stefan-Boltzmann constant, that
have predefined values in SI units. These constants are referenced by a few
specific models, and where a model uses one of these constants this fact is
noted. Use the :ref:`PHYSICAL_CONSTANTS <PHYSICAL_CONSTANTS_Namelist>`
namelist to redefine the value of these constants where necessary.

.. _UM_Work_With_Output:

Working With Output Files
--------------------------
As described earlier, Truchas writes its output files to the directory named
in the -o option, or if omitted, to a directory whose name is generated
from the base name of the input file: "myprob_output" if "myprob.inp" is
the input file, for example. Two primary files are written, a .log file
that is a copy of the terminal output, and a .h5 HDF5 file that contains
the simulation results. HDF5 is a widely-used format for storing and
managing data, and thereare a great many freely-available tools for working
with these files. Various additional human-readable text files may also be
written when requested by the input.

Visualization
^^^^^^^^^^^^^^^^^^^^^
`Paraview <https://www.paraview.org>`_ is the recommended tool for visualizing
the simulation results written to the .h5 file. After opening the file with
paraview, a pop-up window will appear prompting the user to select a reader;
choose the "TRUCHAS Reader" selection.

The legacy visualization software GMV (General Mesh Viewer) may also be used.
Originally developed by LANL, GMV was licensed to CPFD Software, LLC who
assumed its maintenance and commercialized it. Since then they have released
the `source <https://github.com/CPFDSoftware/gmv>`_ under the GPLv3 license on
GitHub. Documentation for GMV may be found at http://www.gmv-barracuda.com/index.html.
Note that support for GMV is deprecated and may be dropped in the future.
To use GMV, the program write-gmv.py must first be used to create GMV-format
input files from the .h5 file. The command syntax is

   write-gmv.py [`options`] `H5FILE`

where `H5FILE` is the .h5 output file. Use the option -h to get a full list
of the available options. This will create a collection of files: one for the
mesh and a numbered sequence of files for the solution fields at the output
times.

Restarts
^^^^^^^^^^^^^^^^^^^^^
The program write-restart.py is used to create Truchas restart files using
data from an .h5 output file. The command syntax is

   write-restart.py [`options`] `H5FILE`

where `H5FILE` is the .h5 output file and the possible `options` are:

   -h, --help    Display usage help and exit.
   -l            Print a list of the available cycles from which the restart
                 file can be created. No restart file is written.
   -n N          Data from cycle `N` is used to create the restart file;
                 if not specified, the last cycle is used.
   -o FILE       Write restart data to `FILE`. If not specified, `FILE` is
                 taken to be the `H5FILE` name with its ".h5" suffix replaced
                 by ".restart.N" where N is the cycle number.
   -m FILE       Create a mapped restart file using the specified ExodusII
                 mesh `FILE` as the target mesh.
   -s FLOAT      Scale the mapped restart mesh by the factor `FLOAT`.
   --use-portage   Use the Portage grid mapper backend. Truchas must have
                   been built with Portage support. The default is Truchas'
                   built-in Kuprat mapper.

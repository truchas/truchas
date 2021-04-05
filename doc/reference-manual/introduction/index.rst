.. _um_introduction:

.. toctree::
   :maxdepth: 1

Introduction
=============

Invoking TRUCHAS
-----------------
Truchas is executed in serial using a command of the form

``truchas [-h] [-d[:n]] [-o:outdir] [-r:rstfile] infile``

assuming ``truchas`` is the name of the executable. The brackets denote optional arguments that are described in :numref:`Table %s <truchas_command_line_option>`

.. _truchas_command_line_option:
.. csv-table:: Truchas Command Line Options
   :header: "Option", "Description"
   :class: tight-table
   :widths: 1 1

   "-h", "Print a usage summary of the command line options and exit."
   "-d[:n]", "Sets the debug output level **n**. The default level is 0, which produces no debug output, with levels 1 and 2 producing progressively more debug output. **-d** is equivalent to **-d:1**."
   "-o:outdir", "Causes all output files to be written to the directory **outdir** instead of the default directory. The directory is created if necessary."
   "-r:rstfile", "Executes in restart mode, restarting from the data in the file **rstfile**. This file is generated from the output of a previous Truchas simulation using post-processing utilities."

The only required argument is the path to the input file ``infile``. This file name must end with the extension **“.inp”** (without the quotes). The general format of the input file is described in the next section, and the following chapters describe the various Fortran namelists that go into the input file to describe the problem tobe simulated.
All of the output files are written to directory whose name is generated from the base name of the input file. For example, if the input file is ``myprob.inp``, the output directory will be named ``myprob_output``. The name of the output directory can be overridden using the **-o** option. The directory will be created if necessary.

The precise manner of executing Truchas in parallel depends on the MPI implementation being used. This may be as simple as prefixing the serial invocation above with ``mpirun -np n``, where **n** is the number of processes. But this varies widely and providing specific instructions is beyond the scope of this document. There is no difference in the Truchas arguments between serial and parallel, however.

.. _example_input:

::
   
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

Stopping Truchas
-----------------
There are occasions where one would like to gracefully terminate a running Truchas simulation before it has reached the final simulation time given in the input file. This is easily done by sending the running process the **SIGURG** signal.

``kill -s SIGURG pid``
where **pid** is the process id. When Truchas receives this signal, it continues until it reaches the end of the current time step, where it writes the final solution and then exits normally.

Input File Format
------------------
The Truchas input file is composed of a sequence of Fortran namelist inputs. Each namelist input has the form

::
   
   &namelist-group-name
      namelist-input-body
   /

The input begins with a line where the first nonblank is the character ``&`` immediately followed by the name of the namelist group. The input continues until the ``/`` character. The body of the namelist input consists of a sequence of ``name = value`` pairs, separated by
commas or newlines, that assign values to namelist variables. Namelist input is a feature of Fortran and a complete description of the syntax can be found in any Fortran reference, however the basic syntax is very intuitive and a few examples like those in :ref:`the example <example_input>` should suffice to explain its essentials. 

The namelists and the variables they contain are described in the following chapters. A particular namelist will be required or optional, and it may appear only once or multiple times. Not all variables of a namelist need to be specified. Some may not be relevant in agiven context, and others may have acceptable default values; only those that are used andneed to be assigned a value need to be specified.

The order of the namelist inputs in the input file is not significant; they may appear in any order. Any text outside of a namelist input is ignored and can be regarded as a comment; see :ref:`the example <example_input>`.

Fortran is case-insensitive when interpreting the namelist group names and variable names; they may be written in any mixture of upper and lower case. However character string values, which are interpreted by Truchas, are case-sensitive unless documented otherwise. In the event of a namelist **input syntax error**, Truchas will report that it was unable to read the namelist, but unfortunately it is not able to provide any specific information about the error because Fortran does not make such information available. In such cases the user will need to scan the namelist input for syntax errors: look for misspelt variable names, variables that don’t belong to the namelist, blank written in place of an underscore, etc.

.. _physcial_units:

Physical Units
---------------

Truchas does not require the use any particular system of physical units, nor does it provide a means for the user to specify the dimension of a numerical value. The user issimply expected to ensure that all input values are given using a consistent system of units. To assist in this end, the dimension for all dimensional quantities is documented using the following abstract units: mass **M**, length **L**, time **T**, thermodynamic temperature **Θ**, and electric current **I**. Thus mass density, for example, will be documented as having dimension :math:`M/L^3`. The following derived abstract units are also used: force :math:`F (= M*L/T^2)` and energy :math:`E (= M*L^2/T^2)`.

There are a few physical constants, like the **Stefan-Boltzmann** constant, that have predefined values in SI units. These constants are referenced by a few specific models, and where a model uses one of these constants this fact is noted. Use the :ref:`PHYSICAL_CONSTANTS <PHYSICAL_CONSTANTS_Namelist>` namelist to redefine the value of these constants where necessary.

.. _UM_Work_With_Output:

Working With Output Files
--------------------------
As described earlier, Truchas writes its output files to the directory named in the **-o** option, or if omitted, to a directory whose name is generated from the base name of the input file: **myprob_output** if **myprob.inp** is the input file, for example. Two primary files are written, a ``.log`` file that is a copy of the terminal output, and a ``.h5`` HDF5 file that contains all the simulation results. HDF5 is a widely-used format for storing and managing data, and thereare a great many freely-available tools for working with these files. In this release, which is the first to feature HDF5 output, we provide only a few essential tools, described below, for processing the ``.h5`` file. We expect to provide additional tools in future releases.

write_restart
^^^^^^^^^^^^^^
The program ``write_restart`` is used to create Truchas restart files using data from an ``.h5`` output file. The command syntax is

``write_restart [options] H5FILE``

where **H5FILE** is the **.h5** output file and the possible options are

.. _write_restart_options:
.. csv-table:: write_restart Options
   :header: "Option", "Description"
   :widths: 1 1
   :class: tight-table

   "**-h** ", "Display usage help and exit."
   "**-l** ", "Print a list of the available cycles from which the restart file can becreated. No restart file is written."
   "**-n N**", "Data from cycle **N** is used to create the restart file; if not specified, the last cycle is used."
   "**-o FILE**", "Write restart data to **FILE**. If not specified, **FILE** is taken to be the H5FILE name with the.h5 suffix replaced by .restart.N where N is the cycle number."
   "**-m FILE**", "Create a mapped restart file using the specified ExodusII mesh **FILE** as the target mesh."
   "**-s FLOAT**", "Scale the mapped restart mesh by the factor FLOAT."

write_probes
^^^^^^^^^^^^^^
The ``write_probes`` utility extracts probe data (see the :ref:`PROBE <PROBE_Namelist>` namelist) from an **.h5** output file and writes it to the terminal (where it can be redirected as needed) in a multicolumn format suitable for many line plotting programs. The command syntax is

``write_probes { -h | -l | -n N} H5FILE``

where **H5FILE** is the **.h5** output file and the available options are

.. _write_probes_options:
.. csv-table:: write_probes Options
   :header: "Option", "Description"
   :class: longtable
   :widths: 1 1

   "**-h**", "Display usage help and exit."
   "**-l**", "Print a list of the available probes."
   "**-n N**", "Data for probe index N is written."

truchas-gmv-parser.py
^^^^^^^^^^^^^^^^^^^^^^^
The python scriptt ``ruchas-gmv-parser.py`` is used to create input files for the GMV visualization tool. Formerly distributed gratis, GMV has been commercialized (http://www.generalmeshviewer.com). Earlier free versions of the tool can still be found on the internet, however, and it remains available within LANL. The command syntax is

``python truchas-gmv-parser.py [options] H5FILE``

where **H5FILE** is the **.h5** output file. Use the option **-h** to get a full list of the available options.
.. _input_file_description:

Input File
===========
.. sectionauthor:: Narendran Raghavan <naren@lanl.gov>

In this section input file for TRUCHAS will be explained. Input file should have an extension of **inp** ``file_name.inp`` . Example input file is provided :download:`here <input_file_example.inp>`. 

In this section, the input file format will be briefly explained. For detailed information regarding the input file format, please refer to the :ref:`Reference Manual <reference_manual>`. All lines in the input file are case-insensitive.

Input file has multiple namelist blocks starting with ``&`` and ending with ``\`` . Anything outside the block is considered as a comment. Inside the namelist block, texts preceded by ``!`` are considered a comment. 

The OUTPUTS namelist defines the problem end time and various output options.
::
   
   This is a comment

   &OUTPUTS
   output_t = 0.0, 600.0, 3600.0  ! This is a comment too
   output_dt = 60.0, 300.0
   /

The MESH namelist specifies the common mesh used by all physics models other than the induction heating model.
::
   
   &MESH
   mesh_file = 'mesh_file_name.gen'
   /  

The PHYSICS namelist specifies which physics models are active in the simulation.
::
   
   &PHYSICS
   heat_transport = T
   materials = 'graphite'
   /

The THERMAL_BC namelist is used to define boundary conditions for the heat transfer physics.
::
   
   &THERMAL_BC
   name = 'X-positive'
   face_set_ids = 1
   type = 'flux'
   flux = 0.0
   /

The MATERIAL namelist defines a material, either single-phase or multi-phase, that is available to be used ina simulation.
::
   
   &MATERIAL
   name = 'graphite'
   density = 1750
   specific_heat = 1500
   conductivity = 100
   /

The DIFFUSION_SOLVER namelist sets the parameters that are specific to the heat and species transport solver.
::

   &DIFFUSION_SOLVER
   abs_temp_tol = 0.0
   rel_temp_tol = 1.0e-3
   abs_enthalpy_tol = 0.0
   rel_enthalpy_tol = 1.0e-3
   max_nlk_itr = 5
   nlk_tol = 0.02
   nlk_preconditioner = 'hypre_amg'
   pc_amg_cycles = 2
   verbose_stepping = .true.
   /

The NUMERICS namelist specifies general numerical parameters not specific to any particular physics,especially those controlling the overall time stepping of the Truchas model.
::
   
   &NUMERICS
   dt_init = 0.01
   dt_grow = 5.0
   dt_max = 60.0
   /


.. toctree::
   :maxdepth: 1
   

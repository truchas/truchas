# Electromagnetics Code Organization

* `em_properties.F90` provides procedures for evaluating the cell-based EM
  properties. It interfaces with the Truchas material database and material
  volume fraction data.

* `em_heat_driver.F90` provides the interface between the Truchas time cycle
  driver and the computation of the time/temperature-dependent electromagnetic
  heat source for heat transfer simulations. It handles the mapping of fields
  between the HT and EM meshes, and manages the input of EM-associated
  namelists, output of restart data, and choice of EM heat solver. Its
  primary dependencies are:

  - `induction_heat_solver_type.F90` for Joule heat driven by magnetic
    induction at low frequencies.
  - `microwave_heat_solver_type.F90` for dielectric heat at microwave
    frequencies.

* `induction_heat_solver_type.F90` manages the computation of Joule heat in
  problems driven by magnetic induction. It determines whether or not the
  Joule heat needs to be updated, based on changes to the induction forcing
  and temperature-dependent EM properties, and if it needs to be updated,
  side-stepping the computation if possible by simple scaling or zeroing.
  Its primary dependencies are:

  - `ih_hfield_factory_type.F90` manages the step-wise, time-dependent
    magnetic induction field source, and the creation of associated
    functions used in BC. Its dependencies are:
    - `ih_hfield_type.F90` computes the magnetic field produced by a set of
      induction coils.
    - `solenoid_fields.F90` computes the magnetic field due to a current loop.
    - `elliptic_integrals.F90` are some special functions
  - `tdme_joule_heat_sim_type.F90` is a periodic time-domain Maxwell equation
    solver for when actual computation is necessary.
  - `fdme_solver_type.F90` is a frequency-domain Maxwell equation solver
    for when actual computation is necessary.

* `microwave_heat_solver_type.F90` manages the computation of dielectric heat
  in microwave problems. It determines whether or not the dielectric heat needs
  to be updated, based on changes in waveguide power input and temperature-
  dependent EM properties, and if it needs to be updated, side-stepping the
  computation if possible by simple scaling or zeroing. Its primary
  dependencies are:

  - `wg_port_bc_plist_factory_type.F90` manages the step-wise, time-dependent
    waveguide input power and the creation of BC parameter list input for the
    waveguide port BC.
  - `fdme_solver_type.F90` is a frequency-domain Maxwell equation solver
    for when actual computation is necessary.

* `tdme_joule_heat_sim_type.F90` computes the Joule heat by solving the time
  dependent Maxwell equations with periodic magnetic induction field forcing
  to a periodic steady state and averaging the periodic Joule heat over one
  period. Its dependencies are:

  - `em_bc_factory_type.F90` creates the boundary condition functions
  - `tdme_model_type.F90` stores the linear system resulting from the
    discretization of Maxwell equations, together with associated methods
    acting on that data (e.g. residual). Its dependencies are:
    - `mimetic_discretization.F90` provides discrete operators and matrices
  - `tdme_solver_type.F90` advances the discrete system by one time step
  - `tdme_cg_solver_type.F90` is a preconditioned CG solver for the linear
    time-step system.

* `fdme_solver.F90` is a general solver for frequency-domain Maxwell equations,
  providing several options for solution method. Its dependencies are:

  - `em_bc_factory_type.F90` creates the boundary condition functions
  - `fdme_model_type.F90` stores the linear system resulting from the
    discretization of the frequency-domain Maxwell equations and its mixed
    formulation, together with associated methods acting on that data.
    Its dependencies are:
    - `mimetic_discretization.F90` provides discrete operators and matrices
    - `index_corrector_type.F90` is used in the conversion of a complex
      matrix to a double-sized real matrix required by the MUMPS solver.
  - `fdme_minres_solver_type.F90` is a low-level solver that uses a
    preconditioned CS-MINRES linear solver for the complex-symmetric system.
  - `fdme_mumps_solver_type.F90` is a low-level solver that uses the MUMPS
    direct linear solver.
  - `fdme_mixed_minres_solver_type.F90` is a low-level solver that uses a
    preconditioned CS-MINRES linear solver for a mixed formulation of the
    frequency-domain Maxwell equations. Its dependencies are:
    - `fdme_mixed_zvector_type.F90` is an abstract vector type that holds
      the double-sized vector of the mixed formulation.

* `fdme_vtk_graphics_proc` provides a procedure that writes a VTKHDF graphics
  file of the solution and derived quantitites computed by the frequency-
  domain Maxwell equation solver.

* `mimetic_discretization.F90` provides discrete operators and inner product
  matrices for the Whitney edge and face finite element spaces on tetrahedral
  meshes.

* `em_bc_factory_type.F90` is a factory that creates BC functions for the
  various types of EM boundary conditions, for both frequency-domain and
  time-domain EM solvers. Its dependencies are:

  - `pec_bndry_func_type.F90` for PEC conditions
  - `nxh_bndry_func_type.F90` for nxH conditions
  - `fd_robin_bndry_func_type.F90` for general frequency-domain Robin conditions
  - `port_feed_func_factory_type.F90` for generating waveguide port feed
    functions for Robin conditions

## Input

Namelist input is orchestrated by `em_heat_driver.F90` which converts it into
a hierarchical parameter list that is then passed down through the call stack
and parceled out to various components. Its dependencies are:

- `electromagnetic_bc_namelist.F90` reads the `electromagnetic_bc` namelists
  and populates a parameter list that is input for `em_bc_factory_type.F90`.
- `induction_source_field_namelist.F90` reads the `induction_source_field`
  namelist and populates a parameter list that is input for
  `ih_hfield_factory_type.F90`
- `electromagnetics_namelist.F90` reads the `electromagnetics` namelist and
  populates a parameter list contains inputs for all other components.

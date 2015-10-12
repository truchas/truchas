MODULE SETUP_MODULE
  !=======================================================================
  ! Purpose:
  !
  !   Define procedures necessary to set up the problem.
  !
  !   Public Interface(s):
  !
  !     * call SETUP ()
  !
  !       Sets up the problem.
  !
  ! Contains: SETUP
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !            Robert Ferrell (CPCA, Ltd., ferrell@cpca.com)
  !=======================================================================
  implicit none

  ! Private module
  private

  ! public procedures
  public :: SETUP

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE SETUP ()
    !=======================================================================
    ! Purpose:
    !
    !    set up the problem
    !=======================================================================
    use kinds, only: r8
    use base_types_A_module,    only: BASE_TYPES_A_ALLOCATE
    use base_types_B_module,    only: MESH_VERTEX_ALLOCATE, BASE_TYPES_B_ALLOCATE
    use bc_module,              only: ASSIGN_BC_BITS, Conc, Prs, Vel
    use cell_geometry_module,   only: GET_CELL_GEOMETRY
    use gs_module,              only: EE_GS_INIT, EN_GS_INIT, NN_GS_INIT
    use init_module,            only: INITIAL
    use mesh_gen_module,        only: MESH_GEN,                    &
                                      Flag_Face_Neighbors
    use mesh_utilities,         only: MESH_DIAGNOSTICS
    use mesh_module,            only: Mesh, Vertex
    use output_control,         only: next_op
    use restart_variables,      only: restart, ignore_t, ignore_dt, restart_t, restart_dt, &
                                      restart_cycle_number
    use restart_driver,         only: close_restart_file
    use time_step_module,       only: cycle_number, cycle_number_restart, &
                                      dt, t, t1, t2
    use timing_tree
    use tensor_module,          only: TENSOR_MATRIX
    use mesh_manager,           only: init_mesh_manager
    use EM,                     only: initialize_EM
    use truchas_danu_output,    only: TDO_write_default_mesh

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the initialization timer.
    call start_timer("Initialization")

    ! Initialize I/O quantities.
    next_op = 1
    
    ! Assign bc bits.
    call ASSIGN_BC_BITS (Prs%Face_bit, Vel%Face_bit, Conc%Face_bit)

    ! Setup the tensor matrix.
    call TENSOR_MATRIX ()

    ! Allocate and default the Mesh and Vertex derived types.
    call MESH_VERTEX_ALLOCATE (Mesh, Vertex)

    ! Read in or generate the mesh; compute the connectivity.
    call MESH_GEN ()

    ! Allocate the base types and set them to their defaults.
    call BASE_TYPES_A_ALLOCATE ()
    call BASE_TYPES_B_ALLOCATE ()

    ! Setup element <-> element (ee) communication.
    call EE_GS_INIT ()

    ! Setup element <-> node (en) communication.
    call EN_GS_INIT ()

    ! Setup node <-> node (en) communication.
    call NN_GS_INIT ()

    ! Compute mesh diagnostics and check for connectivity errors.
    call MESH_DIAGNOSTICS ()

    ! Get the cell geometry.
    call GET_CELL_GEOMETRY ()

    ! Flag which neighbors touch which faces
    call Flag_Face_Neighbors ()

    ! Setup the distributed tet mesh used by the EM solver.
    call init_mesh_manager ()

    ! Write the primary truchas mesh.
    call TDO_write_default_mesh

    ! NNC, Sep 2014.  This used to be done after the call to INITIAL, but is
    ! needed here because of time-dependent boundary conditions.  I left the
    ! rest (DT, CYCLE_NUMBER) where it was because I do not know the
    ! consequences of moving it here, especially CYCLE_NUMBER  as there are
    ! some "if cycle_number == 0" hacks in the flow code.
    if (restart) then
      if (.not.ignore_t)  t  = restart_t
      if (.not.ignore_dt) dt = restart_dt
    end if

    ! Initialize cell-centered fluid variables and thermodynamic quantities.
    call INITIAL (t, dt)

    ! Set the initial timestep value. If this is a restart, take what is in
    ! restart file; if not, then it was already set by the input file.
    if (restart) then
      if (.not.ignore_t) cycle_number = restart_cycle_number
    end if

    ! Initialize electromagnetics.
    call initialize_EM (t)

    if (restart) call close_restart_file ()

    ! Set previous time, last time, and initial cycle number.
    t1 = MAX(t - dt, 0.0_r8)
    t2 = t
    cycle_number_restart = cycle_number

    ! Stop the initialization timer.
    call stop_timer("Initialization")
    
    return

  END SUBROUTINE SETUP

END MODULE SETUP_MODULE

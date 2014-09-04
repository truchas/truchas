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
    use base_types_module,      only: MESH_VERTEX_ALLOCATE, BASE_TYPES_ALLOCATE
    use bc_module,              only: ASSIGN_BC_BITS, Conc, Prs, Vel
    use cell_geometry_module,   only: GET_CELL_GEOMETRY
    use gs_module,              only: EE_GS_INIT, EN_GS_INIT, NN_GS_INIT
    use init_module,            only: INITIAL
    use linear_module,          only: LINEAR_COEF
    use mesh_gen_module,        only: MESH_GEN,                    &
                                      Flag_Face_Neighbors
    use mesh_utilities,         only: MESH_DIAGNOSTICS
    use mesh_module,            only: ASSIGN_CELL_BITS, ASSIGN_VRTX_BITS, CllNgbr, &
                                      Mesh, Vrtx, Vertex,         &
                                      ASSIGN_CELL_EDGES, CELL_EDGE
    use output_control,         only: next_op
    use restart_variables,      only: restart, ignore_dt, restart_t, restart_dt, &
                                      restart_cycle_number
    use restart_driver,         only: close_restart_file
    use time_step_module,       only: cycle_number, cycle_number_restart, &
                                      dt, t, t1, t2
    use timing_tree
    use tensor_module,          only: TENSOR_MATRIX
    use mesh_broker,            only: init_mesh_broker
    use EM,                     only: initialize_EM
#ifdef USE_DANU
    use truchas_danu_output, only: TDO_write_default_mesh
#endif
#ifdef USE_TBROOK
    use output_data_module,     only: enable_tbrook_output
    use tbrook_module,          only: BaseBrook
    use tbrook_utility,         only: TBU_WriteDefaultMesh
    
    integer :: status
#endif

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the initialization timer.
    call start_timer("Initialization")

    ! Initialize I/O quantities.
    next_op = 1
    
    ! Assign bc bits.
    call ASSIGN_BC_BITS (Prs%Face_bit, Vel%Face_bit, Conc%Face_bit)

    ! Assign vertex bits.
    call ASSIGN_VRTX_BITS (Vrtx%Bit)

    ! Assign cell bits.
    call ASSIGN_CELL_BITS (CllNgbr%Bit)

    ! Assign cell edges.
    call ASSIGN_CELL_EDGES (Cell_Edge)

    ! Setup the tensor matrix.
    call TENSOR_MATRIX ()

    ! Setup linear interpolation coefficients.
    call LINEAR_COEF ()

    ! Allocate and default the Mesh and Vertex derived types.
    call MESH_VERTEX_ALLOCATE (Mesh, Vertex)

    ! Read in or generate the mesh; compute the connectivity.
    call MESH_GEN ()

    ! Allocate the base types and set them to their defaults.
    call BASE_TYPES_ALLOCATE ()

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
    call init_mesh_broker ()

#ifdef USE_TBROOK
    if (enable_tbrook_output) then
    ! Write the primary truchas mesh.
    status = 0
    call TBU_WriteDefaultMesh (BaseBrook, status)
    end if
#endif

#ifdef USE_DANU
    ! Write the primary truchas mesh.
    call TDO_write_default_mesh

#endif
    ! NNC, Sep 2014.  This used to be done after the call to INITIAL, but is
    ! needed here because of time-dependent boundary conditions.  I left the
    ! rest (DT, CYCLE_NUMBER) where it was because I do not know the
    ! consequences of moving it here, especially CYCLE_NUMBER  as there are
    ! some "if cycle_number == 0" hacks in the flow code.
    if (restart) t = restart_t

    ! Initialize cell-centered fluid variables and thermodynamic quantities.
    call INITIAL (t)

    ! Set the initial timestep value. If this is a restart, take what is in
    ! restart file; if not, then it was already set by the input file.
    if (restart) then
      if (.not.ignore_dt) dt = restart_dt
      cycle_number = restart_cycle_number
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    use bc_module,              only: ASSIGN_BC_BITS, Conc, Prs, Vel
    use output_control,         only: next_op
    use restart_variables,      only: restart, ignore_t, ignore_dt, restart_t, restart_dt, &
                                      restart_cycle_number
    use restart_driver,         only: close_restart_file
    use time_step_module,       only: cycle_number, cycle_number_restart, &
                                      dt, t, t1, t2
    use init_module,            only: INITIAL
    use timing_tree
    use tensor_module,          only: TENSOR_MATRIX
    use mesh_manager,           only: init_mesh_manager
    use legacy_mesh_api,        only: init_legacy_mesh_api
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

    ! Initialize the legacy mesh data structures.
    call init_legacy_mesh

    ! Instantiate the new mesh objects.
    call init_mesh_manager

    ! Initialize old_mesh_api module
    call init_legacy_mesh_api

    ! Allocate the base types and set them to their defaults.
    call BASE_TYPES_A_ALLOCATE

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

  contains

    subroutine init_legacy_mesh

      use base_types_B_module,    only: MESH_VERTEX_ALLOCATE, BASE_TYPES_B_ALLOCATE
      use cell_geometry_module,   only: GET_CELL_GEOMETRY
      use gs_module,              only: EE_GS_INIT, EN_GS_INIT, NN_GS_INIT
      use mesh_gen_module,        only: MESH_GEN, Flag_Face_Neighbors
      use mesh_utilities,         only: MESH_DIAGNOSTICS
      use mesh_module,            only: Mesh, Vertex

      call MESH_VERTEX_ALLOCATE (Mesh, Vertex) ! allocate the Mesh and Vertex structure arrays
      call MESH_GEN ! read the mesh and compute the connectivity
      call BASE_TYPES_B_ALLOCATE ! allocate the Cell structure array
      call EE_GS_INIT ! setup element-element communication
      call EN_GS_INIT ! setup element-node communication
      call NN_GS_INIT ! setup node-node communication
      call MESH_DIAGNOSTICS ! compute mesh diagnostics and check for connectivity errors
      call GET_CELL_GEOMETRY ! compute cell geometry
      call Flag_Face_Neighbors ! flag which neighbors touch which faces

    end subroutine init_legacy_mesh

  END SUBROUTINE SETUP

END MODULE SETUP_MODULE

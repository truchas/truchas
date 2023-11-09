!TODO: finish documentation
!! HT_2D_MODEL_TYPE
!!
!! David Neill-Asanza <dhna@lanl.gov>
!! July 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module HT_2d_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use mfd_2d_disc_type
  use bndry_func1_class
  use scalar_mesh_func_class
  use TofH_type
  use data_layout_type
  use matl_mesh_func_type
  use prop_mesh_func_type
  use parallel_communication
  use parameter_list_type
  use truchas_logging_services
  implicit none
  private

  type, public :: HT_2d_model
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(mfd_2d_disc),  pointer :: disc => null()
    !! Variable layout
    type(data_layout) :: layout
    integer :: cell_heat_segid, cell_temp_segid, face_temp_segid
    !! Equation parameters
    type(prop_mesh_func) :: conductivity  ! thermal conductivity
    type(prop_mesh_func) :: H_of_T        ! enthalpy as a function of temperature
    type(TofH) :: T_of_H                  ! inverse of enthalpy-temperature function
    class(scalar_mesh_func), allocatable :: source  ! external heat source
    !! Boundary condition data
    class(bndry_func1), allocatable :: bc_dir   ! Dirichlet
    class(bndry_func1), allocatable :: bc_flux  ! Simple flux
  contains
    procedure :: init
    procedure :: num_dof
    procedure :: get_cell_heat_view
    procedure :: get_cell_temp_view
    procedure :: get_face_temp_view
    procedure :: get_cell_heat_copy
    procedure :: get_cell_temp_copy
    procedure :: get_face_temp_copy
    procedure :: set_cell_heat
    procedure :: set_cell_temp
    procedure :: set_face_temp
    procedure :: new_state_array
    procedure :: compute_f
    procedure :: compute_udot
    procedure :: compute_face_temp
  end type HT_2d_model

contains

  subroutine init(this, disc, mmf, params, stat, errmsg)

    use material_model_driver, only: matl_model
    use material_utilities

    class(HT_2d_model), intent(out), target :: this
    type(mfd_2d_disc), intent(in), target :: disc
    type(matl_mesh_func), intent(in), target :: mmf
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: TofH_max_try
    real(r8) :: TofH_tol, TofH_delta
    character(len=256) :: errmsg2 !TODO: remove when possible
    character(:), allocatable :: context
    type(parameter_list), pointer :: sublist

    context = 'processing ' // params%name() // ': '

    this%disc => disc
    this%mesh => disc%mesh

    !! Create the packed layout of the model variables.
    this%cell_heat_segid = alloc_segment(this%layout, this%mesh%ncell_onP)
    this%cell_temp_segid = alloc_segment(this%layout, this%mesh%ncell_onP)
    this%face_temp_segid = alloc_segment(this%layout, this%mesh%nface_onP)
    call alloc_complete(this%layout)

    !! Enthalpy density.
    call required_property_check(matl_model, 'enthalpy', stat, errmsg)
    if (stat /= 0) return
    call this%H_of_T%init(mmf, 'enthalpy', stat, errmsg2)
    if (global_any(stat /= 0)) then
      stat = -1
      errmsg = context // 'unexpected error defining H_of_T: ' // trim(errmsg2)
      return
    end if

    !! Inverse of enthalpy-temperature relation.
    !TODO is "TofH" the right name for this?
    ! if (params%is_sublist('TofH')) then
    !TODO: ok to make entire list optional? ok to create the empty sublist if it does not exist?
    !      or better to not put these in a sublist?
      sublist => params%sublist('TofH')
      !! absolute temperature convergence tolerance, >= 0
      call sublist%get('tolerance', TofH_tol, default=0.0d0, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      !! initial endpoint shift when seeking bracket, > 0
      call sublist%get('delta', TofH_delta, default=1.0d-3, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      !! max tries at seeking a bracketing interval, >= 0
      call sublist%get('max-try', TofH_max_try, default=50, stat=stat, errmsg=errmsg)
      if (stat /= 0) return

      call this%T_of_H%init(this%H_of_T, eps=TofH_tol, max_try=TofH_max_try, delta=TofH_delta)
    ! else
    !   stat = 1
    !   errmsg = context // 'missing "TofH" sublist parameter'
    !   return
    ! end if

    !! Thermal conductivity.
    call required_property_check(matl_model, 'conductivity', stat, errmsg)
    if (stat /= 0) return
    call this%conductivity%init(mmf, 'conductivity', stat, errmsg2)
    if (global_any(stat /= 0)) then
      stat = -1
      errmsg = context // 'unexpected error defining conductivity: ' // trim(errmsg2)
      return
    end if

    !! Defines the boundary condition components
    if (params%is_sublist('bc')) then
      sublist => params%sublist('bc')
      call init_bc(this, sublist, stat, errmsg)
      if (stat /= 0) return
    else
      stat = 1
      errmsg = context // 'missing "bc" sublist parameter'
      return
    end if

    !! Defines the heat source
    if (params%is_sublist('source')) then
      sublist => params%sublist('source')
      call init_source(this, sublist, stat, errmsg)
      if (stat /= 0) return
    else
      !TODO: should it fail if 'source' is specified but not a sublist?
      call TLS_info('No "source" sublist specified')
    end if

  end subroutine init

  subroutine init_bc(model, params, stat, errmsg)

    use bitfield_type
    use thermal_bc_factory1_type
    use string_utilities, only: i_to_c
    use physical_constants, only: stefan_boltzmann, absolute_zero

    class(HT_2d_model), intent(inout), target :: model
    type(parameter_list), intent(inout), target :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(thermal_bc_factory1) :: bc_fac
    type(bitfield) :: bitmask
    character(160) :: string
    logical, allocatable :: mask(:)
    integer, allocatable :: setids(:)
    integer :: j, n

    allocate(mask(model%mesh%nface), source=.false.)

    call bc_fac%init(model%mesh, stefan_boltzmann, absolute_zero, params)

    !! Define the simple flux boundary conditions.
    call bc_fac%alloc_flux_bc(model%bc_flux, stat, errmsg)
    if (stat /= 0) return
    if (allocated(model%bc_flux)) then
      mask(model%bc_flux%index) = .true. ! mark the simple flux faces
    end if

    !! Define the Dirichlet boundary conditions.
    call bc_fac%alloc_dir_bc(model%bc_dir, stat, errmsg)
    if (stat /= 0) return
    if (allocated(model%bc_dir)) then
      if (global_any(mask(model%bc_dir%index))) then
        stat = -1
        errmsg = 'temperature dirichlet boundary condition overlaps with other conditions'
        return
      end if
      mask(model%bc_dir%index) = .true. ! mark the dirichlet faces
    end if

    !! Finally verify that a condition has been applied to every boundary face.
    mask = mask .neqv. btest(model%mesh%face_set_mask,0)
    if (global_any(mask)) then
      call model%mesh%get_face_set_ids(pack([(j,j=1,model%mesh%nface)], mask), setids)
      if (size(setids) == 0) then
        string = '(none)'
      else
        write(string,'(i0,*(:,", ",i0))') setids
      end if
      errmsg = 'incomplete temperature boundary/interface specification;' // &
          ' remaining boundary faces belong to face sets ' // trim(string)
      ! TODO: should this support link sets?
      ! call model%mesh%get_link_set_ids(mask, setids)
      ! if (size(setids) == 0) then
      !   string2 = '(none)'
      ! else
      !   write(string2,'(i0,*(:,", ",i0))') setids
      ! end if
      ! errmsg = 'incomplete temperature boundary/interface specification;' // &
      !     ' remaining boundary faces belong to face sets ' // trim(string) // &
      !     '; and interface link sets ' // trim(string2)
      bitmask = ibset(ZERO_BITFIELD, 0)
      mask = mask .and. (model%mesh%face_set_mask == bitmask)
      ! mask(model%mesh%lface(1,:)) = .false.
      ! mask(model%mesh%lface(2,:)) = .false.
      n = global_count(mask(:model%mesh%nface_onP))
      if (n > 0) errmsg = errmsg // '; ' // i_to_c(n) // ' faces belong to neither'
      stat = -1
      return
    end if

  end subroutine init_bc

  subroutine init_source(model, params, stat, errmsg)

    use thermal_source_factory_type

    class(HT_2d_model), intent(inout), target :: model
    type(parameter_list), intent(inout), target :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(thermal_source_factory) :: src_fac

    call src_fac%init(model%mesh, params)

    !! Allocated function-based source
    call src_fac%alloc_source_func2(model%source, stat, errmsg)
    if (stat /= 0) return

    !TODO: check all cells have a source, if a source was specified?

  end subroutine init_source


  subroutine get_cell_heat_view(this, array, view)
    class(HT_2d_model), intent(in) :: this
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    call get_segment_view(this%layout, array, this%cell_heat_segid, view)
  end subroutine get_cell_heat_view

  subroutine get_cell_temp_view(this, array, view)
    class(HT_2d_model), intent(in) :: this
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    call get_segment_view(this%layout, array, this%cell_temp_segid, view)
  end subroutine get_cell_temp_view

  subroutine get_face_temp_view(this, array, view)
    class(HT_2d_model), intent(in) :: this
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    call get_segment_view(this%layout, array, this%face_temp_segid, view)
  end subroutine get_face_temp_view

  subroutine get_cell_heat_copy(this, array, copy)
    class(HT_2d_model), intent(in) :: this
    real(r8), intent(in), target :: array(:)
    real(r8), intent(inout) :: copy(:)
    call get_segment_copy(this%layout, array, this%cell_heat_segid, copy)
  end subroutine get_cell_heat_copy

  subroutine get_cell_temp_copy(this, array, copy)
    class(HT_2d_model), intent(in) :: this
    real(r8), intent(in), target :: array(:)
    real(r8), intent(inout) :: copy(:)
    call get_segment_copy(this%layout, array, this%cell_temp_segid, copy)
  end subroutine get_cell_temp_copy

  subroutine get_face_temp_copy(this, array, copy)
    class(HT_2d_model), intent(in) :: this
    real(r8), intent(in), target :: array(:)
    real(r8), intent(inout) :: copy(:)
    call get_segment_copy(this%layout, array, this%face_temp_segid, copy)
  end subroutine get_face_temp_copy

  subroutine set_cell_heat(this, source, array)
    class(HT_2d_model), intent(in) :: this
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout) :: array(:)
    call set_segment(this%layout, source, array, this%cell_heat_segid)
  end subroutine set_cell_heat

  subroutine set_cell_temp(this, source, array)
    class(HT_2d_model), intent(in) :: this
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout) :: array(:)
    call set_segment(this%layout, source, array, this%cell_temp_segid)
  end subroutine set_cell_temp

  subroutine set_face_temp(this, source, array)
    class(HT_2d_model), intent(in) :: this
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout) :: array(:)
    call set_segment(this%layout, source, array, this%face_temp_segid)
  end subroutine set_face_temp


  integer function num_dof(this)
    class(HT_2d_model), intent(in) :: this
    num_dof = layout_size(this%layout)
  end function num_dof

  subroutine new_state_array(this, u, state)
    class(HT_2d_model) :: this
    real(r8), intent(in) :: u(:)
    real(r8), allocatable, intent(out) :: state(:,:)
    allocate(state(this%mesh%ncell,0:0))
    call this%get_cell_temp_copy(u, state(:,0))
    call this%mesh%cell_imap%gather_offp(state(:,0))
  end subroutine new_state_array


  subroutine compute_f(this, t, u, udot, f)

    class(HT_2d_model), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), target :: u(:), udot(:)
    real(r8), intent(out), target :: f(:)

    real(r8), dimension(this%mesh%ncell) :: Tcell, Fcell, Hdot
    real(r8), dimension(this%mesh%nface) :: Tface, Fface
    real(r8), pointer :: Hcell(:), FHcell(:)
    real(r8), allocatable :: Fdir(:), state(:,:)
    real(r8) :: cval(this%mesh%ncell)

    call this%new_state_array(u, state)

    !!!! RESIDUAL OF THE ALGEBRAIC ENTHALPY-TEMPERATURE RELATION !!!!!!!!!!!!!!!!!

      call this%get_cell_heat_view(u, Hcell)
      ! call this%get_cell_temp_view(u, Tcell)
      call this%get_cell_heat_view(f, FHcell)
      call this%H_of_T%compute_value(state, cval)
      FHcell = Hcell - cval(:this%mesh%ncell_onP)

    !!!! RESIDUAL OF THE HEAT EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !! Off-process cell and face temperatures
      call this%get_cell_temp_copy(u, Tcell)
      call this%get_face_temp_copy(u, Tface)
      call this%mesh%cell_imap%gather_offp(Tcell)
      call this%mesh%face_imap%gather_offp(Tface)

      !! Off-process cell enthalpy time derivative
      call this%get_cell_heat_copy(udot, Hdot)
      call this%mesh%cell_imap%gather_offp(Hdot)

      !! Pre-compute the Dirichlet condition residual and
      !! impose the Dirichlet data on the face temperature.
      if (allocated(this%bc_dir)) then
        call this%bc_dir%compute(t)
        allocate(Fdir(size(this%bc_dir%index)))
        associate (index => this%bc_dir%index, value => this%bc_dir%value)
          Fdir = Tface(index) - value
          Tface(index) = value
        end associate
      end if

      !! Compute the generic heat equation residual.
      call this%conductivity%compute_value(state, cval)
      call this%disc%apply_diff(cval, Tcell, Tface, Fcell, Fface)
      Fcell = Fcell + this%mesh%volume*Hdot

      !! Optional source function contribution
      if (allocated(this%source)) then
        call this%source%compute(t)
        Fcell = Fcell - this%mesh%volume*this%source%value
      end if

      !! Dirichlet condition residuals.
      if (allocated(this%bc_dir)) then
        associate (index => this%bc_dir%index)
          Fface(index) = Fdir  ! overwrite with pre-computed values
        end associate
        deallocate(Fdir)
      end if

      !! Simple flux BC contribution.
      if (allocated(this%bc_flux)) then
        call this%bc_flux%compute(t)
        associate (index => this%bc_flux%index, value => this%bc_flux%value)
          Fface(index) = Fface(index) + this%mesh%area(index) * value
        end associate
      end if

      !! Return the on-process part of the heat conduction residuals
      call this%set_cell_temp(Fcell, f)
      call this%set_face_temp(Fface, f)

  end subroutine compute_f

  !TODO: is there a better way to set/get rel_tol, max_itr
  !! Computes UDOT given a consistent U

  subroutine compute_udot(this, t, dt, u, udot, rel_tol, max_itr, stat, errmsg)

    class(HT_2d_model), intent(inout) :: this
    real(r8), intent(in) :: t, dt, rel_tol
    real(r8), intent(in),  target :: u(:)
    real(r8), intent(out), target :: udot(:)
    integer, intent(in) :: max_itr
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8), allocatable, target :: f(:)
    real(r8), pointer :: H(:), FH(:), Hdot(:)
    real(r8), pointer :: Tcell(:), Tface(:), Fcell(:), Fface(:)
    real(r8) :: Tmax, Tmin
    integer :: j

    ASSERT(size(u) == this%num_dof())
    ASSERT(size(udot) == size(u))

    call TLS_info ('')
    call TLS_info ('Computing consistent initial state derivative for HT solver ...')

    !! The DAE system F(t,u,udot) = 0 gives the time derivative of the cell
    !! enthalpy as a a function of the cell and face temperatures.  We back
    !! out what it is by computing F with udot set equal 0.  The info is
    !! contained in the cell temperature section of F.  By construction, the
    !! the remaining sections should be zero (the cell enthalpy section to
    !! round-off, and the face temperature section to the solver tolerance).

    allocate(f(size(u)))
    call this%get_cell_heat_view(f, FH)     ! enthalpy / enthalpy-temp AE
    call this%get_cell_temp_view(f, Fcell)  ! cell temp / heat conduction DE
    call this%get_face_temp_view(f, Fface)  ! face temp / face-cell temp AE

    udot = 0.0_r8
    call this%compute_f(t, u, udot, f)

    call this%get_cell_heat_view(udot, Hdot)
    Hdot = -Fcell / this%mesh%volume(:size(Hdot))

    !! The time derivative of the cell and face temperatures are approximated
    !! by a finite difference.  The enthalpy is advanced by a small time step
    !! using its time derivative (forward Euler), and then associated advanced
    !! cell and face temperatures are solved for using the algebraic relations.

    call this%get_cell_heat_view(u, H)
    call this%get_cell_temp_view(u, Tcell)
    call this%get_face_temp_view(u, Tface)

    FH = H + dt*Hdot ! advance the enthalpy

    !! Compute advanced cell temperatures
    do j = 1, size(Tcell)
      !TODO: special value for void cells?
      Tmin = Tcell(j) - 1
      Tmax = Tcell(j) + 1
      call this%T_of_H%compute(j, H(j), Tmin, Tmax, Fcell(j))
    end do

    Fface = Tface  ! initial guess (probably not half bad)
    call this%compute_face_temp(t+dt, f, rel_tol, max_itr, stat, errmsg)
    if (stat /= 0) return

    Fcell = (Fcell - Tcell) / dt
    call this%set_cell_temp(Fcell, udot)

    Fface = (Fface - Tface) / dt
    call this%set_face_temp(Fface, udot)

  end subroutine compute_udot

  !TODO: is there a better way to set/get rel_tol, max_itr? 4 options:
  !        - Choose whether they belong to model or solver.
  !          1. part of init. Store as components directly.
  !          2. Store reference to parameter list as component. Get them when needed.
  !! Solves for the algebraic equation of the heat transfer model for the face temperatures

  subroutine compute_face_temp(this, t, u, rel_tol, max_itr, stat, errmsg)

    use hypre_hybrid_type
    use mfd_2d_diff_matrix_type

    class(HT_2d_model), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(inout), target :: u(:)
    integer, intent(in) :: max_itr
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(mfd_2d_diff_matrix), target :: dm
    type(parameter_list), target :: params
    type(hypre_hybrid) :: solver
    real(r8), allocatable :: coef(:), udot(:), z(:), state(:,:)
    real(r8), allocatable, target :: r(:)
    real(r8), pointer :: Tface(:), rface(:)
    real(r8) :: norm, rel_tol
    integer :: n, num_itr, num_dscg_itr, num_pcg_itr
    character(80) :: msg

    ASSERT(size(u) == this%num_dof())

    !! Compute the RHS
    n = this%num_dof()
    allocate(udot(n), r(n))
    udot = 0.0_r8
    call this%compute_f(t, u, udot, r)
    call this%get_face_temp_view(r, rface)
    norm = sqrt(global_dot_product(rface, rface))

    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      write (msg,'(a,es10.3)') 'HT_2D_model%compute_face_temp: initial ||rface||_2 = ', norm
      call TLS_info (trim(msg))
    end if
    if (norm == 0.0_r8) return

    !! Compute the matrix (the A22 submatrix)
    allocate(coef(this%mesh%ncell))
    call this%new_state_array(u, state)
    call this%conductivity%compute_value(state, coef)
    call dm%init(this%disc)
    call dm%compute(coef)
    if (allocated(this%bc_dir)) then
      call this%bc_dir%compute(t)
      call dm%set_dir_faces(this%bc_dir%index)
    end if
    deallocate(coef, state)

    !! Setup the linear solver.
    call params%set('krylov-method', 'cg')
    call params%set('max-ds-iter', max_itr)
    call params%set('max-amg-iter', max_itr)
    call params%set('rel-tol', rel_tol)
    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      call params%set('print-level', 1)
      call params%set('logging-level', 1)
    else
      call params%set('print-level', 0)
      call params%set('logging-level', 0)
    end if
    call solver%init(dm%a22, params)
    call solver%setup()

    !! Solve
    call this%get_face_temp_view(u, Tface)
    allocate(z(size(Tface)))
    z = 0.0_r8
    call solver%solve(rface, z, stat)
    Tface = Tface - z

    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      call solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, norm)
      write(msg,'(3(a,i0),a,es9.2)') 'solve: num_itr = ', num_itr, &
          ' (', num_dscg_itr, ', ', num_pcg_itr, '), ||r||/||b|| = ', norm
      call TLS_info(trim(msg))

      !! Check the residual.
      call this%compute_f(t, u, udot, r)
      write(msg,'(a,es9.2)') 'solve: ||rface||_2 = ', sqrt(global_dot_product(rface,rface))
      call TLS_info(trim(msg))
    end if

    if (stat /= 0) then
      call solver%get_metrics(num_itr, num_dscg_itr, num_pcg_itr, norm)
      write(msg,'(3(a,i0),a,es9.2)') 'failed to converge: num_itr = ', num_itr, &
          ' (', num_dscg_itr, ', ', num_pcg_itr, '), ||r||/||b|| = ', norm
      errmsg = trim(msg)
    end if

  end subroutine compute_face_temp

end module HT_2d_model_type

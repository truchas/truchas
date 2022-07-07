!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module flow_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use truchas_timers
  use parameter_list_type
  use bndry_func1_class
  use bndry_vfunc_class
  use flow_surface_tension_bc_type
  use unstr_mesh_type
  use flow_bc_factory_type
  implicit none
  private

  type, public :: flow_bc
    class(bndry_func1), allocatable :: p_dirichlet, dp_dirichlet, p_neumann, v_zero_normal
    class(bndry_vfunc), allocatable :: v_dirichlet
    type(surface_tension_bc) :: surface_tension
    logical :: pressure_d
    type(parameter_list), pointer :: inflow_plist ! OWNED
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_initial
    procedure :: is_p_neumann_fix_pe
    procedure, private :: apply_default
  end type flow_bc

contains

  subroutine init(this, mesh, params, vof, temperature_fc, stat, errmsg)

    use parameter_list_type
    use flow_props_type
    use parallel_communication, only: global_sum

    class(flow_bc), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params
    real(r8), intent(in), target :: vof(:), temperature_fc(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(flow_bc_factory) :: f
    type(parameter_list), pointer :: plist

    plist => params%sublist('bc')
    call f%init(mesh, plist)
    call f%alloc_dir_vel_bc(this%v_dirichlet, stat, errmsg)
    if (stat /= 0) return
    call f%alloc_dir_prs_bc(this%p_dirichlet, stat, errmsg)
    if (stat /= 0) return
    call f%alloc_dir_prs_bc(this%dp_dirichlet, stat, errmsg)
    if (stat /= 0) return
    call f%alloc_neu_prs_bc(this%p_neumann, stat, errmsg)
    if (stat /= 0) return
    call f%alloc_zero_vn_bc(this%v_zero_normal, stat, errmsg)
    if (stat /= 0) return

    !TODO: incorporate this into the bc factory
    call this%surface_tension%init(plist, mesh, vof, temperature_fc)

    ! apply the default boundary condition on faces with no user-specified BC
    call this%apply_default(mesh)

    this%pressure_d = global_sum(size(this%p_dirichlet%index)) > 0

    allocate(this%inflow_plist)
    call make_inflow_plist(plist, this%inflow_plist)

  end subroutine init

  ! Assigns an MPI rank to do the pressure-neumann null space fixup.
  ! The rank must have fluid in the current time step, and otherwise
  ! neumann pressure boundaries.

  logical function is_p_neumann_fix_pe(this, any_real_fluid_onP)

    use parallel_communication

    class(flow_bc), intent(inout) :: this
    logical, intent(in) :: any_real_fluid_onP

    integer :: fix_pe
    logical :: is_valid_pe, has_p_neumann(nPE)

    is_p_neumann_fix_PE = .false.

    if (.not.this%pressure_d) then
      ! find the lowest PE with a pressure neumann boundary
      is_valid_pe = size(this%p_neumann%index) > 0 .and. any_real_fluid_onP
      call gather(is_valid_pe, has_p_neumann)
      if (is_IOP) fix_pe = findloc(has_p_neumann, .true., dim=1)
      call broadcast(fix_pe)
      INSIST(fix_pe > 0 .and. fix_pe <= nPE)
      is_p_neumann_fix_PE = (fix_pe == this_PE)
    end if

  end function is_p_neumann_fix_pe

  ! apply the default boundary condition (slip) to
  ! boundary faces with unspecified conditions

  subroutine apply_default(this, mesh)

    use bndry_face_func_type
    use scalar_func_class
    use scalar_func_factories, only: alloc_const_scalar_func

    class(flow_bc), intent(inout) :: this
    type(unstr_mesh), intent(in) :: mesh

    integer :: f, nf
    integer, allocatable :: faces(:)
    class(scalar_func), allocatable :: func

    ! get list of boundary faces with no user-applied BC
    allocate(faces(mesh%nface))

    nf = 0
    do f = 1, mesh%nface_onP
      if (mesh%fcell(2,f) /= 0) cycle ! skip internal faces

      ! skip faces with user-given BC
      if (any(this%p_dirichlet%index == f)) cycle
      if (any(this%p_neumann%index == f)) cycle ! catches v_dirichlet & surface_tension

      nf = nf + 1
      faces(nf) = f
    end do

    if (nf == 0) return

    ! apply slip bc (p_neumann + v_zero_normal) to faces in faces(:nf)
    call alloc_const_scalar_func(func, 0.0_r8)

    select type (bc => this%p_neumann)
    type is (bndry_face_func)
      call bc%add_face_list(func, faces(:nf))
    class default
      ASSERT(.false.)
    end select

    select type (bc => this%v_zero_normal)
    type is (bndry_face_func)
      call bc%add_face_list(func, faces(:nf))
    class default
      ASSERT(.false.)
    end select

  end subroutine apply_default

  subroutine compute(this, t, dt)
    class(flow_bc), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    call this%p_dirichlet%compute(t)
    call this%p_neumann%compute(t)
    call this%v_dirichlet%compute(t)
    call this%surface_tension%compute(t)
    call this%dp_dirichlet%compute(t+dt)
    this%dp_dirichlet%value = this%dp_dirichlet%value - this%p_dirichlet%value
  end subroutine compute

  !! Compute BC data for calculating the initial flow state. Only the Dirichlet
  !! data for dp is different. Normally the data is the difference in pressure
  !! data over the time step, but in this case we are not taking an actual time
  !! step, so dp is 0 on pressure Dirichlet boundaries.

  subroutine compute_initial(this, t)
    class(flow_bc), intent(inout) :: this
    real(r8), intent(in) :: t
    call this%p_dirichlet%compute(t)
    call this%p_neumann%compute(t)
    call this%v_dirichlet%compute(t)
    call this%surface_tension%compute(t)
    this%dp_dirichlet%value = 0.0_r8
  end subroutine compute_initial

  !! This auxiliary subroutine filters the BC parameter list, extracting extra
  !! data associated with material inflow boundaries into a new parameter list.
  !! The input parameter list is a list of BC sublists. Those prescribing
  !! conditions where material inflow may occur are examined for the presence
  !! of parameters with names beginning with 'inflow'. If any are found, a
  !! sublist in the output parameter list is created with the same name, and
  !! the 'face sets' and all 'inflow*' parameters copied to the sublist. Only
  !! scalar integer, 8-byte real, and character valued inflow parameters are
  !! handled.

  subroutine make_inflow_plist(bc_params, inflow_params)

    use string_utilities, only: lower_case

    type(parameter_list), intent(inout) :: bc_params, inflow_params

    type(parameter_list_iterator) :: piter, ibc_piter
    type(parameter_list), pointer :: ibc, obc
    character(:), allocatable :: bc_type, pname
    logical :: found_inflow
    integer, allocatable :: setids(:)
#ifdef INTEL_BUG20180115
    class(*), pointer :: pval
#endif

    piter = parameter_list_iterator(bc_params, sublists_only=.true.)
    do while (.not.piter%at_end())
      ibc => piter%sublist()
      call ibc%get('type', bc_type)
      select case (lower_case(bc_type))
      case ('velocity', 'pressure')
        found_inflow = .false.
        ibc_piter = parameter_list_iterator(ibc)
        do while (.not.ibc_piter%at_end())
          pname = ibc_piter%name()
          if (pname(:min(6,len(pname))) == 'inflow') then
            if (.not.found_inflow) then
              found_inflow = .true.
              obc => inflow_params%sublist(piter%name())
              call ibc%get('face-set-ids', setids)
              call obc%set('face-set-ids', setids)
            end if
            if (ibc_piter%is_scalar()) then
#ifdef INTEL_BUG20180115
              pval => ibc_piter%scalar()
              select type (pval)
#else
              select type (pval => ibc_piter%scalar())
#endif
              type is (integer)
                call obc%set(pname, pval)
              type is (real(r8))
                call obc%set(pname, pval)
              type is (character(*))
                call obc%set(pname, pval)
              class default
                INSIST(.false.)
              end select
            else
              INSIST(.false.)
            end if
          end if
          call ibc_piter%next
        end do
      end select
      call piter%next()
    end do

  end subroutine make_inflow_plist

end module flow_bc_type

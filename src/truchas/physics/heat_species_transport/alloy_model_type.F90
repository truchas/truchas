!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module alloy_model_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use alloy_vector_type
  use alloy_lever_rule_type
  use alloy_back_diff_type
  use mfd_disc_type
  use scalar_func_class
  use bndry_func1_class
  use bndry_func2_class
  use scalar_mesh_multifunc_type
  use truchas_timers
  use parameter_list_type
  implicit none
  private

  type, public :: alloy_model
    type(unstr_mesh), pointer :: mesh => null() ! unowned reference
    type(mfd_disc) :: disc
    !! Equation parameters
    type(alloy_lever_rule) :: lever
    type(alloy_back_diff) :: back_diff
    class(scalar_func), allocatable :: k_sol, k_liq ! thermal conductivity
    real(r8), allocatable :: q_adv(:) ! advective source
    type(scalar_mesh_multifunc), allocatable :: src ! external heat source
    !! Boundary condition data
    class(bndry_func1), allocatable :: bc_dir  ! Dirichlet
    class(bndry_func1), allocatable :: bc_flux ! simple flux
    class(bndry_func2), allocatable :: bc_htc  ! external HTC
    integer :: model_type, num_comp
  contains
    procedure :: init
    procedure :: init_vector
    procedure :: compute_f
    procedure :: get_conductivity
    procedure :: set_heat_source
    procedure :: get_liq_conc, get_sol_conc
  end type alloy_model

contains

  subroutine init_vector(this, vec)
    class(alloy_model), intent(in) :: this
    type(alloy_vector), intent(out) :: vec
    select case (this%model_type)
    case (1) ! lever rule
      call vec%init(this%mesh)
    case (2) ! Wang-Beckermann
      call vec%init(this%mesh, num_comp=this%back_diff%num_comp)
    end select
  end subroutine

  subroutine init(this, mesh, matl, params, stat, errmsg)

    use thermal_bc_factory1_type
    use material_class

    class(alloy_model), intent(out), target :: this
    type(unstr_mesh), intent(in), target :: mesh
    class(material), intent(in) :: matl
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable :: errmsg

    integer :: j
    type(phase), pointer :: phi
    real(r8), allocatable :: C0(:)

    this%mesh => mesh
    call this%disc%init(this%mesh, use_new_mfd=.true.)

    this%model_type = merge(2, 1, params%is_parameter('gamma'))

    select case (this%model_type)
    case (1) ! lever rule
      call this%lever%init(matl, params, stat, errmsg)
      if (stat /= 0) return
      this%num_comp = this%lever%num_comp
    case (2) ! Wang-Beckermann
      call this%back_diff%init(matl, params, stat, errmsg)
      if (stat /= 0) return
      this%num_comp = this%back_diff%num_comp
    end select

    if (matl%has_prop('conductivity')) then
      phi => matl%phase_ref(1)
      call phi%get_prop('conductivity', this%k_sol)
      phi => matl%phase_ref(2)
      call phi%get_prop('conductivity', this%k_liq)
    else
      stat = 1
      errmsg = 'conductivity property is undefined for material "' // matl%name // '"'
      return
    end if

    !! External heat source.
    block
      use thermal_source_factory_type
      type(thermal_source_factory) :: src_fac
      type(parameter_list), pointer :: plist
      plist => params%sublist('sources')
      call src_fac%init(this%mesh, plist)
      call src_fac%alloc_source_funcs(this%src, stat, errmsg)
      if (stat /= 0) return
    end block

    !! Defines the boundary condition components of MODEL.
    call define_system_bc(stat, errmsg)
    if (stat /= 0) return

  contains

    subroutine define_system_bc(stat, errmsg)

      use bitfield_type
      use parallel_communication, only: global_any, global_count
      use string_utilities, only: i_to_c

      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      integer :: j, n
      logical, allocatable :: mask(:)
      integer, allocatable :: setids(:)
      character(160) :: string
      type(bitfield) :: bitmask
      type(parameter_list), pointer :: plist
      type(thermal_bc_factory1) :: bc_fac

      plist => params%sublist('bc')
      call bc_fac%init(this%mesh, 0.0_r8, 0.0_r8, plist) ! dummy values for stefan_boltzmann and absolute_zero

      allocate(mask(this%mesh%nface))
      mask = .false.

      !! Simple flux boundary conditions
      call bc_fac%alloc_flux_bc(this%bc_flux, stat, errmsg)
      if (stat /= 0) return
      if (allocated(this%bc_flux)) mask(this%bc_flux%index) = .true.

      !! External HTC boundary conditions
      call bc_fac%alloc_htc_bc(this%bc_htc, stat, errmsg)
      if (stat /= 0) return
      if (allocated(this%bc_htc)) mask(this%bc_htc%index) = .true.

      !! Dirichlet boundary conditions
      call bc_fac%alloc_dir_bc(this%bc_dir, stat, errmsg)
      if (stat /= 0) return
      if (allocated(this%bc_dir)) then
        if (global_any(mask(this%bc_dir%index))) then
          stat = -1
          errmsg = 'temperature dirichlet boundary condition overlaps with other conditions'
          return
        end if
        mask(this%bc_dir%index) = .true. ! mark the dirichlet faces
      end if

      !! Finally verify that a condition has been applied to every boundary face.
      mask = mask .neqv. btest(this%mesh%face_set_mask,0)
      if (global_any(mask)) then
        call this%mesh%get_face_set_ids(pack([(j,j=1,this%mesh%nface)], mask), setids)
        if (size(setids) == 0) then
          string = '(none)'
        else
          write(string,'(i0,*(:,", ",i0))') setids
        end if
        errmsg = 'incomplete temperature boundary specification;' // &
            ' remaining boundary faces belong to face sets ' // trim(string)
        bitmask = ibset(ZERO_BITFIELD, 0)
        mask = mask .and. (this%mesh%face_set_mask == bitmask)
        n = global_count(mask(:this%mesh%nface_onP))
        if (n > 0) errmsg = errmsg // '; ' // i_to_c(n) // ' faces belong to none'
        stat = -1
        return
      end if

    end subroutine define_system_bc

  end subroutine

  subroutine set_heat_source(this, q)
    class(alloy_model), intent(inout) :: this
    real(r8), intent(in) :: q(:)
    ASSERT(size(q) == this%mesh%ncell)
    if (.not.allocated(this%q_adv)) allocate(this%q_adv(this%mesh%ncell))
    this%q_adv = q
  end subroutine

  subroutine compute_f(this, C, Cdot, t, u, udot, f)

    use mfd_disc_type

    class(alloy_model), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in) :: C(:,:), Cdot(:,:)
    type(alloy_vector), intent(inout) :: u, udot ! data is intent(in)
    type(alloy_vector), intent(inout) :: f       ! data is intent(out)
    target :: u

    integer :: j, n
    real(r8), dimension(this%mesh%ncell) :: value
    real(r8), allocatable :: Tdir(:)

    call start_timer('ht-function')

    call u%gather_offp
    call udot%gather_offp

    select case (this%model_type)
    case (1) ! lever rule
      call this%lever%compute_f(C, u%lf, u%hc, u%tc, f%lf, f%hc)
    case (2) ! Wang-Beckermann
      do j = 1, this%mesh%ncell
        call this%back_diff%compute_f(C(:,j), Cdot(:,j), &
                               u%lsf(:,j), u%lf(j), u%hc(j), u%tc(j), &
                               udot%lsf(:,j), udot%lf(j), udot%hc(j), udot%tc(j), &
                               f%lsf(:,j), f%lf(j), f%hc(j))
      end do
    end select

  !!!! RESIDUAL OF THE HEAT EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Overwrite the temperature on Dirichlet faces with the boundary
    !! data, saving the original values to restore them later.
    if (allocated(this%bc_dir)) then
      call this%bc_dir%compute(t)
      allocate(Tdir(size(this%bc_dir%index)))
      do j = 1, size(this%bc_dir%index)
        n = this%bc_dir%index(j)
        Tdir(j) = u%tf(n)
        u%tf(n) = this%bc_dir%value(j)
      end do
    end if

    !! Compute the generic heat equation residual.
    call this%get_conductivity(u%lf, u%tc, value)
    call this%disc%apply_diff(value, u%tc, u%tf, f%tc, f%tf)
    call this%mesh%cell_imap%gather_offp(udot%hc)
    if (allocated(this%q_adv)) then
      f%tc = f%tc + this%mesh%volume*(udot%hc - this%q_adv)
    else
      f%tc = f%tc + this%mesh%volume*udot%hc
    end if

    !! Additional heat source
    if (allocated(this%src)) then
      call this%src%compute(t, u%tc)
      f%tc = f%tc - this%mesh%volume*this%src%value
    end if

    !! Overwrite face residuals with the Dirichlet BC residual and
    !! restore the face temperatures to their original input values.
    if (allocated(this%bc_dir)) then
      do j = 1, size(this%bc_dir%index)
        n = this%bc_dir%index(j)
        f%tf(n) = Tdir(j) - this%bc_dir%value(j)
        u%tf(n) = Tdir(j)
      end do
      deallocate(Tdir)
    end if

    !! Simple flux BC contribution.
    if (allocated(this%bc_flux)) then
      call this%bc_flux%compute(t)
      do j = 1, size(this%bc_flux%index)
        n = this%bc_flux%index(j)
        f%tf(n) = f%tf(n) + this%mesh%area(n) * this%bc_flux%value(j)
      end do
    end if

    !! External HTC flux contribution.
    if (allocated(this%bc_htc)) then
      call this%bc_htc%compute(t, u%tf)
      do j = 1, size(this%bc_htc%index)
        n = this%bc_htc%index(j)
        f%tf(n) = f%tf(n) + this%bc_htc%value(j)
      end do
    end if

    !TODO: is this necessary? Off-process values are not needed, but may be
    !      used in dummy vector operations, and we don't want fp exceptions.
    call f%gather_offp

    call stop_timer('ht-function')

  end subroutine compute_f

  !! This auxiliary subroutine computes the phase-averaged thermal conductivity
  subroutine get_conductivity(this, lfrac, temp, k_avg)
    class(alloy_model), intent(in) :: this
    real(r8), intent(in)  :: lfrac(:), temp(:)
    real(r8), intent(out) :: k_avg(:)
    integer :: j
    do j = 1, size(k_avg)
      k_avg(j) = (1-lfrac(j))*this%k_sol%eval([temp(j)]) + lfrac(j)*this%k_liq%eval([temp(j)])
    end do
  end subroutine

  subroutine get_liq_conc(this, C, u, n, C_liq)
    class(alloy_model), intent(inout) :: this
    real(r8), intent(in) :: C(:,:)
    type(alloy_vector), intent(in) :: u
    integer, intent(in) :: n
    real(r8), intent(out) :: C_liq(:)
    integer :: j
    ASSERT(size(C_liq) >= this%mesh%ncell_onp)
    select case (this%model_type)
    case (1) ! lever rule
      call this%lever%compute_C_liq(C, u%hc, n, C_liq)
    case (2) ! Wang-Beckermann
      do j = 1, this%mesh%ncell_onp
        if (u%lf(j) > 1d-6) then
          C_liq(j) = u%lsf(n,j) / u%lf(j)
        else
          C_liq(j) = 0.0_r8 ! really undefined
        end if
      end do
    end select
  end subroutine

  subroutine get_sol_conc(this, C, u, n, C_sol)
    class(alloy_model), intent(inout) :: this
    real(r8), intent(in) :: C(:,:)
    type(alloy_vector), intent(in) :: u
    integer, intent(in) :: n
    real(r8), intent(out) :: C_sol(:)
    integer :: j
    ASSERT(size(C_sol) >= this%mesh%ncell_onp)
    select case (this%model_type)
    case (1) ! lever rule
      call this%lever%compute_C_sol(C, u%hc, n, C_sol)
    case (2) ! Wang-Beckermann
      do j = 1, this%mesh%ncell_onp
        if (1 - u%lf(j) > 1d-6) then
          C_sol(j) = (C(n,j) - u%lsf(n,j)) / (1 - u%lf(j))
        else
          C_sol(j) = 0.0_r8 ! really undefined
        end if
      end do
    end select
  end subroutine

end module alloy_model_type

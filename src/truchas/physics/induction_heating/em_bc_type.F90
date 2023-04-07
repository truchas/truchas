!!
!! EM_BC_type
!!
!! This module provides a type encapsulating EM BC data.
!!
!! Zach Jibben <zjibben@lanl.gov>
!! April 2023
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module em_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use,intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use truchas_logging_services
  use simpl_mesh_type
  use scalar_func_class
  use vector_func_class
  use EM_boundary_data, only: efield_bc, hfield_bc
  implicit none
  private

  type, public :: em_bc
    private
    logical, allocatable, public :: is_ebc_edge(:)
    real(r8), allocatable, public :: efield(:)
    real(r8), allocatable, public :: hsource(:)

    logical :: use_custom_bcs
    integer, allocatable :: hbface(:)
    real(r8), allocatable :: hfield_polarization(:)
    type(simpl_mesh), pointer :: mesh ! unowned reference
    class(scalar_func), allocatable :: hfield_func
    class(vector_func), allocatable :: hfield_vfunc
  contains
    procedure :: init
    procedure :: compute
    procedure, private :: init_custom
    procedure, private :: compute_hb_custom
    procedure, private :: face_edge_projection
    procedure, private :: hfield
  end type em_bc

contains

  subroutine init(this, mesh, use_custom_bcs, params)

    use parameter_list_type

    class(em_bc), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    logical, intent(in) :: use_custom_bcs
    type(parameter_list), intent(inout), pointer :: params

    this%mesh => mesh
    this%use_custom_bcs = use_custom_bcs
    allocate(this%efield(mesh%nedge), this%hsource(mesh%nedge))

    if (this%use_custom_bcs) then
      ! using custom BCs
      call this%init_custom(params)
    else
      ! using automatic cylinder BCs
      ASSERT(allocated(efield_bc%ebedge))
      this%is_ebc_edge = efield_bc%ebedge /= 0
      this%efield = 0
      this%hsource = 0
    end if

  end subroutine init


  !! For now, this is hardwired to only accept one E BC and/or one B BC.
  !! I.e., at least one must be provided, and multiple of either one is
  !! not allowed.
  subroutine init_custom(this, params)

    use parameter_list_type
    use bitfield_type
    use scalar_func_factories, only: alloc_scalar_func
    use vector_func_factories, only: alloc_vector_func

    class(em_bc), intent(inout) :: this
    type(parameter_list), intent(inout), pointer :: params

    integer :: i, j, stat
    integer, allocatable :: setids(:)
    type(parameter_list), pointer :: pliste => null(), plistb => null()
    character(:), allocatable :: errmsg
    type(bitfield) :: bitmask

    call set_plists

    if (associated(pliste)) then
      call pliste%get('face-set-ids', setids)
      allocate(this%is_ebc_edge(this%mesh%nedge))
      this%is_ebc_edge = .false.
      do i = 1, size(setids)
        call this%mesh%get_face_set_bitmask([setids(i)], bitmask, stat)
        ASSERT(stat == 0)
        do j = 1, this%mesh%nface
          if (popcnt(iand(bitmask, this%mesh%face_set_mask(j))) /= 0) &
              this%is_ebc_edge(this%mesh%fedge(:,j)) = .true.
        end do
      end do
      this%efield = 0 ! NB: currently hardwired to zero
    end if

    if (associated(plistb)) then
      call plistb%get('face-set-ids', setids)
      call plistb%get('bfield-polarization', this%hfield_polarization)
      if (plistb%is_scalar('bfield-vfunc')) then
        call alloc_vector_func(plistb, 'bfield-vfunc', this%hfield_vfunc, stat, errmsg)
      else
        call alloc_scalar_func(plistb, 'bfield', this%hfield_func, stat, errmsg)
      end if
      if (stat /= 0) call tls_fatal("ELECTROMAGNETICS_BC: bfield function failure " // errmsg)

      allocate(this%hbface(this%mesh%nface))
      this%hbface = 0
      do i = 1, size(setids)
        call this%mesh%get_face_set_bitmask([setids(i)], bitmask, stat)
        ASSERT(stat == 0)
        do j = 1, this%mesh%nface
          if (popcnt(iand(bitmask, this%mesh%face_set_mask(j))) /= 0) then
            if (this%hbface(j) == 0) this%hbface(j) = merge(i, -i, this%mesh%fcell(2,j) == 0)
          end if
        end do
      end do
    end if

  contains

    subroutine set_plists()

      type(parameter_list), pointer :: plist
      type(parameter_list_iterator) :: piter
      character(:), allocatable :: type

      piter = parameter_list_iterator(params, sublists_only=.true.)
      if (piter%at_end()) call tls_fatal("No ELECTROMAGNETICS_BC namelists found!")

      plist => piter%sublist()
      call plist%get('type', type)
      if (type == 'bfield') then
        plistb => plist
      else
        pliste => plist
      end if

      call piter%next
      if (piter%at_end()) return
      plist => piter%sublist()
      call plist%get('type', type)
      if (type == 'bfield') then
        if (associated(plistb)) call tls_fatal("ELECTROMAGNETICS_BC: Provided two E-type BCs")
        plistb => plist
      else
        if (associated(pliste)) call tls_fatal("ELECTROMAGNETICS_BC: Provided two B-type BCs")
        pliste => plist
      end if

    end subroutine set_plists

  end subroutine init_custom


  subroutine compute(this, t)
    class(em_bc), intent(inout) :: this
    real(r8), intent(in) :: t
    ASSERT(allocated(this%is_ebc_edge))
    if (this%use_custom_bcs) then
      call this%compute_hb_custom(t)
    else
      call efield_bc%set_Eb_values(coef=[0.0_r8], e=this%efield)
      call hfield_bc%get_Hb_source(coef=[1.0_r8, 1.0_r8], bsrc=this%hsource)
    end if
    where (this%is_ebc_edge) this%hsource = 0
    ASSERT(all(ieee_is_finite(this%efield)))
    ASSERT(all(ieee_is_finite(this%hsource)))
  end subroutine compute


  subroutine compute_hb_custom(this, t)

    class(em_bc), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: j, e(3)
    real(r8) :: c, h(3)

    ASSERT(allocated(this%hbface))

    this%hsource = 0
    do j = 1, this%mesh%nface
      if (this%hbface(j) == 0) cycle
      e = this%mesh%fedge(:,j)
      c = sign(1, this%hbface(j)) / 6.0_r8
      h = this%face_edge_projection(j, t)
      this%hsource(e(1)) = this%hsource(e(1)) + c * (h(3) - h(2))
      this%hsource(e(2)) = this%hsource(e(2)) - c * (h(1) - h(3))
      this%hsource(e(3)) = this%hsource(e(3)) + c * (h(2) - h(1))
    end do
    ! this%hsource = 2 * this%hsource

  end subroutine compute_hb_custom


  function face_edge_projection(this, f, t) result(p)

    class(em_bc), intent(in) :: this
    integer, intent(in) :: f
    real(r8), intent(in) :: t
    real(r8) :: p(3)

    integer :: k
    real(r8) :: fv(3,3), args(0:3), tangent(3)

    associate (x => this%mesh%x(:,this%mesh%fnode(:,f)))
      !! Length-weighted edge projections using trapezoid rule
      args(0) = t
      do k = 1, 3
        args(1:) = x(:,k)
        fv(:,k) = this%hfield(args)
      end do
      p(1) = dot_product(x(:,3)-x(:,2), fv(:,3)+fv(:,2)) / 2
      p(2) = dot_product(x(:,1)-x(:,3), fv(:,1)+fv(:,3)) / 2
      p(3) = dot_product(x(:,2)-x(:,1), fv(:,2)+fv(:,1)) / 2
    end associate

  end function face_edge_projection

  function hfield(this, args)
    class(em_bc), intent(in) :: this
    real(r8), intent(in) :: args(:)
    real(r8) :: hfield(3)
    if (allocated(this%hfield_vfunc)) then
      hfield = this%hfield_vfunc%eval(args)
    else
      hfield = this%hfield_func%eval(args) * this%hfield_polarization
    end if
  end function hfield

end module em_bc_type

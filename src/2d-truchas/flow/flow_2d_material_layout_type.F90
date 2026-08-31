!!
!! FLOW_2D_MATERIAL_LAYOUT_TYPE
!!
!! Defines the flow-facing reduction of simulation material composition. Real
!! fluids occupy distinct slots, followed by optional VOID and lumped SOLID
!! slots. The material-slot order is fixed; geometric reconstruction uses the
!! separate full PRIORITY permutation.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module flow_2d_material_layout_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use material_model_type
  use material_class
  use material_composition_type
  implicit none
  private

  type, public :: flow_2d_material_layout
    private
    integer, allocatable :: fluid_material_ids(:), solid_material_ids(:)
    character(:), allocatable :: fluid_material_names(:)
    integer :: void_material_id = 0
    integer, allocatable :: priority(:)
  contains
    procedure :: init
    procedure :: set_priority
    procedure :: num_real_fluid
    procedure :: num_fluid
    procedure :: num_material
    procedure :: get_real_fluid_material_ids
    procedure :: get_priority
    procedure :: get_reduced_volume_fractions
    procedure :: put_reduced_volume_fractions
  end type

contains

  subroutine init(this, matl_model, stat, errmsg)

    class(flow_2d_material_layout), intent(out) :: this
    type(material_model), intent(in) :: matl_model
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, j, mid, nchar, nsolid
    class(material), pointer :: matl

    allocate(this%fluid_material_ids(matl_model%nmatl_real), this%solid_material_ids(matl_model%nmatl_real))
    this%fluid_material_ids = 0
    this%solid_material_ids = 0
    j = 0
    nsolid = 0
    do mid = 1, matl_model%nmatl_real
      if (matl_model%num_matl_phase(mid) /= 1) then
        stat = 1
        errmsg = 'flow volume tracking currently requires single-phase materials'
        return
      end if
      call matl_model%get_matl_ref(mid, matl)
      if (matl%has_attr('fluid')) then
        j = j + 1
        this%fluid_material_ids(j) = mid
      else
        nsolid = nsolid + 1
        this%solid_material_ids(nsolid) = mid
      end if
    end do
    this%fluid_material_ids = this%fluid_material_ids(:j)
    this%solid_material_ids = this%solid_material_ids(:nsolid)
    nchar = 1
    do i = 1, size(this%fluid_material_ids)
      nchar = max(nchar, len(matl_model%matl_name(this%fluid_material_ids(i))))
    end do
    allocate(character(nchar) :: this%fluid_material_names(size(this%fluid_material_ids)))
    do i = 1, size(this%fluid_material_ids)
      this%fluid_material_names(i) = matl_model%matl_name(this%fluid_material_ids(i))
    end do
    if (matl_model%have_void) this%void_material_id = matl_model%nmatl

    allocate(this%priority(this%num_material()))
    this%priority = [(i, i=1,size(this%priority))]

    stat = 0
    errmsg = ''

  end subroutine


  subroutine set_priority(this, params, stat, errmsg)

    class(flow_2d_material_layout), intent(inout) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: names(:)
    integer :: i

    ASSERT(allocated(this%priority))
    if (.not.params%is_vector('material-priority')) then
      stat = 0
      errmsg = ''
      return
    end if

    call params%get('material-priority', names, stat, errmsg)
    if (stat /= 0) return
    if (size(names) /= size(this%priority)) then
      stat = 1
      errmsg = 'material-priority must name every fluid, VOID, and SOLID slot exactly once'
      return
    end if
    do i = 1, size(names)
      this%priority(i) = slot_index(this, names(i))
      if (this%priority(i) == 0) then
        stat = 1
        errmsg = 'invalid material-priority entry: ' // names(i)
        return
      end if
    end do
    if (.not.all([(count(this%priority == i), i=1,size(this%priority))] == 1)) then
      stat = 1
      errmsg = 'material-priority entries must be unique'
      return
    end if
    stat = 0
    errmsg = ''
  end subroutine


  integer function slot_index(this, name) result(slot)

    class(flow_2d_material_layout), intent(in) :: this
    character(*), intent(in) :: name
    integer :: k

    slot = 0
    do k = 1, size(this%fluid_material_ids)
      if (name == this%fluid_material_names(k)) then
        slot = k
        return
      end if
    end do
    if (this%void_material_id /= 0 .and. name == 'VOID') then
      slot = this%num_fluid()
      return
    end if
    if (size(this%solid_material_ids) > 0 .and. name == 'SOLID') slot = this%num_material()

  end function


  integer function num_real_fluid(this)
    class(flow_2d_material_layout), intent(in) :: this

    num_real_fluid = size(this%fluid_material_ids)
  end function


  integer function num_fluid(this)
    class(flow_2d_material_layout), intent(in) :: this

    num_fluid = this%num_real_fluid() + merge(1, 0, this%void_material_id /= 0)
  end function


  integer function num_material(this)
    class(flow_2d_material_layout), intent(in) :: this

    num_material = this%num_fluid() + merge(1, 0, size(this%solid_material_ids) > 0)
  end function


  subroutine get_real_fluid_material_ids(this, material_ids)
    class(flow_2d_material_layout), intent(in) :: this
    integer, intent(out) :: material_ids(:)

    ASSERT(size(material_ids) == size(this%fluid_material_ids))
    material_ids = this%fluid_material_ids
  end subroutine


  subroutine get_priority(this, priority)
    class(flow_2d_material_layout), intent(in) :: this
    integer, intent(out) :: priority(:)

    ASSERT(size(priority) == size(this%priority))
    priority = this%priority
  end subroutine


  subroutine get_reduced_volume_fractions(this, composition, vfrac)
    class(flow_2d_material_layout), intent(in) :: this
    type(material_composition), intent(in) :: composition
    real(r8), intent(out) :: vfrac(:,:)

    integer :: m

    ASSERT(size(vfrac,1) == this%num_material())
    ASSERT(size(vfrac,2) >= size(composition%vfrac,2))
    vfrac = 0.0_r8
    do m = 1, this%num_real_fluid()
      vfrac(m,:size(composition%vfrac,2)) = composition%vfrac(this%fluid_material_ids(m),:)
    end do
    if (this%void_material_id /= 0) then
      vfrac(this%num_fluid(),:size(composition%vfrac,2)) = composition%vfrac(this%void_material_id,:)
    end if
    if (size(this%solid_material_ids) > 0) then
      do m = 1, size(this%solid_material_ids)
        vfrac(this%num_material(),:size(composition%vfrac,2)) = &
            vfrac(this%num_material(),:size(composition%vfrac,2)) + composition%vfrac(this%solid_material_ids(m),:)
      end do
    end if
  end subroutine


  subroutine put_reduced_volume_fractions(this, vfrac, composition)
    class(flow_2d_material_layout), intent(in) :: this
    real(r8), intent(in) :: vfrac(:,:)
    type(material_composition), intent(inout) :: composition

    integer :: m

    ASSERT(size(vfrac,1) == this%num_material())
    ASSERT(size(vfrac,2) >= size(composition%vfrac,2))
    do m = 1, this%num_real_fluid()
      composition%vfrac(this%fluid_material_ids(m),:) = vfrac(m,:size(composition%vfrac,2))
    end do
    if (this%void_material_id /= 0) then
      composition%vfrac(this%void_material_id,:) = vfrac(this%num_fluid(),:size(composition%vfrac,2))
    end if
  end subroutine

end module flow_2d_material_layout_type

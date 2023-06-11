#include "f90_assert.fpp"

module pbf_material_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use material_model_type
  use matl_mesh_func_type
  implicit none
  private

  type, public :: pbf_material
    private
    type(material_model), pointer :: matl_model => null() ! reference only
    integer :: m1, m2
    real(r8), allocatable :: sfrac(:) ! saved state; solid (aka non-liquid) fraction
    real(r8), allocatable :: pfrac(:) ! persistent temporary
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: update_mmf
  end type

contains

  subroutine init(this, matl_model, mmf, params, stat, errmsg)

    use parameter_list_type

    class(pbf_material), intent(out) :: this
    type(material_model), intent(in), target :: matl_model
    type(matl_mesh_func), intent(in) :: mmf
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    character(:), allocatable :: string
    real(r8), pointer :: vfrac(:,:)

    this%matl_model => matl_model

    call params%get('material1', string, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    this%m1 = matl_model%matl_index(string)
    if (this%m1 == 0) then
      stat = -1
      errmsg = 'unknown material for "material1": ' // string
      return
    else if (matl_model%num_matl_phase(this%m1) < 2) then
      stat = -1
      errmsg = '"material1" is not a multiphase material'
      return
    end if

    call params%get('material2', string, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    this%m2 = matl_model%matl_index(string)
    if (this%m2 == 0) then
      stat = -1
      errmsg = 'unknown material for "material2": ' // string
      return
    else if (this%m2 == this%m1) then
      stat = -1
      errmsg = 'same value assigned to "material1" and "material2"'
    else if (matl_model%num_matl_phase(this%m2) < 2) then
      stat = -1
      errmsg = '"material2" is not a multiphase material'
      return
    end if

    n = matl_model%num_matl_phase(this%m2)
    allocate(this%pfrac(n))

    ASSERT(mmf%num_reg() == 1)
    vfrac => mmf%reg_vol_frac(1)
    allocate(this%sfrac(size(vfrac,1)))

  end subroutine init


  subroutine set_initial_state(this, temp)
    class(pbf_material), intent(inout) :: this
    real(r8), intent(in) :: temp(:)
    integer :: j
    ASSERT(size(temp) == size(this%sfrac))
    do j = 1, size(this%sfrac)
      call this%matl_model%get_matl_phase_frac(this%m2, temp(j), this%pfrac)
      this%sfrac(j) = 1.0_r8 - this%pfrac(size(this%pfrac))
    end do
  end subroutine


  subroutine update_mmf(this, temp, mmf)

    class(pbf_material), intent(inout) :: this
    real(r8), intent(in) :: temp(:)
    type(matl_mesh_func), intent(inout) :: mmf

    integer :: j
    real(r8) :: sfrac, vfrac1
    real(r8), pointer :: vfrac(:,:)

    ASSERT(mmf%num_reg() == 1)
    vfrac => mmf%reg_vol_frac(1)
    ASSERT(size(temp) == size(vfrac,1))
    do j = 1, size(vfrac,1)
      if (vfrac(j,this%m1) == 0.0_r8) cycle ! no conversion needed
      call this%matl_model%get_matl_phase_frac(this%m2, temp(j), this%pfrac)
      sfrac = 1.0_r8 - this%pfrac(size(this%pfrac))
      if (sfrac == this%sfrac(j)) cycle ! no phase change occurred
      if (sfrac == 0.0_r8) then
        vfrac1 = 0.0_r8
      else if (sfrac < this%sfrac(j)) then
        vfrac1 = max(0.0_r8, vfrac(j,this%m1) - ((this%sfrac(j)-sfrac)/sfrac)*vfrac(j,this%m2))
      else if (sfrac > this%sfrac(j)) then
        vfrac1 = (this%sfrac(j)/sfrac)*vfrac(j,this%m1)
      end if
      vfrac(j,this%m2) = vfrac(j,this%m1) + vfrac(j,this%m2) - vfrac1
      vfrac(j,this%m1) = vfrac1
      this%sfrac(j) = sfrac
    end do

  end subroutine update_mmf

end module pbf_material_type

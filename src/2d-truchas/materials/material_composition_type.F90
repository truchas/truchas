!!
!! MATERIAL_COMPOSITION_TYPE
!!
!! This module defines MATERIAL_COMPOSITION, the owned-cell material volume
!! fraction state of a two-dimensional simulation.  The simulation owns this
!! state and provides non-owning references to physics models that need it.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module material_composition_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use material_model_type
  use region_func_type
  use vol_frac_init_procs
  use parallel_communication, only: global_any
  implicit none
  private

  public :: get_material_region_names

  type :: region_material
    integer :: material_index
  end type region_material

  type, public :: material_composition
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    type(region_material), allocatable :: region(:)
    real(r8), allocatable, public :: volume_fraction(:,:) ! (material, owned cell)
  contains
    procedure :: init
  end type material_composition

contains

  !! Return the unique material names referenced by PARAMS, in first-reference
  !! order.  PARAMS is the MATERIAL-REGIONS sublist.
  subroutine get_material_region_names(params, names, stat, errmsg)

    use parameter_list_type

    type(parameter_list), intent(inout) :: params
    character(:), allocatable, intent(out) :: names(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
    type(parameter_list_iterator) :: piter
    character(:), allocatable :: name
    character(:), allocatable :: all_names(:)
    integer :: i, j, maxlen, nregion, nname
    logical :: found

    stat = 0
    errmsg = ''
    piter = parameter_list_iterator(params, sublists_only=.true.)
    nregion = piter%count()
    if (nregion == 0) then
      stat = 1
      errmsg = 'at least one material region is required'
      return
    end if

    maxlen = 0
    do i = 1, nregion
      plist => piter%sublist()
      call plist%get('material', name, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // plist%path() // ': ' // errmsg
        return
      end if
      if (len_trim(name) == 0) then
        stat = 1
        errmsg = 'processing ' // plist%path() // ': "material" must not be empty'
        return
      end if
      maxlen = max(maxlen, len(name))
      call piter%next()
    end do

    allocate(character(maxlen) :: all_names(nregion))
    piter = parameter_list_iterator(params, sublists_only=.true.)
    nname = 0
    do i = 1, nregion
      plist => piter%sublist()
      call plist%get('material', name, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // plist%path() // ': ' // errmsg
        return
      end if
      found = .false.
      do j = 1, nname
        if (all_names(j) == name) then
          found = .true.
          exit
        end if
      end do
      if (.not.found) then
        nname = nname + 1
        all_names(nname) = name
      end if
      call piter%next()
    end do
    names = all_names(:nname)

  end subroutine get_material_region_names


  subroutine init(this, mesh, matl_model, params, rlev, stat, errmsg)

    use parameter_list_type

    class(material_composition), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(material_model), intent(in) :: matl_model
    type(parameter_list), intent(inout) :: params
    integer, intent(in) :: rlev
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(region_func) :: rfunc
    type(parameter_list), pointer :: plist
    type(parameter_list_iterator) :: piter
    character(:), allocatable :: name
    real(r8), allocatable :: region_fraction(:,:)
    integer :: i, mid

    this%mesh => mesh
    call rfunc%init(mesh, params, stat, errmsg)
    if (stat /= 0) return

    allocate(this%region(rfunc%num_region()))
    piter = parameter_list_iterator(params, sublists_only=.true.)
    do i = 1, size(this%region)
      plist => piter%sublist()
      call plist%get('material', name, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // plist%path() // ': ' // errmsg
        return
      end if
      mid = matl_model%matl_index(name)
      if (mid == 0) then
        stat = 1
        errmsg = 'processing ' // plist%path() // ': unknown material: ' // name
        return
      end if
      this%region(i)%material_index = mid
      call piter%next()
    end do

    allocate(region_fraction(size(this%region), mesh%ncell_onP))
    call compute_volume_fractions(mesh, rfunc, rlev, region_fraction, stat)
    if (stat /= 0) then
      errmsg = 'computing material-region volume fractions failed'
      return
    end if

    allocate(this%volume_fraction(matl_model%nmatl, mesh%ncell_onP), source=0.0_r8)
    do i = 1, size(this%region)
      mid = this%region(i)%material_index
      this%volume_fraction(mid,:) = this%volume_fraction(mid,:) + region_fraction(i,:)
    end do

    if (global_any(any(this%volume_fraction < 0.0_r8) .or. any(this%volume_fraction > 1.0_r8) .or. &
        any(abs(sum(this%volume_fraction, dim=1) - 1.0_r8) > 16.0_r8 * epsilon(1.0_r8)))) then
      stat = 1
      errmsg = 'invalid material volume fractions'
      return
    end if
    stat = 0
    errmsg = ''

  end subroutine init

end module material_composition_type

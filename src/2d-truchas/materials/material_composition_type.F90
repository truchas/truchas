!!
!! MATERIAL_COMPOSITION_TYPE
!!
!! This module defines the MATERIAL_COMPOSITION type, which stores the
!! on-process cell material volume fractions of a two-dimensional simulation.
!! The simulation owns this state and provides non-owning references to
!! physics models that need it.
!!
!! Neil Carlson <neil.n.carlson@gmail.com>, August 2026
!! SPDX-License-Identifier: BSD-3-Clause
!!

#include "f90_assert.fpp"

module material_composition_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_2d_mesh_type
  use material_model_type
  implicit none
  private

  public :: get_material_region_names

  type, public :: material_composition
    private
    type(unstr_2d_mesh), pointer :: mesh => null()
    real(r8), allocatable, public :: vfrac(:,:) ! (material, on-process cell)
  contains
    procedure :: init
    procedure :: init_uniform
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
    use region_func_type
    use vol_frac_init_procs
    use parallel_communication, only: global_any

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
    integer, allocatable :: region_mid(:)
    integer :: i, mid

    this%mesh => mesh
    call rfunc%init(mesh, params, stat, errmsg)
    if (stat /= 0) return

    allocate(region_mid(rfunc%num_region()))
    piter = parameter_list_iterator(params, sublists_only=.true.)
    do i = 1, size(region_mid)
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
      region_mid(i) = mid
      call piter%next()
    end do

    allocate(region_fraction(size(region_mid), mesh%ncell_onP))
    call compute_volume_fractions(mesh, rfunc, rlev, region_fraction, stat)
    if (stat /= 0) then
      errmsg = 'computing material-region volume fractions failed'
      return
    end if

    allocate(this%vfrac(matl_model%nmatl, mesh%ncell_onP), source=0.0_r8)
    do i = 1, size(region_mid)
      mid = region_mid(i)
      this%vfrac(mid,:) = this%vfrac(mid,:) + region_fraction(i,:)
    end do

    if (global_any(any(this%vfrac < 0.0_r8) .or. any(this%vfrac > 1.0_r8) .or. &
        any(abs(sum(this%vfrac, dim=1) - 1.0_r8) > 16.0_r8 * epsilon(1.0_r8)))) then
      stat = 1
      errmsg = 'invalid material volume fractions'
      return
    end if
    stat = 0
    errmsg = ''

  end subroutine init


  subroutine init_uniform(this, mesh, matl_model, material_index, stat, errmsg)

    class(material_composition), intent(out) :: this
    type(unstr_2d_mesh), target, intent(in) :: mesh
    type(material_model), intent(in) :: matl_model
    integer, intent(in) :: material_index
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    if (material_index < 1 .or. material_index > matl_model%nmatl) then
      stat = 1
      errmsg = 'invalid uniform material index'
      return
    end if
    this%mesh => mesh
    allocate(this%vfrac(matl_model%nmatl, mesh%ncell_onP), source=0.0_r8)
    this%vfrac(material_index,:) = 1.0_r8
    stat = 0
    errmsg = ''

  end subroutine init_uniform

end module material_composition_type

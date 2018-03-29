!!
!! BC_FACTORY_TYPE
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module flow_bc_factory_type

  use flow_mesh_type
  use parameter_list_type
  use truchas_logging_services
  implicit none
  private

  type, public :: flow_bc_factory
    private
    type(flow_mesh),     pointer :: mesh   => null() ! reference only - do not own
    type(parameter_list), pointer :: params => null() ! reference only - do not own
  contains
    procedure :: init
    generic :: alloc_bc => alloc_vector_bc, alloc_scalar_bc
  end type flow_bc_factory

contains

  subroutine init(this, mesh, params)
    class(flow_bc_factory), intent(out) :: this
    type(flow_mesh), intent(in), target :: mesh
    type(parameter_list), intent(in), target :: params
    this%mesh => mesh
    this%params => params
  end subroutine init

  subroutine alloc_scalar_bc(this, condition, bf)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use scalar_func_class
    use bndry_func_class
    use bndry_face_func_type
    use string_utilities, only: raise_case, i_to_c
    !use parameter_list_type !BUG: not necessary but Intel needs it (should be reported)

    class(flow_bc_factory), intent(in) :: this
    character(*), intent(in) :: condition
    class(bndry_func), allocatable, intent(out) :: bf

    type(bndry_face_func), allocatable :: bff
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: bc_params
    character(:), allocatable :: bc_type, data_func
    class(scalar_func), allocatable :: f
    integer, allocatable :: setids(:)
    integer :: seq, stat
    character(:), allocatable :: errmsg
    real(r8) :: const

!    call TLS_info('  Generating "' // trim(condition) // '" boundary condition for "' &
!                                   // trim(variable)  // '" variable')

    allocate(bff)
    call bff%init(this%mesh%mesh)

    piter = parameter_list_iterator(this%params, sublists_only=.true.)
    do while (.not.piter%at_end())
      bc_params => piter%sublist()
      call bc_params%get('condition', bc_type)
      if (raise_case(bc_type) == raise_case(condition)) then
        call bc_params%get('face-sets', setids)
        call alloc_sf(bc_params, 'data', f)
        call bff%add(f, setids, stat, errmsg)
        if (stat /= 0) call TLS_fatal('error generating boundary condition: ' // errmsg)
      end if
      call piter%next
    end do

    call bff%add_complete
    call move_alloc(bff, bf)

  end subroutine alloc_scalar_bc


  subroutine alloc_vector_bc(this, condition, bf)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use vector_func_class
    use bndry_vfunc_class
    use bndry_face_vfunc_type
    use string_utilities, only: raise_case, i_to_c
    !use parameter_list_type !BUG: not necessary but Intel needs it (should be reported)

    class(flow_bc_factory), intent(in) :: this
    character(*), intent(in) :: condition
    class(bndry_vfunc), allocatable, intent(out) :: bf

    type(bndry_face_vfunc), allocatable :: bff
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: bc_params
    character(:), allocatable :: bc_type, data_func
    class(vector_func), allocatable :: f
    integer, allocatable :: setids(:)
    integer :: seq, stat
    character(:), allocatable :: errmsg
    real(r8) :: const

!    call TLS_info('  Generating "' // trim(condition) // '" boundary condition for "' &
!                                   // trim(variable)  // '" variable')

    allocate(bff)
    call bff%init(this%mesh%mesh)

    piter = parameter_list_iterator(this%params, sublists_only=.true.)
    do while (.not.piter%at_end())
      bc_params => piter%sublist()
      call bc_params%get('condition', bc_type)
      if (raise_case(bc_type) == raise_case(condition)) then
        call bc_params%get('face-sets', setids)
        call alloc_vf(bc_params, 'data', f)
        call bff%add(f, setids, stat, errmsg)
        if (stat /= 0) call TLS_fatal('error generating boundary condition: ' // errmsg)
      end if
      call piter%next
    end do

    call bff%add_complete
    call move_alloc(bff, bf)

  end subroutine alloc_vector_bc


  !! Return the SCALAR_FUNC class function F specified by PARAM, which is either
  !! a real scalar or the name of a function stored in the table of functions.

  subroutine alloc_sf(plist, param, f)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use scalar_func_factories
    use scalar_func_table, only: lookup_func

    type(parameter_list), intent(inout) :: plist
    character(*), intent(in) :: param
    class(scalar_func), allocatable, intent(out) :: f

    integer :: stat
    real(r8) :: const
    character(:), allocatable :: name

    call plist%get(param, const, stat=stat)
    if (stat == 0) then ! it is a real constant
      call alloc_const_scalar_func(f, const)
    else  ! it is the name of a function
      call plist%get(param, name)
      call lookup_func(name, f)
    end if

  end subroutine alloc_sf


  !! Return the VECTOR_FUNC class function F specified by PARAM, which is either
  !! a real vector or the name of a function stored in the table of functions.

  subroutine alloc_vf(plist, param, f)

    use,intrinsic :: iso_fortran_env, only: r8 => real64
    use vector_func_factories
    use vector_func_table, only: lookup_func

    type(parameter_list), intent(inout) :: plist
    character(*), intent(in) :: param
    class(vector_func), allocatable, intent(out) :: f

    integer :: stat
    real(r8) :: const
    character(:), allocatable :: name

    call plist%get(param, const, stat=stat)
    if (stat == 0) then ! it is a real constant
      call alloc_const_vector_func(f, const)
    else  ! it is the name of a function
      call plist%get(param, name)
      call lookup_func(name, f)
    end if

  end subroutine alloc_vf

end module flow_bc_factory_type

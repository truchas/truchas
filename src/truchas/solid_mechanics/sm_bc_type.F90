!!
!! Zach Jibben <zjibben@lanl.gov>
!! September 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use truchas_logging_services
  use bndry_func1_class
  use bndry_face_func_type
  implicit none
  private

  type :: bndry_func1_box
    class(bndry_func1), allocatable :: p
  end type bndry_func1_box

  type, public :: sm_bc
    type(bndry_func1_box) :: displacement(3), traction(3)
    class(bndry_func1), allocatable :: displacementn, tractionn
  contains
    procedure :: init
  end type sm_bc

contains

  subroutine init(this, params, mesh, stat, errmsg)

    use parameter_list_type
    use unstr_mesh_type

    class(sm_bc), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    stat = 0
    call alloc_bc('displacement', this%displacement)
    if (stat /= 0) return
    call alloc_bc('traction', this%traction)
    if (stat /= 0) return
    ! call alloc_displacementn_bc
    ! call alloc_tractionn_bc
    ! call alloc_contact_gaps
    ! call alloc_normalconstraint_gaps
    ! call alloc_freeinterface_gaps

    ! TODO: ensure consistency

  contains

    subroutine alloc_bc(prefix, bc)

      character(*), intent(in) :: prefix
      class(bndry_func1_box), intent(out) :: bc(:)

      character(1), parameter :: dirstr(3) = ['x','y','z']
      integer :: d
      type(bndry_face_func), allocatable :: bff

      do d = 1, 3
        allocate(bff)
        call bff%init(mesh, bndry_only=.false.)
        call iterate_list(params, prefix//'-'//dirstr(d), prefix, bff, stat, errmsg)
        if (stat /= 0) return
        call bff%add_complete
        call move_alloc(bff, bc(d)%p)
      end do

    end subroutine alloc_bc

  end subroutine init


  !! This auxiliary subroutine iterates over the parameter list and for each
  !! BC sublist that matches the given TYPE, it calls the supplied subroutine
  !! PROC, which is expected to incrementally construct the BC using that data.
  subroutine iterate_list(params, type, data_label, bff, stat, errmsg)

    use scalar_func_class
    use scalar_func_factories, only: alloc_scalar_func
    use string_utilities, only: lower_case

    type(parameter_list), intent(inout) :: params
    character(*), intent(in) :: type, data_label
    type(bndry_face_func), intent(inout) :: bff
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    integer, allocatable :: setids(:)
    character(:), allocatable :: this_type
    class(scalar_func), allocatable :: f

    stat = 0
    piter = parameter_list_iterator(params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      call plist%get('type', this_type, stat=stat, errmsg=errmsg)
      if (stat /= 0) exit
      if (lower_case(this_type) == type) then  ! use this sublist
        call TLS_info('  using SM_BC[' // piter%name() // ']')
        call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
        if (stat /= 0) exit
        call alloc_scalar_func(plist, data_label, f, stat, errmsg)
        if (stat /= 0) return
        call bff%add(f, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do
    if (stat /= 0) errmsg = 'SM_BC[' // piter%name() // ']: ' // errmsg

  end subroutine iterate_list

end module sm_bc_type

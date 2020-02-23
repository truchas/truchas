!!
!! SOURCE_FACTORY_TYPE
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module source_factory_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use base_mesh_class
  use parameter_list_type
  use scalar_mesh_func_class
  implicit none
  private

  type, public :: source_factory
    private
    class(base_mesh), pointer :: mesh => null()       ! reference only
    type(parameter_list), pointer :: params => null() ! reference only
  contains
    procedure :: init
    procedure :: alloc_source
  end type

contains

  subroutine init(this, mesh, params)
    class(source_factory), intent(out) :: this
    class(base_mesh), intent(in), target :: mesh
    type(parameter_list), intent(in), target :: params
    this%mesh => mesh
    this%params => params
  end subroutine init

  subroutine alloc_source(this, smf, stat, errmsg)

    use scalar_cell_func1_type
    use scalar_func_class

    class(source_factory), intent(inout) :: this
    class(scalar_mesh_func), allocatable, intent(out) :: smf
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    type(scalar_cell_func1), allocatable :: scf
    class(scalar_func), allocatable :: f
    character(:), allocatable :: file

    piter = parameter_list_iterator(this%params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      call plist%get('data-file', file, stat=stat)
      if (stat == 0) then ! use this sublist
        call get_scalar_func(plist, 'prefactor', f, stat, errmsg)
        if (stat /= 0) return
        if (.not.allocated(scf)) then
          allocate(scf)
          call scf%init(this%mesh)
        end if
        call scf%add(file, f, stat, errmsg)
        if (stat /= 0) return
      end if
      call piter%next
    end do
    stat = 0
    if (.not.allocated(scf)) return
    call scf%assemble
    call move_alloc(scf, smf)

  end subroutine alloc_source

  !! This auxiliary subroutine gets the scalar function specified by the value
  !! of the parameter PARAM in the parameter list PLIST. The parameter value is
  !! either a real acalar or a character string that is the name of a function
  !! in the function table.

  subroutine get_scalar_func(plist, param, f, stat, errmsg)

    use scalar_func_factories
    use scalar_func_table, only: lookup_func

    type(parameter_list), intent(inout) :: plist
    character(*), intent(in) :: param
    class(scalar_func), allocatable, intent(out) :: f
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    real(r8) :: const
    character(:), allocatable :: fname

    call plist%get(param, fname, stat=stat)
    if (stat == 0) then ! name of a function
      call lookup_func(fname, f)
      if (.not.allocated(f)) then
        stat = 1
        errmsg = 'unknown function name: ' // fname
        return
      end if
    else  ! it must be a constant value (default 1)
      call plist%get(param, const, default=1.0_r8, stat=stat, errmsg=errmsg)
      if (stat /= 0) return
      call alloc_const_scalar_func(f, const)
    end if

  end subroutine

end module source_factory_type

!!
!! THERMAL_SOURCE_FACTORY_TYPE
!!
!! A derived type that implements methods for creating the abstract volumetric
!! source objects used by the heat conduction model. This factory creates
!! boundary conditions over a UNSTR_MESH object that are defined by a
!! parameter list.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! David Neill-Asanza <dhna@lanl.gov>
!! August 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module thermal_source_factory_type

  use unstr_mesh_type
  use parameter_list_type
  use scalar_mesh_func_class
  use scalar_func_class
  use scalar_func_factories, only: alloc_scalar_func
  use truchas_logging_services
  implicit none
  private

  type, public :: thermal_source_factory
    private
    type(unstr_mesh), pointer :: mesh => null() ! reference only
    type(parameter_list), pointer :: params => null() ! reference only
  contains
    procedure :: init
    procedure :: alloc_source_funcs
    procedure, private :: alloc_source_func1
    procedure, private :: alloc_source_func2
    procedure, private :: iterate_list
  end type

  !! Interface for the source call back subroutine
  abstract interface
    subroutine src_cb(plist, stat, errmsg)
      import parameter_list
      implicit none
      type(parameter_list), intent(inout) :: plist
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
  end interface

  integer, parameter :: TYPE_NONE = 0
  integer, parameter :: TYPE_SCF1 = 1
  integer, parameter :: TYPE_SCF2 = 2

contains

  subroutine init(this, mesh, params)
    class(thermal_source_factory), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(in), target :: params
    this%mesh => mesh
    this%params => params
  end subroutine init


  subroutine alloc_source_funcs(this, src, stat, errmsg)

    use scalar_mesh_multifunc_type
    use scalar_cell_func1_type
    use scalar_cell_func2_type

    class(thermal_source_factory), intent(inout) :: this
    type(scalar_mesh_multifunc), allocatable, intent(out) :: src
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(scalar_cell_func1), allocatable :: scf1
    type(scalar_cell_func2), allocatable :: scf2

    call this%alloc_source_func1(scf1, stat, errmsg)
    if (stat /= 0) return

    call this%alloc_source_func2(scf2, stat, errmsg)
    if (stat /= 0) return

    if (allocated(scf1) .or. allocated(scf2)) allocate(src)

    if (allocated(scf1)) then
      call scf1%assemble
      call move_alloc(scf1, src%f1)
    end if

    if (allocated(scf2)) then
      call scf2%assemble
      call move_alloc(scf2, src%f2)
    end if

  end subroutine alloc_source_funcs


  subroutine alloc_source_func1(this, src, stat, errmsg)

    use scalar_cell_func1_type
    use bndry_face_func_type

    class(thermal_source_factory), intent(inout) :: this
    type(scalar_cell_func1), allocatable, intent(inout) :: src
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    !TODO (08/12/2020): needs a better name than "scalar_cell_func1"?
    call TLS_info('  generating "scalar_cell_func1" thermal source')
    call this%iterate_list(TYPE_SCF1, proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(src)) call TLS_info('    none specified')

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to a scalar_cell_func1 source specification and incrementally builds the
    !! source object accordingly. NB: The SCF and MESH objects are accessed from
    !! the parent subroutine through host association.

    subroutine proc(plist, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      character(:), allocatable :: file
      call plist%get('data-file', file, stat, errmsg)
      if (stat /= 0) return
      call alloc_scalar_func(plist, 'prefactor', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(src)) then
        allocate(src)
        call src%init(this%mesh)
      end if
      call src%add(file, f, stat, errmsg)
   end subroutine proc

  end subroutine alloc_source_func1


  subroutine alloc_source_func2(this, src, stat, errmsg)

    use scalar_cell_func2_type
    use bndry_face_func_type

    class(thermal_source_factory), intent(inout) :: this
    type(scalar_cell_func2), allocatable, intent(inout) :: src
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    !TODO (08/12/2020): needs a better name than "scalar_cell_func2"?
    call TLS_info('  generating "scalar_cell_func2" thermal source')
    call this%iterate_list(TYPE_SCF2, proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(src)) call TLS_info('    none specified')

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to a scalar_cell_func2 source specification and incrementally builds the
    !! source object accordingly. NB: The SCF and MESH objects are accessed from
    !! the parent subroutine through host association.

    subroutine proc(plist, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      integer, allocatable :: setids(:)
      call plist%get('cell-set-ids', setids, stat, errmsg)
      if (stat /= 0) return
      call alloc_scalar_func(plist, 'source', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(src)) then
        allocate(src)
        call src%init(this%mesh)
      end if
      call src%add(f, setids, stat, errmsg)
   end subroutine proc

  end subroutine alloc_source_func2

  !! This auxiliary function determines the 'type' of a sublist from the
  !! parameters present (or absent) in the sublist.

  integer function sublist_type(plist)
    type(parameter_list), intent(in) :: plist
    if (plist%is_parameter('data-file')) then
      sublist_type = TYPE_SCF1
    else if (plist%is_parameter('source')) then
      sublist_type = TYPE_SCF2
    else
      sublist_type = TYPE_NONE
    end if
  end function

  !! This auxiliary subroutine iterates over the parameter list and for each
  !! SOURCE sublist of the given TYPE, it calls the supplied subroutine PROC,
  !! which is expected to incrementally construct the source using that data.

  subroutine iterate_list(this, type, proc, stat, errmsg)

    class(thermal_source_factory), intent(in) :: this
    integer, intent(in) :: type
    procedure(src_cb) :: proc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist

    stat = 0
    piter = parameter_list_iterator(this%params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      if (type == sublist_type(plist)) then  ! use this sublist
        call TLS_info('    using THERMAL_SOURCE[' // piter%name() // ']')
        call proc(plist, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do
    if (stat /= 0) errmsg = 'THERMAL_SOURCE[' // piter%name() // ']: ' // errmsg

  end subroutine iterate_list

end module thermal_source_factory_type

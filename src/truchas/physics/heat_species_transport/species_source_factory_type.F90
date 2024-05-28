!!
!! SPECIES_SOURCE_FACTORY_TYPE
!!
!! A derived type that implements methods for creating the abstract volumetric
!! source objects used by the species advection-diffusion model. This factory
!! creates mesh functions over a UNSTR_MESH object that are defined by a
!! parameter list.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module species_source_factory_type

  use unstr_base_mesh_class
  use parameter_list_type
  use scalar_mesh_func_class
  use scalar_func_class
  use scalar_func_factories, only: alloc_scalar_func
  use string_utilities, only: i_to_c
  use truchas_logging_services
  implicit none
  private

  type, public :: species_source_factory
    private
    class(unstr_base_mesh), pointer :: mesh => null() ! reference only
    type(parameter_list), pointer :: params => null() ! reference only
  contains
    procedure :: init
    procedure :: alloc_source_func
    procedure, private :: iterate_list
  end type

  !! Interface for the source call back subroutine
  abstract interface
    subroutine src_cb(plist, setids, stat, errmsg)
      import parameter_list
      implicit none
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
  end interface

contains

  subroutine init(this, mesh, params)
    class(species_source_factory), intent(out) :: this
    class(unstr_base_mesh), intent(in), target :: mesh
    type(parameter_list), intent(in), target :: params
    this%mesh => mesh
    this%params => params
  end subroutine


  subroutine alloc_source_func(this, comp, src, stat, errmsg)

    use scalar_cell_func_type

    class(species_source_factory), intent(inout) :: this
    integer, intent(in) :: comp
    class(scalar_mesh_func), allocatable, intent(out) :: src
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    character(:), allocatable :: var
    type(scalar_cell_func), allocatable :: scf

    var = 'source' // i_to_c(comp)
    
    call TLS_info('  generating source for species component ' // i_to_c(comp))
    call this%iterate_list(var, proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(scf)) call TLS_info('    none specified')

    if (allocated(scf)) then
      call scf%assemble
      call move_alloc(scf, src)
    end if

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to a scalar_cell_func source specification and incrementally builds the
    !! source object accordingly. NB: The SCF and MESH objects are accessed from
    !! the parent subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      call alloc_scalar_func(plist, var, f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(scf)) then
        allocate(scf)
        call scf%init(this%mesh)
      end if
      call scf%add(f, setids, stat, errmsg)
    end subroutine proc

  end subroutine alloc_source_func

  !! This auxiliary subroutine iterates over the parameter list and for each
  !! sublist where VAR is defined, it calls the supplied subroutine PROC,
  !! which is expected to incrementally construct the source using that data.

  subroutine iterate_list(this, var, proc, stat, errmsg)

    class(species_source_factory), intent(in) :: this
    character(*), intent(in) :: var
    procedure(src_cb) :: proc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    integer, allocatable :: setids(:)

    stat = 0
    piter = parameter_list_iterator(this%params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      if (plist%is_parameter(var)) then
        call TLS_info('    using ' // piter%name())
        call plist%get('cell-set-ids', setids, stat, errmsg)
        if (stat /= 0) exit
        call proc(plist, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do
    if (stat /= 0) errmsg = piter%name() // ': ' // errmsg

  end subroutine iterate_list

end module species_source_factory_type

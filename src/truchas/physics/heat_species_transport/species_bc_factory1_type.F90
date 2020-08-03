!!
!! SPECIES_BC_FACTORY1_TYPE
!!
!! A concrete factory class that implements methods for creating the abstract
!! boundary condition objects used by the species advection-diffusion model.
!! This factory creates boundary conditions over a UNSTR_MESH object that are
!! defined by a parameter list.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! March 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  INIT(MESH, PARAMS) initializes the factory, which will hold a reference to
!!    the specified UNSTR_MESH object MESH and parameter list PARAMS.
!!
!! The parameter list is expected to be a collection of sublist parameters
!! with each sublist defining a specific boundary condition. Any parameter
!! that is not a sublist is ignored. The names of the sublists are arbitrary
!! but should be unique to avoid possible confusion in diagnostic output.
!! Each sublist is expected to define these parameters:
!!
!!    type -- the type of boundary condition
!!    comp -- the species component to which this applies
!!    face-set-ids -- array of mesh face set IDs
!!
!! The allowed types and the additional type-specific parameters that are
!! expected are as follows:
!!
!!    "concentration"
!!      conc -- the (real) constant value of the concentration, or the name of
!!        a defined function that evaluates to the value of the concentration.
!!
!!    "flux"
!!      flux -- the (real) constant value of the concentration flux, or the
!!        name of a defined function that evaluates to the value of the
!!        concentration flux.
!!
!! Once the factory object is initialized, the following subroutines can be
!! called to create specific boundary condition objects used by the species
!! advection-diffusion model. Each subroutine has output arguments STAT and
!! ERRMSG: if an error occurs, STAT return returns a non-zero integer value,
!! and the deferred-length allocatable character ERRMSG returns an explanatory
!! message.
!!
!!  ALLOC_DIR_BC(BC, STAT, ERRMSG) allocates the BNDRY_FUNC1 class object BC
!!    that describes the species concentration Dirichlet boundary condition.
!!    Sublists with 'type' equal to 'concentration' are used to define this object.
!!
!!  ALLOC_FLUX_BC(BC, STAT, ERRMSG) allocates the BNDRY_FUNC1 class object BC
!!    that describes the species concentration flux boundary condition.
!!    Sublists with 'type' equal to 'flux' are used to define this object.
!!

module species_bc_factory1_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use species_bc_factory_class
  use unstr_mesh_type
  use parameter_list_type
  use scalar_func_class
  use scalar_func_factories, only: alloc_scalar_func
  use string_utilities, only: lower_case, i_to_c
  use truchas_logging_services
  implicit none
  private

  type, extends(species_bc_factory), public :: species_bc_factory1
    type(unstr_mesh), pointer :: mesh => null() ! reference only
    type(parameter_list), pointer :: params => null() ! reference only
  contains
    procedure :: init
    procedure :: alloc_dir_bc
    procedure :: alloc_flux_bc
    procedure, private :: iterate_list
  end type

  !! Interface for the boundary condition call back subroutine
  abstract interface
    subroutine bc_cb(plist, setids, stat, errmsg)
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

    class(species_bc_factory1), intent(out) :: this
    type(unstr_mesh), intent(in), target:: mesh
    type(parameter_list), intent(in), target :: params

    this%mesh => mesh
    this%params => params

  end subroutine init


  subroutine alloc_dir_bc(this, comp, bc, stat, errmsg)

    use bndry_func1_class
    use bndry_face_func_type

    class(species_bc_factory1), intent(inout) :: this
    integer, intent(in) :: comp
    class(bndry_func1), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_func), allocatable :: bff

    call TLS_info('Generating "concentration" species-'//i_to_c(comp)//' boundary condition')
    call this%iterate_list('concentration', comp, proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(bff)) call TLS_info('  none specified')

    if (allocated(bff)) then
      call bff%add_complete
      call move_alloc(bff, bc)
    end if

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to a concentration BC specification and incrementally builds the BC object
    !! accordingly. NB: The BFF and MESH objects are accessed from the parent
    !! subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      call alloc_scalar_func(plist, 'conc', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(f, setids, stat, errmsg)
    end subroutine proc

  end subroutine alloc_dir_bc


  subroutine alloc_flux_bc(this, comp, bc, stat, errmsg)

    use bndry_func1_class
    use bndry_face_func_type

    class(species_bc_factory1), intent(inout) :: this
    integer, intent(in) :: comp
    class(bndry_func1), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_func), allocatable :: bff

    call TLS_info('Generating "flux" species-'//i_to_c(comp)//' boundary condition')
    call this%iterate_list('flux', comp, proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(bff)) call TLS_info('  none specified')

    if (allocated(bff)) then
      call bff%add_complete
      call move_alloc(bff, bc)
    end if

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to a flux BC specification and incrementally builds the BC object
    !! accordingly. NB: The BFF and MESH objects are accessed from the parent
    !! subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      call alloc_scalar_func(plist, 'flux', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(f, setids, stat, errmsg)
    end subroutine proc

  end subroutine alloc_flux_bc

  !! This auxiliary subroutine iterates over the parameter list and for each BC
  !! sublist that matches the given TYPE, it calls the supplied subroutine
  !! PROC, which is expected to incrementally construct the BC using that data.

  subroutine iterate_list(this, type, comp, proc, stat, errmsg)

    class(species_bc_factory1), intent(in) :: this
    character(*), intent(in) :: type
    integer, intent(in) :: comp
    procedure(bc_cb) :: proc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    integer, allocatable :: setids(:)
    character(:), allocatable :: this_type
    integer :: this_comp

    stat = 0
    piter = parameter_list_iterator(this%params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      call plist%get('type', this_type, stat=stat, errmsg=errmsg)
      if (stat /= 0) exit
      call plist%get('comp', this_comp, default=1, stat=stat, errmsg=errmsg)
      if (stat /= 0) exit
      if (lower_case(this_type) == type .and. this_comp == comp) then  ! use this sublist
        call TLS_info('  using species_BC[' // piter%name() // ']')
        call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
        if (stat /= 0) exit
        call proc(plist, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do
    if (stat /= 0) errmsg = 'species_BC[' // piter%name() // ']: ' // errmsg

  end subroutine iterate_list

end module species_bc_factory1_type

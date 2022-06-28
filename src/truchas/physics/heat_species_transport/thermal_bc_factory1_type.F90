!!
!! THERMAL_BC_FACTORY1_TYPE
!!
!! A concrete factory class that implements methods for creating the abstract
!! boundary condition objects used by the heat conduction model. This factory
!! creates boundary conditions over a UNSTR_BASE_MESH object that are defined
!! by a parameter list.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! November 2018
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  INIT(MESH, SIGMA, ABSZERO, PARAMS) initializes the factory, which will hold
!!    a reference to the specified UNSTR_BASE_MESH object MESH and parameter
!!    list PARAMS. SIGMA is the value of the Stefan-Boltzmann constant, and
!!    ABSZERO the temperature of absolute zero, which are needed for the
!!    radiation BC.
!!
!! The parameter list is expected to be a collection of sublist parameters
!! with each sublist defining a specific boundary condition. Any parameter
!! that is not a sublist is ignored. The names of the sublists are arbitrary
!! but should be unique to avoid possible confusion in diagnostic output.
!! Each sublist is expected to define these parameters:
!!
!!    type -- the type of boundary condition
!!    face-set-ids -- array of mesh face set IDs
!!
!! The allowed types and the additional type-specific parameters that are
!! expected are as follows:
!!
!!    "temperature"
!!      temp -- the (real) constant value of the temperature, or the name of
!!        a defined function of (t,x,y,z).
!!
!!    "flux"
!!      flux -- the (real) constant value of the heat flux, or the name of a
!!        defined function of (t,x,y,z).
!!
!!    "oriented-flux"
!!      flux -- the (real) constant vector value of the heat flux, or the name
!!        of a defined vector function of (t,x,y,z).
!!
!!    "htc"
!!      htc -- the (real) constant heat transfer coefficient, or the name of
!!        a defined function that evaluates to the coefficient.
!!      ambient-temp -- (real) constant reference temperature, or the name of
!!        a defined function of (t,x,y,z).
!!
!!    "radiation"
!!      emissivity -- the (real) constant value of emissivity, or the name of
!!        a defined function that evaluates to the emissivity.
!!      ambient-temp -- (real) constant reference temperature, or the name of a
!!        defined function of (t).
!!
!!    "interface-htc"
!!      htc -- the (real) constant heat transfer coefficient, or the name of
!!        a defined function of (t,x,y,z).
!!
!!    "gap-radiation"
!!      emissivity -- the (real) constant value of emissivity, or the name of
!!        a defined function of (t,x,y,z).
!!
!! Once the factory object is initialized, the following subroutines can be
!! called to create specific boundary condition objects used by the heat
!! conduction model. Each subroutine has output arguments STAT and ERRMSG:
!! if an error occurs, STAT return returns a non-zero integer value, and the
!! deferred-length allocatable character ERRMSG returns an explanatory message.
!!
!!  ALLOC_DIR_BC(BC, STAT, ERRMSG) allocates the BNDRY_FUNC1 class object BC
!!    that describes the temperature Dirichlet boundary condition. Sublists
!!    with 'type' equal to 'temperature' are used to define this object.
!!
!!  ALLOC_FLUX_BC(BC, STAT, ERRMSG) allocates the BNDRY_FUNC1 class object BC
!!    that describes the heat flux boundary condition. Sublists with 'type'
!!    equal to 'flux' are used to define this object.
!!
!!  ALLOC_VFLUX_BC(BC, STAT, ERRMSG) allocates the BNDRY_VFUNC class object BC
!!    that describes the heat flux boundary condition. Sublists with 'type'
!!    equal to 'oriented-flux' are used to define this object.
!!
!!  ALLOC_HTC_BC(BC, STAT, ERRMSG) allocates the BNDRY_FUNC2 class object BC
!!    that describes the HTC boundary condition. Sublists with 'type' equal
!!    to 'htc' are used to define this object.
!!
!!  ALLOC_RAD_BC(BC, STAT, ERRMSG) allocates the BNDRY_FUNC2 class object BC
!!    that describes the simple radiation boundary condition. Sublists with
!!    'type' equal to 'radiation' are used to define this object.
!!
!!  ALLOC_HTC_IC(IC, STAT, ERRMSG) allocates the INTFC_FUNC2 class object IC
!!    that describes the interface HTC interface condition. Sublists with
!!    'type' equal to 'interface-htc' are used to define this object.
!!
!!  ALLOC_RAD_IC(IC, STAT, ERRMSG) allocates the INTFC_FUNC2 class object IC
!!    that describes the gap radiation interface condition. Sublists with
!!    'type' equal to 'gap-radiation' are used to define this object.
!!

module thermal_bc_factory1_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use thermal_bc_factory_class
  use unstr_base_mesh_class
  use parameter_list_type
  use scalar_func_class
  use scalar_func_factories, only: alloc_scalar_func
  use string_utilities, only: lower_case
  use truchas_logging_services
  implicit none
  private

  type, extends(thermal_bc_factory), public :: thermal_bc_factory1
    class(unstr_base_mesh), pointer :: mesh => null() ! reference only
    type(parameter_list), pointer :: params => null() ! reference only
    real(r8) :: sigma, abszero
  contains
    procedure :: init
    procedure :: alloc_dir_bc
    procedure :: alloc_flux_bc
    procedure :: alloc_vflux_bc
    procedure :: alloc_htc_bc
    procedure :: alloc_rad_bc
    procedure :: alloc_htc_ic
    procedure :: alloc_rad_ic
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

  subroutine init(this, mesh, sigma, abszero, params)

    class(thermal_bc_factory1), intent(out) :: this
    class(unstr_base_mesh), intent(in), target:: mesh
    real(r8), intent(in) :: sigma, abszero
    type(parameter_list), intent(in), target :: params

    this%mesh => mesh
    this%sigma = sigma
    this%abszero = abszero
    this%params => params

  end subroutine init


  subroutine alloc_dir_bc(this, bc, stat, errmsg)

    use bndry_func1_class
    use bndry_face_func_type

    class(thermal_bc_factory1), intent(inout) :: this
    class(bndry_func1), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_func), allocatable :: bff

    call TLS_info('  generating "temperature" thermal boundary condition')
    call this%iterate_list('temperature', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(bff)) call TLS_info('    none specified')

    if (allocated(bff)) then
      call bff%add_complete
      call move_alloc(bff, bc)
    end if

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to a temperature BC specification and incrementally builds the BC object
    !! accordingly. NB: The BFF and MESH objects are accessed from the parent
    !! subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      call alloc_scalar_func(plist, 'temp', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(f, setids, stat, errmsg)
    end subroutine proc

  end subroutine alloc_dir_bc


  subroutine alloc_flux_bc(this, bc, stat, errmsg)

    use bndry_func1_class
    use bndry_face_func_type

    class(thermal_bc_factory1), intent(inout) :: this
    class(bndry_func1), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_face_func), allocatable :: bff

    call TLS_info('  generating "flux" thermal boundary condition')
    call this%iterate_list('flux', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(bff)) call TLS_info('    none specified')

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


  subroutine alloc_vflux_bc(this, bc, stat, errmsg)
    
    use bndry_func2_class
    use vflux_bndry_func_type

    class(thermal_bc_factory1), intent(inout) :: this
    class(bndry_func2), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(vflux_bndry_func), allocatable :: bff

    call TLS_info('  generating "oriented-flux" thermal boundary condition')
    call this%iterate_list('oriented-flux', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(bff)) call TLS_info('    none specified')

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
      call alloc_scalar_func(plist, 'absorptivity', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(bff)) then
        allocate(bff)
        call bff%init(this%mesh)
      end if
      call bff%add(f, setids, stat, errmsg)
    end subroutine proc

  end subroutine alloc_vflux_bc


  subroutine alloc_htc_bc(this, bc, stat, errmsg)

    use bndry_func2_class
    use htc_bndry_func_type

    class(thermal_bc_factory1), intent(inout) :: this
    class(bndry_func2), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(htc_bndry_func), allocatable :: htc

    call TLS_info('  generating "htc" thermal boundary condition')
    call this%iterate_list('htc', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(htc)) call TLS_info('    none specified')

    if (allocated(htc)) then
      call htc%add_complete
      call move_alloc(htc, bc)
    end if

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to an HTC BC specification and incrementally builds the BC object
    !! accordingly. NB: The HTC and MESH objects are accessed from the parent
    !! subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f1, f2
      call alloc_scalar_func(plist, 'htc', f1, stat, errmsg)
      if (stat /= 0) return
      call alloc_scalar_func(plist, 'ambient-temp', f2, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(htc)) then
        allocate(htc)
        call htc%init(this%mesh)
      end if
      call htc%add(f1, f2, setids, stat, errmsg)
    end subroutine proc

  end subroutine alloc_htc_bc


  subroutine alloc_rad_bc(this, bc, stat, errmsg)

    use bndry_func2_class
    use rad_bndry_func_type

    class(thermal_bc_factory1), intent(inout) :: this
    class(bndry_func2), allocatable, intent(out) :: bc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(rad_bndry_func), allocatable :: rad

    call TLS_info('  generating "radiation" thermal boundary condition')
    call this%iterate_list('radiation', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(rad)) call TLS_info('    none specified')

    if (allocated(rad)) then
      call rad%add_complete
      call move_alloc(rad, bc)
    end if

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to a radiation BC specification and incrementally builds the BC object
    !! accordingly. NB: The RAD and MESH objects are accessed from the parent
    !! subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f1, f2
      call alloc_scalar_func(plist, 'emissivity', f1, stat, errmsg)
      if (stat /= 0) return
      call alloc_scalar_func(plist, 'ambient-temp', f2, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(rad)) then
        allocate(rad)
        call rad%init(this%mesh, this%sigma, this%abszero)
      end if
      call rad%add(f1, f2, setids, stat, errmsg)
    end subroutine proc

  end subroutine alloc_rad_bc


  subroutine alloc_htc_ic(this, ic, stat, errmsg)

    use intfc_func2_class
    use htc_intfc_func_type

    class(thermal_bc_factory1), intent(inout) :: this
    class(intfc_func2), allocatable, intent(out) :: ic
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(htc_intfc_func), allocatable :: htc

    call TLS_info('  generating "interface-htc" thermal interface condition')
    call this%iterate_list('interface-htc', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(htc)) call TLS_info('    none specified')

    if (allocated(htc)) then
      call htc%add_complete
      call move_alloc(htc, ic)
    end if

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to an interface HTC IC specification and incrementally builds the IC object
    !! accordingly. NB: The HTC and MESH objects are accessed from the parent
    !! subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      call alloc_scalar_func(plist, 'htc', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(htc)) then
        allocate(htc)
        call htc%init(this%mesh)
      end if
      call htc%add(f, setids, stat, errmsg)
    end subroutine proc

  end subroutine alloc_htc_ic


  subroutine alloc_rad_ic(this, ic, stat, errmsg)

    use intfc_func2_class
    use rad_intfc_func_type

    class(thermal_bc_factory1), intent(inout) :: this
    class(intfc_func2), allocatable, intent(out) :: ic
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(rad_intfc_func), allocatable :: rad

    call TLS_info('  generating "gap-radiation" thermal interface condition')
    call this%iterate_list('gap-radiation', proc, stat, errmsg)
    if (stat /= 0) return
    if (.not.allocated(rad)) call TLS_info('    none specified')

    if (allocated(rad)) then
      call rad%add_complete
      call move_alloc(rad, ic)
    end if

  contains

    !! This call-back subroutine processes parameter list data that is specific
    !! to a gap radiation IC specification and incrementally builds the IC object
    !! accordingly. NB: The rad and MESH objects are accessed from the parent
    !! subroutine through host association.

    subroutine proc(plist, setids, stat, errmsg)
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(scalar_func), allocatable :: f
      call alloc_scalar_func(plist, 'emissivity', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(rad)) then
        allocate(rad)
        call rad%init(this%mesh, this%sigma, this%abszero)
      end if
      call rad%add(f, setids, stat, errmsg)
    end subroutine proc

  end subroutine alloc_rad_ic

  !! This auxiliary subroutine iterates over the parameter list and for each BC
  !! sublist that matches the given TYPE, it calls the supplied subroutine
  !! PROC, which is expected to incrementally construct the BC using that data.

  subroutine iterate_list(this, type, proc, stat, errmsg)

    class(thermal_bc_factory1), intent(in) :: this
    character(*), intent(in) :: type
    procedure(bc_cb) :: proc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    integer, allocatable :: setids(:)
    character(:), allocatable :: this_type

    stat = 0
    piter = parameter_list_iterator(this%params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      call plist%get('type', this_type, stat=stat, errmsg=errmsg)
      if (stat /= 0) exit
      if (lower_case(this_type) == lower_case(type)) then  ! use this sublist
        call TLS_info('    using THERMAL_BC[' // piter%name() // ']')
        call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
        if (stat /= 0) exit
        call proc(plist, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do
    if (stat /= 0) errmsg = 'THERMAL_BC[' // piter%name() // ']: ' // errmsg

  end subroutine iterate_list

end module thermal_bc_factory1_type

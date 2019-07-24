!!
!! PROBES_TYPE
!!
!! This module defines the derived type PROBES and supporting abstract base
!! classes PROBE_FIELD and PROBE_FIELD_FACTORY which implement a generic
!! framework for solution probes that is independent of application details.
!!
!! Michael Hall <hall@lanl.gov>
!! Neil N. Carlson <nnc@lanl.gov>
!! July 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The solution probe framework is implemented using three classes:
!!
!!    * The abstract base class PROBE_FIELD defines the interface that will
!!      be used to access data from cell or node-based solution fields. This
!!      isolates probes from the details of where the data is coming from.
!!      Clients will extend this class with concrete implementations tailored
!!      to their application. The base class contains two components, KIND
!!      which indicates whether the field is cell or node based, and LABEL,
!!      which is a string used to label the data column(s) in the output file.
!!      It also has a single type bound function VALUE which returns the value
!!      of the solution field on a given cell or node.
!!
!!    * The abstract factory class PROBE_FIELD_FACTORY defines the interface
!!      ALLOC_PROBE_FIELD that will be used to create PROBE_FIELD class objects
!!      as directed by the input. Clients will extend this class with a concrete
!!      factory tailored to their application. ALLOC_PROBE_FIELD takes a data
!!      name as input and allocates a PROBE_FIELD class object of the corresponding
!!      type
!!
!!    * The derived type PROBES encapsulates a collection of solution probes.
!!      It has no public components and only two type bound subroutines: INIT
!!      that initializes the object, and WRITE that writes the probe data to
!!      the output file.
!!

#include "f90_assert.fpp"

module probes_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use parameter_list_type
  implicit none
  private

  !! Abstract class that defines the interface probes will use for getting
  !! field data. Clients will define concrete implementations of this class.
  type, abstract, public :: probe_field
    character(:), allocatable :: label    ! column label for output file
    character(:), allocatable :: kind     ! cell/node toggle
  contains
    procedure(value), deferred :: value
  end type probe_field

  abstract interface
    function value(this, index)
      import probe_field, r8
      class(probe_field), intent(in) :: this
      integer, intent(in) :: index
      real(r8), allocatable :: value(:)
    end function
  end interface

  !! Abstract factory that probes will use to instantiate PROBE_FIELD class
  !! objects. Clients will pass an instance of an extension they define.
  type, abstract, public :: probe_field_factory
  contains
    procedure(alloc_field), deferred :: alloc_field
  end type probe_field_factory

  abstract interface
    subroutine alloc_field(this, field, data)
      import probe_field_factory, probe_field
      class(probe_field_factory), intent(in) :: this
      class(probe_field), allocatable, intent(out) :: field
      character(*), intent(in) :: data
    end subroutine
  end interface

  type :: point_probe
    integer :: lun    ! logical unit for output
    integer :: index  ! probe cell/node index
    class(probe_field), allocatable :: field
  contains
    procedure :: init  => point_probe_init
    procedure :: write => point_probe_write
  end type point_probe

  type, public :: probes
    private
    type(point_probe), allocatable :: probe(:)
  contains
    procedure :: init  => probes_init
    procedure :: write => probes_write
  end type probes

contains

  subroutine probes_init(this, mesh, pf_fac, outdir, params, stat, errmsg)

    class(probes), intent(out) :: this
    type(unstr_mesh), intent(in) :: mesh
    class(probe_field_factory), intent(in) :: pf_fac
    character(*), intent(in) :: outdir
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n
    type(parameter_list), pointer :: plist
    type(parameter_list_iterator) :: piter

    piter = parameter_list_iterator(params, sublists_only=.true.)

    n = piter%count()
    allocate(this%probe(n))

    n = 0
    do while (.not.piter%at_end())
      n = n + 1
      plist => piter%sublist()
      call this%probe(n)%init(mesh, pf_fac, outdir, plist, stat, errmsg)
      if (stat /= 0) return
      call piter%next
    end do
    stat = 0

  end subroutine

  subroutine probes_write(this, time)
    class(probes), intent(in) :: this
    real(r8), intent(in) :: time
    integer :: n
    do n = 1, size(this%probe)
      call this%probe(n)%write(time)
    end do
  end subroutine probes_write

  subroutine point_probe_init(this, mesh, pf_fac, outdir, params, stat, errmsg)

    use parallel_communication, only: global_maxval

    class(point_probe), intent(out) :: this
    type(unstr_mesh), intent(in) :: mesh
    class(probe_field_factory), intent(in) :: pf_fac
    character(*), intent(in) :: outdir  ! with trailing /
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: ios
    character(:), allocatable :: string, data_file, outfile
    real(r8), allocatable :: coord(:)
    real(r8) :: s

    !! Allocate the probe field.
    call params%get('data', string, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    call pf_fac%alloc_field(this%field, string)
    if (.not.allocated(this%field)) then
      stat = 1
      errmsg = 'probe data type "' // string // '" not available'
      return
    end if

    !! Identify the cell/node closest to the probe location. The process that
    !! owns the cell/node, signified by INDEX > 0, is responsible for output.
    call params%get('coord', coord, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (size(coord) /= 3) then
      stat = 1
      errmsg = 'coord not a 3-vector'
      return
    end if
    call params%get('coord-scale-factor', s, default=1.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (s <= 0.0_r8) then
      stat = 1
      errmsg = 'coord-scale-factor must be > 0'
      return
    end if
    coord = s * coord
    select case (this%field%kind)
    case ('cell')
      this%index = mesh%nearest_cell(coord)
    case ('node')
      this%index = mesh%nearest_node(coord)
    case default
      INSIST(.false.)
    end select

    !! Process that owns the closest cell/node opens the output file.
    call params%get('data-file', data_file, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    outfile = outdir // data_file
    if (this%index > 0) then
      open(newunit=this%lun,file=outfile,status='replace',iostat=ios)
    else
      ios = 0
    end if
    ios = global_maxval(ios)
    if (ios /= 0) then
      stat = 1
      errmsg = 'error opening probe output file "' // outfile // '"'
      return
    end if

    !! Write probe output file header.
    if (this%index > 0) then
      call params%get('description', string, default='')
      write(this%lun,'(4a)') '# ', data_file, ': ', string
      write(this%lun,'(a,3(es13.6,:,","))',advance='no') '# location (', coord
      select case (this%field%kind)
      case ('cell')
        write(this%lun,'(a,i0)') ') mapped to cell ', mesh%xcell(this%index)
      case ('node')
        write(this%lun,'(a,i0)') ') mapped to node ', mesh%xnode(this%index)
      end select
      write(this%lun,'(a,a)') '# time', this%field%label
    end if

  end subroutine point_probe_init

  subroutine point_probe_write(this, time)
    class(point_probe), intent(in) :: this
    real(r8), intent(in) :: time
    if (this%index > 0) then
      write(this%lun,'(es13.6,*(es18.10))') time, this%field%value(this%index)
    end if
  end subroutine point_probe_write

end module probes_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The all_probes_module contains the all_probes_type, the probe_type, the 
! probe_field type, and a probe_field_factory type. A brief description is:
!
!   all_probes_type : a vector of probes and a pointer to the mesh
!   probe_type : an output unit, a cell/node number, and a probe field
!   probe_field : data names and labels, and a way to get the field value
!   probe_field_factory : an initialization method for a probe_field
!
! The all_probes_type and the probe_type have initialization and write 
! procedures defined here.
!
! The probe_field has a value procedure (to access the probe's value) which is
! defined later in the probes_driver module. The probe_field_factory has an 
! alloc_field procedure (to initialize a probe_field) which is defined later 
! in the probes_driver module.

module all_probes_module

  use kinds
  use unstr_mesh_type
  use parameter_list_type

  implicit none
  private

  ! Public types.

  public :: probe_field, probe_field_factory, all_probes_type

  ! Abstract probe_field type which will be extended in probes_driver. This
  ! contains the field data, or a way to generate the field data.

  type, abstract :: probe_field
    character(:), allocatable :: name     ! Field data type (e.g. Temperature).
    character(:), allocatable :: label    ! Column label for outfile header.
    character(:), allocatable :: kind     ! Cell / Node toggle.
  contains
    procedure(value), deferred :: value
  end type

  abstract interface
    function value (this, index)
      import r8, probe_field
      class(probe_field), intent(in) :: this
      integer, intent(in) :: index
      real(r8), allocatable :: value(:)
    end function
  end interface

  ! Abstract probe_field_factory type which will be extended in probes_driver.
  ! This provides a way to initialize a probe_field variable.

  type, abstract :: probe_field_factory
  contains
    procedure(alloc_field), deferred :: alloc_field
  end type

  abstract interface
    subroutine alloc_field (this, field, data)
      import probe_field
      import probe_field_factory
      class(probe_field_factory), intent(in) :: this
      class(probe_field), allocatable, intent(out) :: field
      character(*), intent(in) :: data
    end subroutine
  end interface

  ! Probe type which contains a single probe.

  type :: probe_type
    integer :: index                         ! Cell/node number for the probe.
    integer :: lun                           ! Output unit for the probe's data.
    class(probe_field), allocatable :: field ! Field data to be output.
  contains
    procedure :: init  => probe_init
    procedure :: write => probe_write
  end type

  ! All_probes type, which is essentially a vector of probes with pass thru 
  ! calls to the probe procedures.

  type :: all_probes_type
    type(probe_type), allocatable :: probe(:)    ! Vector of probes.
    type(unstr_mesh), pointer :: mesh => NULL()  ! Pointer to the mesh.
  contains
    procedure :: init  => all_probes_init
    procedure :: write => all_probes_write
  end type

contains

  subroutine probe_write (this, time)

    ! Input variables.

    class(probe_type), intent(in) :: this   ! The probe to be output.
    real(r8), intent(in) :: time            ! The time for the current cycle.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Output a line of data for this probe.
    ! Only output if this PE contains the probe data.

    if (this%index /= 0) then
      ! This format is necessary (specifically, the .9 part) to pass some of
      ! the more restrictive probe tests without modification of the tests.
      write(this%lun,'(es13.6,*(es18.9))') time, this%field%value(this%index)
    end if

  end subroutine

  subroutine probe_init (this, mesh, pf_fac, sublist, stat, errmsg)
    use truchas_env,            only: output_dir
    use input_utilities,        only: NULL_C

    ! Input variables.

    type(unstr_mesh), intent(in) :: mesh
    class(probe_field_factory), intent(in) :: pf_fac
    type(parameter_list), intent(inout) :: sublist

    ! Output variables.

    class(probe_type), intent(out) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    ! Internal variables.

    character(:), allocatable :: output_file     ! Output file name.
    integer :: unit                              ! Unit for output.
    character(:), allocatable :: name            ! Name for the probe.
    real(r8), allocatable :: coord(:)            ! Spatial coordinates.
    real(r8) :: coord_scale_factor               ! Coordinates scale factor.
    character(:), allocatable :: description     ! Description of the probe.
    character(:), allocatable :: data            ! Field name to be output.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Start with zero problems.

    stat = 0
    errmsg = NULL_C

    ! Set probe internal values from parameter list.

    call sublist%get ('name',               name)
    call sublist%get ('coord',              coord)
    call sublist%get ('data',               data)
    call sublist%get ('coord-scale-factor', coord_scale_factor, &
                                            default=1.0_r8)
    call sublist%get ('description',        description, &
                                            default=NULL_C)

    ! Set probe field data.

    call pf_fac%alloc_field (this%field, data)
    if (.not.allocated(this%field)) then
      stat = -1
      errmsg = 'Probe Field Data ' // trim(data) // ' not available.'
      return
    end if

    ! Scale coordinates and set cell or node closest to the result.
    ! A zero value indicates that the location is not on this PE.

    coord = coord * coord_scale_factor
    if (this%field%kind == "Cell") then
      this%index = mesh%nearest_cell(coord)
    else if (this%field%kind == "Node") then
      this%index = mesh%nearest_node(coord)
    end if

    ! Return if this PE does not contain the data.

    if (this%index == 0) return

    ! Set output file name and unit. Note that the unit is only opened
    ! on the PE that contains the probe cell or node.

    output_file = trim(output_dir) // name // '.dat'
    open (newunit=unit, file=output_file, status='replace')
    this%lun = unit

    ! Write header information.

    write (unit,'(a,a)') '# Probe name: ', name
    write (unit,'(a,a)') '# Probe description: ', description
    write (unit,'(a,1pe13.6,a,1pe13.6,a,1pe13.6,a)') &
      '# Probe coordinates: (', coord(1), ',', &
                                coord(2), ',', &
                                coord(3), ')'
    write (unit,'(a,a)') '# Probe kind: ', this%field%kind
    if (this%field%kind == "Cell") then
      write (unit,'(a,i12)') '# Cell number (local):  ', this%index
      write (unit,'(a,i12)') '# Cell number (global): ', mesh%xcell(this%index)
    else if (this%field%kind == "Node") then
      write (unit,'(a,i12)') '# Node number (local):  ', this%index
      write (unit,'(a,i12)') '# Node number (global): ', mesh%xnode(this%index)
    end if
    write (unit,'(a)')   '#'
    write (unit,'(a,a)') '# Time          ', trim(this%field%label)
    write (unit,'(a)')   '#'

  end subroutine

  subroutine all_probes_write (this, time)

    ! Input variables.

    class(all_probes_type), intent(in) :: this ! All_probes, the probe vector.
    real(r8), intent(in) :: time               ! The time for the current cycle.

    ! Internal variable.

    integer p   ! Loop variable for probes.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Pass through call to write out data from each individual probe.

    do p = 1, size(this%probe)
      call this%probe(p)%write(time)
    end do

  end subroutine

  subroutine all_probes_init (this, mesh, pf_fac, params, stat, errmsg)

    use input_utilities,        only: NULL_C

    ! Input variables.

    type(parameter_list), intent(inout) :: params   ! Probe info parameter list.
    type(unstr_mesh), intent(in) :: mesh             ! Mesh object.
    class(probe_field_factory), intent(in) :: pf_fac ! Probe field factory.

    ! Output variables.

    class(all_probes_type), intent(out) :: this  ! All_probes, the probe vector.
    integer, intent(out) :: stat                     ! Status flag.
    character(:), allocatable, intent(out) :: errmsg ! Error message.

    ! Internal variables.

    integer p                                 ! Loop variable for probes.
    integer nprobes                           ! Number of probes.
    type(parameter_list), pointer :: sublist  ! params sublist pointer. 
    type(parameter_list_iterator) :: piter    ! Iterator for parameter list.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Start with zero problems.

    stat = 0
    errmsg = NULL_C

    ! Initialize the iterator over the params parameter list.

    piter = parameter_list_iterator(params, sublists_only=.true.)

    ! Allocate all_probes to the correct size.

    nprobes = piter%count()
    allocate (this%probe(nprobes))

    ! Iterate over the probe sublists and initialize each corresponding probe.

    p = 1
    do while (.not.piter%at_end())
      sublist => piter%sublist()
      call this%probe(p)%init (mesh, pf_fac, sublist, stat, errmsg)
      if (stat /= 0) return
      call piter%next
      p = p+1
    end do

  end subroutine

end module all_probes_module

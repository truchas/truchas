!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The probes_driver module provides the public interface for the probes to the
! rest of Truchas. 
!
! It provides pass-through calls to initialize and write out the probes 
! (probes_init and probes_write). It extends the alloc_field procedure in a 
! probe_field_factory to initialize probe_field variables. It provides three 
! different methods to extend the value procedure to access Truchas field data
! in a probe_field. 

module probes_driver

  use all_probes_module
  use, intrinsic :: iso_fortran_env, only: r8 => real64
  use zone_module, only: zone
  use flow_driver, only: flow_enabled, flow_vel_cc_view

  ! Diagnostic use:
  use truchas_logging_services    ! Includes TLS_* procedures.

  implicit none
  private
  public :: probes_init, probes_write 

  ! Global module variable, the private object is stored here.

  type(all_probes_type) :: all_probes

  ! Extensions of classes for this module.

  type, extends(probe_field_factory) :: truchas_probe_field_factory
  contains
    procedure :: alloc_field
  end type

  type, extends(probe_field) :: probe_field_vector
    real(r8), pointer :: vector(:) => NULL()   ! Pointer to data for this field.
  contains
    procedure :: value => value_vector
  end type

  type, extends(probe_field) :: probe_field_array
    real(r8), pointer :: array(:,:) => NULL()  ! Pointer to data for this field.
  contains
    procedure :: value => value_array
  end type

  type, extends(probe_field) :: probe_field_zone_vc
  contains
    procedure :: value => value_zone_vc
  end type

contains

  function value_vector (this, index)

    ! Input variables.

    class(probe_field_vector), intent(in) :: this   ! probe_field as a vector.
    integer, intent(in) :: index                    ! Cell/node number.
 
    ! Output variable.

    real(r8), allocatable :: value_vector(:)        ! The vector of values.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Set value from internal pointer to vector.

    value_vector = [this%vector(index)]   ! Brackets used to handle scalar case.

  end function

  function value_array (this, index)

    ! Input variables.

    class(probe_field_array), intent(in) :: this    ! probe_field as an array.
    integer, intent(in) :: index                    ! Cell/node number.
  
    ! Output variable.

    real(r8), allocatable :: value_array(:)         ! The vector of values.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Set value from internal pointer to array.

    value_array = this%array(:,index)

  end function

  function value_zone_vc (this, index)

    ! Input variables.

    class(probe_field_zone_vc), intent(in) :: this  ! probe_field as zone%vc.
    integer, intent(in) :: index                    ! Cell/node number.
 
    ! Output variable.

    real(r8), allocatable :: value_zone_vc(:)       ! The vector of values.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Set value from zone%vc via use association.

    value_zone_vc = zone(index)%Vc

  end function

  subroutine alloc_field (this, field, data)

    ! Input variables.

    class(truchas_probe_field_factory), intent(in) :: this
    character(*), intent(in) :: data
 
    ! Output variable.

    class(probe_field), allocatable, intent(out) :: field

    ! Internal variables.

    type(probe_field_vector),  allocatable :: pf_field_vector
    type(probe_field_array),   allocatable :: pf_field_array
    type(probe_field_zone_vc), allocatable :: pf_field_zone_vc

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Toggle on type of data requested for the probe.

    select case (data)
    case ('Temperature')

      allocate(pf_field_vector)
      pf_field_vector%name   = 'Temperature'
      pf_field_vector%label  = 'Temperature'
      pf_field_vector%vector => Zone%Temp
      pf_field_vector%kind   = 'Cell'
      call move_alloc (pf_field_vector, field)

    case ('Pressure')

      allocate(pf_field_vector)
      pf_field_vector%name   = 'Pressure'
      pf_field_vector%label  = 'Pressure'
      pf_field_vector%vector => Zone%P
      pf_field_vector%kind   = 'Cell'
      call move_alloc (pf_field_vector, field)

    case ('Density')

      allocate(pf_field_vector)
      pf_field_vector%name   = 'Density'
      pf_field_vector%label  = 'Density'
      pf_field_vector%vector => Zone%Rho
      pf_field_vector%kind   = 'Cell'
      call move_alloc (pf_field_vector, field)

    case ('Velocity')

      if (flow_enabled()) then

        ! New flow method.
        allocate(pf_field_array)
        pf_field_array%name  = 'Velocity'
        pf_field_array%label = 'Velocity'
        pf_field_array%array => flow_vel_cc_view()
        pf_field_array%kind  = 'Cell'
        call move_alloc (pf_field_array, field)

      else

        ! Legacy flow method.
        allocate(pf_field_zone_vc)
        pf_field_zone_vc%name  = 'Legacy Velocity'
        pf_field_zone_vc%label = 'Legacy Velocity'
        pf_field_zone_vc%kind  = 'Cell'
        call move_alloc (pf_field_zone_vc, field)

      end if

    end select
  end subroutine

  subroutine probes_init ()
    use string_utilities, only: i_to_c
    use probe_namelist,   only: params
    use mesh_manager,     only: unstr_mesh_ptr
    use unstr_mesh_type

    type(unstr_mesh), pointer :: mesh => null() ! Reference only -- not owned.
    integer :: stat
    character(:), allocatable :: errmsg
    type(truchas_probe_field_factory) :: fac

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Set pointer to the mesh.

    mesh => unstr_mesh_ptr('MAIN')

    ! Initialize all probes.

    call all_probes%init (mesh, fac, params, stat, errmsg)
    if (stat /= 0) then
      call TLS_info ('MLH : all_probes_init, stat = ' // i_to_c(stat))
      call TLS_info ('MLH : all_probes_init, errmsg = ' // errmsg)
    end if

  end subroutine probes_init

  subroutine probes_write (time)

    ! Input variable.

    real(r8), intent(in) :: time      ! The time for the current cycle.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Write out data for all the probes.

    call all_probes%write (time)

  end subroutine probes_write

end module probes_driver

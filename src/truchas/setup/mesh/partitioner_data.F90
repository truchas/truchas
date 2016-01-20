!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module partitioner_data
   !----------------------------------------------------------------------------
   ! purpose:
   !
   !    data for determining which partitioning strategy to use
   !
   ! contains:
   !
   ! author: Bryan Lally (LANL, lally@lanl.gov)
   !----------------------------------------------------------------------------
   use mesh_parameter_module, only: ndim
   use truchas_logging_services
   implicit none
   private

   ! public procedures
   public :: SET_PARTITIONER,             &
             GET_PARTITIONER,             &
             GET_PARTITIONER_DEFAULT,     &
             SET_PROCESSOR_ARRAY,         &
             GET_PROCESSOR_ARRAY,         &
             GET_PROCESSOR_ARRAY_DEFAULT, &
             PARTITIONER_INIT

   ! valid partitioners
   integer, parameter, public :: PART_NONE      = 0
   integer, parameter, public :: PART_CARTESIAN = 1
   integer, parameter, public :: PART_CHACO     = 2
   integer, parameter, public :: PART_METIS     = 3
   integer, parameter, public :: PART_AUTOMATIC = 4
   integer, parameter, public :: PART_ERROR     = 5

   ! the selected partitioner and processor array
   integer, save :: partitioner
   integer, dimension(ndim), public, save :: processor_array

   !----------------------------------------------------------------------------

contains

   !----------------------------------------------------------------------------

   function SET_PARTITIONER (partitioner_arg)
      !-------------------------------------------------------------------------
      ! purpose:
      !
      !    set the partitioner value
      !
      !    returns 0 for success
      !    returns 1 for failure (invalid partitioner)
      !-------------------------------------------------------------------------
      use string_utilities, only: lower_case

      ! arguments
      character(*) :: partitioner_arg

      ! return value
      integer :: SET_PARTITIONER

      !-------------------------------------------------------------------------

      SET_PARTITIONER = 0

      select case (lower_case(trim(partitioner_arg)))
      case ('')
         partitioner = PART_NONE
      case ('none')
         partitioner = PART_NONE
      case ('cartesian')
         partitioner = PART_CARTESIAN
      case ('chaco')
         partitioner = PART_CHACO
      case ('metis')
         partitioner = PART_METIS
      case ('automatic')
         partitioner = PART_AUTOMATIC
      case default
         partitioner = PART_ERROR
         SET_PARTITIONER = 1
      end select

   end function SET_PARTITIONER

   !----------------------------------------------------------------------------

   function GET_PARTITIONER ()
      !-------------------------------------------------------------------------
      ! purpose:
      !
      !    return the partitioner value
      !-------------------------------------------------------------------------
      integer :: GET_PARTITIONER

      !-------------------------------------------------------------------------

      GET_PARTITIONER = partitioner

   end function GET_PARTITIONER

   !----------------------------------------------------------------------------

   function GET_PARTITIONER_DEFAULT ()
      !-------------------------------------------------------------------------
      ! purpose:
      !
      !    return the default partitioner value
      !-------------------------------------------------------------------------
      use parameter_module, only: string_len

      character(string_len) :: GET_PARTITIONER_DEFAULT

      !-------------------------------------------------------------------------

      GET_PARTITIONER_DEFAULT = 'automatic'

   end function GET_PARTITIONER_DEFAULT

   !----------------------------------------------------------------------------

   function SET_PROCESSOR_ARRAY (processor_array_arg)
      !-------------------------------------------------------------------------
      ! purpose:
      !
      !    set the processor array
      !
      !    returns 0 for success
      !    returns 1 for failure
      !-------------------------------------------------------------------------
      use mesh_parameter_module, only: ndim

      integer, dimension(ndim) :: processor_array_arg
      integer :: SET_PROCESSOR_ARRAY

      !-------------------------------------------------------------------------

      processor_array = processor_array_arg
      SET_PROCESSOR_ARRAY = 0

      if (ANY(processor_array < 0)) then
         SET_PROCESSOR_ARRAY = 1
      end if

   end function SET_PROCESSOR_ARRAY

   !----------------------------------------------------------------------------

   function GET_PROCESSOR_ARRAY ()
      !-------------------------------------------------------------------------
      ! purpose:
      !
      !    return the processor array
      !-------------------------------------------------------------------------
      use mesh_parameter_module, only: ndim

      integer, dimension(ndim) :: GET_PROCESSOR_ARRAY

      !-------------------------------------------------------------------------

      GET_PROCESSOR_ARRAY = processor_array

   end function GET_PROCESSOR_ARRAY

   !----------------------------------------------------------------------------

   function GET_PROCESSOR_ARRAY_DEFAULT ()
      !-------------------------------------------------------------------------
      ! purpose:
      !
      !    return the default processor array
      !-------------------------------------------------------------------------
      use mesh_parameter_module, only: ndim

      integer, dimension(ndim) :: GET_PROCESSOR_ARRAY_DEFAULT

      !-------------------------------------------------------------------------

      GET_PROCESSOR_ARRAY_DEFAULT = 0

   end function GET_PROCESSOR_ARRAY_DEFAULT

   !----------------------------------------------------------------------------

   subroutine PARTITIONER_INIT ()
      !-------------------------------------------------------------------------
      ! purpose:
      !
      !    check the partitioner info for consistency
      !    finish initialization
      !-------------------------------------------------------------------------
      use parameter_module,     only: nx_tot
      use mesh_parameter_module, only: ndim
      use parallel_info_module, only: p_info
      use mesh_gen_data,        only: GENERATED_MESH, partitions_total

      integer :: d
      integer :: processor_array_product
      logical :: fatal
      character(128), allocatable :: message(:)

      !-------------------------------------------------------------------------

      processor_array_product = 1
      do d = 1, ndim
         processor_array_product = processor_array_product * processor_array(d)
      end do

      ! logic to automatically select a partitioner depending on what we know about our environment
      if (partitioner == PART_AUTOMATIC) then
         if (partitions_total == 1) then
            partitioner = PART_NONE
            call TLS_info (' Automatic partitioner selection: None')
         else if (p_info%npe == processor_array_product .and. GENERATED_MESH()) then
            partitioner = PART_CARTESIAN
            call TLS_info (' Automatic partitioner selection: Cartesian')
         else
            partitioner = PART_CHACO
            call TLS_info (' Automatic partitioner selection: Chaco')
         end if
      end if

      select case (partitioner)

      case (PART_NONE)
         if (partitions_total > 1) then
            allocate(message(3))
            message(1) = 'simple block partitioner selected with multiple processes'
            message(2) = '*** this could be extremely iefficient ***'
            message(3) = 'cartesian, chaco, or metic could be a better choice'
            call TLS_warn (message)
            deallocate(message)
         end if

      case (PART_CARTESIAN)
         if (.not. GENERATED_MESH()) then
            call TLS_fatal ('PARTITIONER_INIT: cartesian partitioner incompatible with external mesh')
         end if

         if (ANY(processor_array < 1)) then
            call TLS_fatal ('PARTITIONER_INIT: cartesian partitioner selected, but an element of processor_array is less than 1')
         end if

         if (processor_array_product /= p_info%npe) then
            call TLS_fatal ('PARTITIONER_INIT: cartesian partitioner selected, and product of processor_array not equal to number of processors')
         end if

         fatal = .false.
         do d = 1, ndim
            if (MOD(nx_tot(d),processor_array(d)) /= 0) then
               allocate(message(1))
               write(message,'(a,i0)') 'processor_array does not divide evenly into number of cells in dimension ', d
               call TLS_error (message)
               deallocate(message)
               fatal = .true.
            end if
         end do
         call TLS_fatal_if_any (fatal, 'PARTITIONER_INIT: processor_array number-of-cells mismatch')

      case (PART_CHACO)
#ifndef USE_CHACO 
         call TLS_fatal ('PARTITIONER_INIT: Chaco requested in input file, but executable was built without Chaco')
#endif

      case (PART_METIS)
#ifndef USE_METIS
         call TLS_fatal ('PARTITIONER_INIT: Metis requested in input file, but executable was built without Metis')
#endif

      case default
         call TLS_fatal ('PARTITIONER_INIT: invalid partitioner - internal error')

      end select

      !-------------------------------------------------------------------------

   end subroutine PARTITIONER_INIT

   !----------------------------------------------------------------------------

end module partitioner_data

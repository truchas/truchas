!!
!! USTRUC_COMP_FACTORY
!!
!! This module provides a parameter list-driven procedure for instantiating a
!! USTRUC_COMP object that encapsulates the low level microstructure modeling.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014
!!
!! PROGRAMMING INTERFACE
!!
!!  NEW_USTRUC_COMP(N, PARAMS) returns a pointer to a new CLASS(USTRUC_COMP)
!!    object carries out the microstructure analysis. N is the number of points
!!    where the analysis will be applied, and PARAMS is a parameter list that
!!    contains the specification of the analysis to be done.
!!
!!  The parameters that are understood are documented in the source files for
!!  the concrete implementations of the USTRUC_PLUGIN class.  Currently, these
!!  are USTRUC_VEL1, USTRUC_TIME, and USTRUC_GV0.
!!
!! NOTES
!!
!!  1. This is a work in progress.  The present implementation is a minimal
!!  stand-in that merely enables all the analysis components, drawing data
!!  from a flat shared parameter list.  A more flexible implementation that
!!  allows for a finer grained control is planned (through a hierarchical
!!  parameter list, for example), but this will probably need to wait until
!!  the input file format is upgraded.
!!
!!  2. As currently implemented, the USTRUC_OBJECT cannot be instantiated until
!!  the number of points N is known.  Yet the specification of the analysis
!!  components is completely orthogonal to the size of the internal state
!!  arrays (N).  This coupling is undesirable and these two concepts ought
!!  to be separated if possible (FIXME).
!!

#include "f90_assert.fpp"

module ustruc_comp_factory

  use ustruc_comp_class
  use ustruc_core_type
  use ustruc_vel1_type
  use ustruc_time_type
  use ustruc_gv0_type
  use ustruc_gv1_type
  use parameter_list_type
  use parameter_list_json
  implicit none
  private

  public :: new_ustruc_comp

contains

  function new_ustruc_comp (n, params) result (comp)
  
    use truchas_logging_services

    integer, intent(in) :: n
    type(parameter_list) :: params
    class(ustruc_comp), pointer :: comp
    
    integer :: lun, stat
    character(:), allocatable :: filename, errmsg, gv_model_type
    type(parameter_list), pointer :: plist

    comp => new_ustruc_core(n)
    comp => new_ustruc_vel1(comp, params)
    
    if (params%is_parameter('gv-model-file')) then
      call params%get ('gv-model-file', filename)
      plist => params%sublist('gv-params')
      !FIXME: every process is reading -- better to read on one and broadcast(?)
      open(newunit=lun,file=filename,action='read',access='stream',form='unformatted')
      call parameters_from_json_stream (lun, plist, stat, errmsg)
      if (stat /= 0) call TLS_fatal ('error reading "' // filename // '": ' // errmsg)
      close(lun)
      call plist%get ('gv-model-type', gv_model_type)
      select case (gv_model_type)
      case ('gv0')
        comp => new_ustruc_gv0(comp, plist)
      case ('gv1')
        comp => new_ustruc_gv1(comp, plist)
      case default
        call TLS_fatal ('unknown value for gv-model-type: "' // gv_model_type // '"')
      end select
    else  ! default to the basic GV0
      comp => new_ustruc_gv0(comp, params)
    end if

    !NB: need to make sure the included components match up with the data
    !    that is being output from USTRUC_DRIVER:USTRUC_OUTPUT.

  end function new_ustruc_comp

end module ustruc_comp_factory

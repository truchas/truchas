!!
!! USTRUC_ANALYSIS_FACTORY
!!
!! This module provides a parameter list-driven procedure for instantiating a
!! USTRUC_ANALYSIS object that encapsulates the low level microstructure
!! modeling.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NOTES
!!
!!  1. The present implementation is a minimal one that merely enables a single
!!  analysis component (on top of the required core component). By default that
!!  is the basic GL analysis component, which will probably be sufficient for
!!  nearly all users. When the "model-file" parameter is defined, a different
!!  analysis component will be enabled instead, using input read from the file
!!  specified by that parameter. Currently that is only the LDRD component;
!!  this is where a custom analysis component would be enabled. A more flexible
!!  implementation that allows for a finer grained control is possible (through
!!  a hierarchical parameter list, for example), but this will probably need to
!!  wait until the Truchas input file format is upgraded.
!!
!!  2. As currently implemented, the USTRUC_OBJECT cannot be instantiated until
!!  the number of points N is known.  Yet the specification of the analysis
!!  components is completely orthogonal to the size of the internal state
!!  arrays (N).  This coupling is undesirable and these two concepts ought
!!  to be separated if possible (FIXME?).
!!

#include "f90_assert.fpp"

module ustruc_analysis_factory

  use ustruc_analysis_class
  use ustruc_core_type
  use ustruc_gl_type
  use ustruc_ldrd_type
  use parameter_list_type
  use parameter_list_json
  implicit none
  private

  public :: new_ustruc_analysis

contains

  function new_ustruc_analysis(n, params) result(this)

    use truchas_logging_services

    integer, intent(in) :: n
    type(parameter_list) :: params
    class(ustruc_analysis), pointer :: this

    integer :: lun, stat
    character(:), allocatable :: filename, errmsg, model_type
    type(parameter_list), pointer :: plist

    this => new_ustruc_core(n, params)

    if (params%is_parameter('model-file')) then  ! custom analysis module
      call params%get('model-file', filename)
      plist => params%sublist('model')
      !FIXME: every process is reading -- better to read on one and broadcast(?)
      open(newunit=lun,file=filename,action='read',access='stream',form='unformatted')
      call parameters_from_json_stream(lun, plist, stat, errmsg)
      if (stat /= 0) call TLS_fatal('error reading "' // filename // '": ' // errmsg)
      close(lun)
      call plist%get('model-type', model_type)
      select case (model_type)
      case ('ldrd')
        this => new_ustruc_ldrd(this, plist)
      case default
        call TLS_fatal('unknown value for model-type: "' // model_type // '"')
      end select
    else  ! default to the basic GL analysis module
      this => new_ustruc_gl(this, params)
    end if

    !NB: need to make sure the included components match up with the data
    !    that is being output from USTRUC_DRIVER:USTRUC_OUTPUT.

  end function new_ustruc_analysis

end module ustruc_analysis_factory

!!
!! EXODUS_ERRORS
!!
!! Neil N. Carlson <nnc@lanl.gov> 21 Sep 2004
!! Last revised 10 Nov 2004
!!
!! This module provides the catalog of error messages for the Exodus II mesh
!! read/write package.  It also provides a function that returns an error
!! string given an error number.
!!
!! This is a private module; application code should only use the top-level
!! module EXODUS.
!!

module exodus_errors

  implicit none
  private

  public :: exo_err_str

  !! Error codes
  integer, parameter, public :: &
    EX_BADMESH =  1, &! mesh data structure is not well-defined
    EX_EOPENW  =  2, &! failed to open exodus file for writing
    EX_EOPENR  =  3, &! failed to open exodus file for reading
    EX_ECLOSE  =  4, &! failed to close exodus file
    EX_EWGLOB  =  5, &! failed to write global parameters
    EX_EWNODE  =  6, &! failed to write coordinate data
    EX_EWEBLK  =  7, &! failed to write element block data
    EX_EWNSET  =  8, &! failed to write node set data
    EX_EWSSET  =  9, &! failed to write side set data
    EX_EWQA    = 10, &! failed to write QA record
    EX_ERGLOB  = 11, &! failed to read global data
    EX_ERNODE  = 12, &! failed to read coordinate data
    EX_EREBLK  = 13, &! failed to read element block data
    EX_ERNSET  = 14, &! failed to read node set data
    EX_ERSSET  = 15, &! failed to read side set data
    EX_ELFILE  = 16   ! no Exodus II large-model support

contains

  function exo_err_str (errno)
    integer, intent(in) :: errno
    character(len=len_trim(error_string(errno))) :: exo_err_str
    exo_err_str = trim(error_string(errno))
  end function exo_err_str

  pure function error_string (errno)
    integer, intent(in) :: errno
    character(len=64) :: error_string
    select case (errno)
    case (EX_BADMESH)
      error_string = 'mesh data structure is not well-defined'
    case (EX_EOPENW)
      error_string = 'failed to open exodus file for writing'
    case (EX_EOPENR)
      error_string = 'failed to open exodus file for reading'
    case (EX_ECLOSE)
      error_string = 'failed to close exodus file'
    case (EX_EWGLOB)
      error_string = 'failed to write global parameters'
    case (EX_EWNODE)
      error_string = 'failed to write coordinate data'
    case (EX_EWEBLK)
      error_string = 'failed to write element block data'
    case (EX_EWNSET)
      error_string = 'failed to write node set data'
    case (EX_EWSSET)
      error_string = 'failed to write side set data'
    case (EX_EWQA)
      error_string = 'failed to write QA record'
    case (EX_ERGLOB)
      error_string = 'failed to read global parameters'
    case (EX_ERNODE)
      error_string = 'failed to read coordinate data'
    case (EX_EREBLK)
      error_string = 'failed to read element block data'
    case (EX_ERNSET)
      error_string = 'failed to read node set data'
    case (EX_ERSSET)
      error_string = 'failed to read side set data'
    case (EX_ELFILE)
      error_string = 'Exodus II large-models are not supported'
    case default
      error_string = 'no such error number!'
    end select
  end function error_string

end module exodus_errors

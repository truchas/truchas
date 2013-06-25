!!
!! FGMVWRITE -- Fortran 90 interface to GMVWRITE functions.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!

#include "f90_assert.fpp"

module fgmvwrite

  implicit none
  public

  integer, parameter :: CELLDATA = 0, NODEDATA = 1

  character(len=3), parameter :: TRICELLTYPE  = 'tri'
  character(len=4), parameter :: QUADCELLTYPE = 'quad'
  character(len=3), parameter :: TETCELLTYPE  = 'tet'
  character(len=3), parameter :: HEXCELLTYPE  = 'hex'

  !! Interfaces to the external Fortran-compatible C-procedures in f_gmvwrite.c
  interface
    subroutine gmvwrite_openfile_ir_f (filenam, isize, rsize)
      character(len=*) :: filenam
      integer :: isize, rsize
      end subroutine
    subroutine gmvwrite_openfile_ir_ascii_f (filenam, isize, rsize)
      character(len=*) :: filenam
      integer :: isize, rsize
      end subroutine
    subroutine gmvwrite_closefile_f ()
      end subroutine
    subroutine gmvwrite_node_data_f (nndes, x, y, z)
      use kinds, only: r8
      integer :: nndes
      real(r8) :: x(*), y(*), z(*)
      end subroutine
    subroutine gmvwrite_nodeids_f (nodeids)
      integer :: nodeids(*)
      end subroutine
    subroutine gmvwrite_cell_header_f (ncells)
      integer :: ncells
      end subroutine
    subroutine gmvwrite_cell_type_f (cell_type, nverts, nodes)
      character(len=*) :: cell_type
      integer :: nverts, nodes(*)
      end subroutine
    subroutine gmvwrite_cellids_f (cellids)
      integer :: cellids(*)
      end subroutine
    subroutine gmvwrite_material_header_f (nmats, data_type)
      integer :: nmats, data_type
      end subroutine
    subroutine gmvwrite_material_name_f (matname)
      character(len=*) :: matname
      end subroutine
    subroutine gmvwrite_material_ids_f (matids, data_type)
      integer :: matids(*), data_type
      end subroutine
    subroutine gmvwrite_flag_header_f ()
      end subroutine gmvwrite_flag_header_f
    subroutine gmvwrite_flag_name_f (flagname, numtypes, data_type)
      character(len=*) :: flagname
      integer :: numtypes, data_type
      end subroutine
    subroutine gmvwrite_flag_subname_f (subname)
      character(len=*) :: subname
      end subroutine
    subroutine gmvwrite_flag_data_f (data_type, flag_data)
      integer :: data_type, flag_data(*)
      end subroutine
    subroutine gmvwrite_flag_endflag_f ()
      end subroutine gmvwrite_flag_endflag_f
    subroutine gmvwrite_probtime_f (ptime)
      use kinds, only: r8
      real(r8) :: ptime
      end subroutine gmvwrite_probtime_f
    subroutine gmvwrite_cycleno_f (cyclenum)
      integer :: cyclenum
      end subroutine gmvwrite_cycleno_f
    subroutine gmvwrite_variable_header_f ()
      end subroutine gmvwrite_variable_header_f
    subroutine gmvwrite_variable_name_data_f(data_type, varname, vids)
      use kinds, only: r8
      integer :: data_type
      character(len=*) :: varname
      real(r8) :: vids(*)
      end subroutine
    subroutine gmvwrite_variable_endvars_f()
      end subroutine
  end interface

end module fgmvwrite

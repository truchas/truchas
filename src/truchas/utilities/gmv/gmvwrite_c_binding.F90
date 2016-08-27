!!
!! GMVWRITE_C_BINDING
!!
!! Bindings to a subset of the GMVWRITE C library functions.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2016, modernization of fgmvwrite.F90
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gmvwrite_c_binding

  use,intrinsic :: iso_c_binding, only: c_int, c_float, c_double, c_char, c_null_char, c_ptr, c_loc
  implicit none
  private
  
  public :: gmvwrite_openfile_ir_f, gmvwrite_openfile_ir_ascii_f, gmvwrite_closefile_f
  public :: gmvwrite_node_data_f, gmvwrite_nodeids_f
  public :: gmvwrite_cell_header_f, gmvwrite_cell_type_f, gmvwrite_cellids_f
  public :: gmvwrite_material_header_f, gmvwrite_material_name_f, gmvwrite_material_ids_f
  public :: gmvwrite_flag_header_f, gmvwrite_flag_name_f, gmvwrite_flag_subname_f, &
            gmvwrite_flag_data_f, gmvwrite_flag_endflag_f
  public :: gmvwrite_probtime_f, gmvwrite_cycleno_f
  public :: gmvwrite_variable_header_f, gmvwrite_variable_name_data_f, gmvwrite_variable_endvars_f

  integer, parameter, public :: CELLDATA = 0, NODEDATA = 1

  character(len=3), parameter, public :: TRICELLTYPE  = 'tri'
  character(len=4), parameter, public :: QUADCELLTYPE = 'quad'
  character(len=3), parameter, public :: TETCELLTYPE  = 'tet'
  character(len=3), parameter, public :: HEXCELLTYPE  = 'hex'

  !! Interfaces to gmvwrite.c functions.  Those with a _f suffix are directly
  !! callable by client code.  Those with a _c suffix are wrapped by module
  !! procedures that handle some argument conversions.
  
  interface

    subroutine gmvwrite_openfile_ir_c(filenam, isize, rsize) &
        bind(c,name='gmvwrite_openfile_ir')
      import c_char, c_int
      character(kind=c_char) :: filenam(*)
      integer(c_int), value :: isize, rsize
    end subroutine

    subroutine gmvwrite_openfile_ir_ascii_c(filenam, isize, rsize) &
        bind(c,name='gmvwrite_openfile_ir_ascii')
      import c_char, c_int
      character(kind=c_char) :: filenam(*)
      integer(c_int), value :: isize, rsize
    end subroutine

    subroutine gmvwrite_closefile_f() bind(c,name='gmvwrite_closefile')
    end subroutine

    subroutine gmvwrite_node_data_c(nndes, x, y, z) bind(c,name='gmvwrite_node_data')
      import c_ptr
      type(c_ptr), value :: nndes, x, y, z
    end subroutine

    subroutine gmvwrite_nodeids_c(nodeids) bind(c,name='gmvwrite_nodeids')
      import c_ptr
      type(c_ptr), value :: nodeids
    end subroutine

    subroutine gmvwrite_cell_header_c(ncells) bind(c,name='gmvwrite_cell_header')
      import c_ptr
      type(c_ptr), value :: ncells
    end subroutine
  
    subroutine gmvwrite_cell_type_c(cell_type, nverts, nodes) &
        bind(c,name='gmvwrite_cell_type')
      import c_char, c_int, c_ptr
      character(kind=c_char) :: cell_type(*)
      integer(c_int), value :: nverts
      type(c_ptr), value :: nodes
    end subroutine

    subroutine gmvwrite_cellids_c(cellids) bind(c,name='gmvwrite_cellids')
      import c_ptr
      type(c_ptr), value :: cellids
    end subroutine

    subroutine gmvwrite_material_header_f(nmats, data_type) &
        bind(c,name='gmvwrite_material_header')
      import c_int
      integer(c_int), value :: nmats, data_type
    end subroutine
  
    subroutine gmvwrite_material_name_c(matname) bind(c,name='gmvwrite_material_name')
      import c_char
      character(kind=c_char) :: matname(*)
    end subroutine
    
    subroutine gmvwrite_material_ids_f(matids, data_type) &
        bind(c,name='gmvwrite_material_ids')
      import c_int
      integer(c_int) :: matids(*)
      integer(c_int), value :: data_type
    end subroutine
    
    subroutine gmvwrite_flag_header_f() bind(c,name='gmvwrite_flag_header')
    end subroutine

    subroutine gmvwrite_flag_name_c(flagname, numtypes, data_type) &
        bind(c,name='gmvwrite_flag_name')
      import c_char, c_int
      character(kind=c_char) :: flagname(*)
      integer(c_int), value :: numtypes, data_type
    end subroutine

    subroutine gmvwrite_flag_subname_c(subname) bind(c,name='gmvwrite_flag_subname')
      import c_char
      character(kind=c_char) :: subname(*)
    end subroutine

    subroutine gmvwrite_flag_data_f(data_type, flag_data) &
        bind(c,name='gmvwrite_flag_data')
      import c_int
      integer(c_int), value :: data_type
      integer(c_int) :: flag_data(*)
    end subroutine

    subroutine gmvwrite_flag_endflag_f() bind(c,name='gmvwrite_flag_endflag')
    end subroutine
    
    subroutine gmvwrite_probtime_c(ptime) bind(c,name='gmvwrite_probtime')
      import c_double
      real(c_double), value :: ptime
    end subroutine
    
    subroutine gmvwrite_cycleno_f(cyclenum) bind(c,name='gmvwrite_cycleno')
      import c_int
      integer(c_int), value :: cyclenum
    end subroutine
    
    subroutine gmvwrite_variable_header_f() bind(c,name='gmvwrite_variable_header')
    end subroutine

    subroutine gmvwrite_variable_name_data_c(data_type, varname, vids) &
        bind(c,name='gmvwrite_variable_name_data')
      import c_char, c_int, c_ptr
      integer, value :: data_type
      character(kind=c_char) :: varname(*)
      type(c_ptr), value :: vids
    end subroutine
    
    subroutine gmvwrite_variable_endvars_f() bind(c,name='gmvwrite_variable_endvars')
    end subroutine

  end interface

  !! Public generic procedures for some of the wrapper module procedures
  
  interface gmvwrite_node_data_f
    procedure gmvwrite_node_data_i4r4, gmvwrite_node_data_i4r8
  end interface

  interface gmvwrite_nodeids_f
    procedure gmvwrite_nodeids_i4
  end interface
  
  interface gmvwrite_cell_header_f
    procedure gmvwrite_cell_header_i4
  end interface

  interface gmvwrite_cell_type_f
    procedure gmvwrite_cell_type_i4
  end interface

  interface gmvwrite_cellids_f
    procedure gmvwrite_cellids_i4
  end interface

  interface gmvwrite_probtime_f
    procedure gmvwrite_probtime_r4, gmvwrite_probtime_c
  end interface

  interface gmvwrite_variable_name_data_f
    procedure gmvwrite_variable_name_data_r4, gmvwrite_variable_name_data_r8
  end interface

contains

  subroutine gmvwrite_openfile_ir_f(filenam, isize, rsize)
    character(*,kind=c_char), intent(in) :: filenam
    integer(c_int), intent(in) :: isize, rsize
    call gmvwrite_openfile_ir_c(trim(filenam)//c_null_char, isize, rsize)
  end subroutine gmvwrite_openfile_ir_f

  subroutine gmvwrite_openfile_ir_ascii_f(filenam, isize, rsize)
    character(*,kind=c_char), intent(in) :: filenam
    integer(c_int), intent(in) :: isize, rsize
    call gmvwrite_openfile_ir_ascii_c(trim(filenam)//c_null_char, isize, rsize)
  end subroutine gmvwrite_openfile_ir_ascii_f

  subroutine gmvwrite_node_data_i4r4(nndes, x, y, z)
    integer(c_int), intent(in), target :: nndes
    real(c_float), intent(in), target :: x(*), y(*), z(*)
    call gmvwrite_node_data_c(c_loc(nndes), c_loc(x), c_loc(y), c_loc(z))
  end subroutine gmvwrite_node_data_i4r4

  subroutine gmvwrite_node_data_i4r8(nndes, x, y, z)
    integer(c_int), intent(in), target :: nndes
    real(c_double), intent(in), target :: x(*), y(*), z(*)
    call gmvwrite_node_data_c(c_loc(nndes), c_loc(x), c_loc(y), c_loc(z))
  end subroutine gmvwrite_node_data_i4r8

  subroutine gmvwrite_nodeids_i4(nodeids)
    integer(c_int), intent(in), target :: nodeids(*)
    call gmvwrite_nodeids_c(c_loc(nodeids))
  end subroutine gmvwrite_nodeids_i4

  subroutine gmvwrite_cell_header_i4(ncells)
    integer(c_int), intent(in), target :: ncells
    call gmvwrite_cell_header_c(c_loc(ncells))
  end subroutine gmvwrite_cell_header_i4

  subroutine gmvwrite_cell_type_i4(cell_type, nverts, nodes)
    character(*,kind=c_char), intent(in) :: cell_type
    integer(c_int), intent(in) :: nverts
    integer(c_int), intent(in), target :: nodes(*)
    character(8,kind=c_char) :: string
    string = cell_type
    call gmvwrite_cell_type_c(string//c_null_char, nverts, c_loc(nodes))
  end subroutine gmvwrite_cell_type_i4

  subroutine gmvwrite_cellids_i4(cellids)
    integer(c_int), intent(in), target :: cellids(*)
    call gmvwrite_cellids_c(c_loc(cellids))
  end subroutine gmvwrite_cellids_i4

  subroutine gmvwrite_material_name_f(matname)
    character(*,kind=c_char), intent(in) :: matname
    character(32,kind=c_char) :: string
    string = matname
    call gmvwrite_material_name_c(string//c_null_char)
  end subroutine gmvwrite_material_name_f

  subroutine gmvwrite_flag_name_f(flagname, numtypes, data_type)
    character(*,kind=c_char), intent(in) :: flagname
    integer(c_int), intent(in) :: numtypes, data_type
    character(32,kind=c_char) :: string
    string = flagname
    call gmvwrite_flag_name_c(string//c_null_char, numtypes, data_type)
  end subroutine gmvwrite_flag_name_f
  
  subroutine gmvwrite_flag_subname_f(subname)
    character(*,kind=c_char), intent(in) :: subname
    character(32,kind=c_char) :: string
    string = subname
    call gmvwrite_flag_subname_c(string//c_null_char)
  end subroutine gmvwrite_flag_subname_f

  subroutine gmvwrite_probtime_r4(ptime)
    real(c_float), intent(in) :: ptime
    call gmvwrite_probtime_c(real(ptime,kind=c_double))
  end subroutine gmvwrite_probtime_r4

  subroutine gmvwrite_variable_name_data_r4(data_type, varname, vids)
    integer(c_int), intent(in) :: data_type
    character(*,kind=c_char), intent(in) :: varname
    real(c_float), intent(in), target :: vids(*)
    character(32,kind=c_char) :: string
    string = varname
    call gmvwrite_variable_name_data_c(data_type, string//c_null_char, c_loc(vids))
  end subroutine gmvwrite_variable_name_data_r4

  subroutine gmvwrite_variable_name_data_r8(data_type, varname, vids)
    integer(c_int), intent(in) :: data_type
    character(*,kind=c_char), intent(in) :: varname
    real(c_double), intent(in), target :: vids(*)
    character(32,kind=c_char) :: string
    string = varname
    call gmvwrite_variable_name_data_c(data_type, string//c_null_char, c_loc(vids))
  end subroutine gmvwrite_variable_name_data_r8

end module gmvwrite_c_binding

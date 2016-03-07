!
!  DANU
!     Fortran module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
module danu_module

    use iso_c_binding

    use danu_iface

    implicit none
    private

! ==============================================================================
! Public 
! ==============================================================================

    ! Parameters

    ! Return status
    integer, public, parameter :: Danu_SUCCESS = 1
    integer, public, parameter :: Danu_FAILURE = -1

    ! Dataset type
    integer, public, parameter ::DANU_DATASET_INT     = Z'0001'
    integer, public, parameter ::DANU_DATASET_FLOAT   = Z'0002'
    integer, public, parameter ::DANU_DATASET_DOUBLE  = Z'0003'
    integer, public, parameter ::DANU_DATASET_STRING  = Z'0004'
    integer, public, parameter ::DANU_DATASET_UNKNOWN = Z'FFFF'

    ! Mesh element order
    integer(C_INT), public, parameter :: LINE_ELEM_ORDER = 2
    integer(C_INT), public, parameter :: TRI_ELEM_ORDER  = 3
    integer(C_INT), public, parameter :: QUAD_ELEM_ORDER = 4
    integer(C_INT), public, parameter :: TET_ELEM_ORDER  = 4
    integer(C_INT), public, parameter :: HEX_ELEM_ORDER  = 8

    ! Enumerated Types

    ! Mesh types
    enum, bind(c)
          enumerator :: INVALID_MESH = -3
          enumerator :: UNSTRUCTURED_MESH
          enumerator :: STRUCTURED_MESH
    end enum  

    public INVALID_MESH
    public UNSTRUCTURED_MESH
    public STRUCTURED_MESH

    ! Element types
    enum, bind(c)
          enumerator :: INVALID_ELEM  = 0
          enumerator :: LINE_ELEM
          enumerator :: TRI_ELEM
          enumerator :: QUAD_ELEM
          enumerator :: TET_ELEM
          enumerator :: HEX_ELEM
    end enum

    public INVALID_ELEM
    public LINE_ELEM
    public TRI_ELEM
    public QUAD_ELEM
    public TET_ELEM
    public HEX_ELEM

    ! Functions
    public ::                    &
        output_file_is_valid

    ! Subroutines

    ! Basic file control 
    public ::                    &
        danu_file_create,        &
        danu_file_close,         &
        danu_file_open_rdonly,   &
        danu_file_open_rdwr

    ! Attributes 
    public ::                    &
        attribute_exists,        &
        attribute_count,         &
        attribute_list,          & 
        attribute_write,         &    
        attribute_read

    ! Basic output file control
    public ::                    &
        output_file_create,      &
        output_file_open_rdonly, &
        output_file_open_rdwr,   &
        output_file_open_append, &
        output_file_close      

    ! Group control
    public ::                    &
        group_exists,            &
        group_create,            &
        group_open,              &
        group_close 

    ! Meshes
    public ::                         &
        mesh_count,                   &
        mesh_list,                    &
        mesh_exists,                  &
        mesh_open,                    &
        mesh_add_unstructured,        &
        mesh_create_hex_unstruct,     &
        mesh_write_coordinates,       &
        mesh_read_coordinates,        &
        mesh_write_connectivity,      &
        mesh_read_connectivity,       &
        mesh_get_type,                &
        mesh_get_elementtype,         &
        mesh_get_dimension,           &
        mesh_get_nnodes,              &
        mesh_get_nelem,               &
        mesh_get_elem_order

    ! Simulations
    public ::                         &
        simulation_count,             &
        simulation_exists,            &
        simulation_list,              &
        simulation_add,               &
        simulation_open,              &
        simulation_link_mesh,         &
        simulation_open_mesh_link,    &
        simulation_mesh_link_exists
   
    ! Non-series datasets
    public ::                         &
        data_open_dataset,            &
        data_exists,                  &
        data_count,                   &
        data_type,                    &
        data_rank,                    &
        data_dimensions,              &
        data_list,                    &
        data_write,                   &
        data_read

    ! Series group
    public ::                         &
        sequence_exists,              &
        sequence_count,               &
        sequence_list,                &
        sequence_next_id,             &
        sequence_get_id

    ! Series data
    public ::                         &
        simulation_open_data,         &
        simulation_data_exists,       &
        simulation_data_count,        &
        simulation_data_list,         &
        simulation_data_write,        &
        simulation_data_read,         &
        simulation_data_type,         &
        simulation_data_dimensions,   &
        simulation_data_rank

    ! Probes
    public ::                         &
        probe_data_exists,            &
        probe_data_count,             &
        probe_data_list,              &
        probe_data_open,              &
        probe_data_dimensions,        &
        probe_create_data,            &
        probe_data_write,             &
        probe_data_read

! ==============================================================================
! Interfaces
! ==============================================================================

interface danu_file_create
    module procedure df90_file_create
end interface danu_file_create  

interface danu_file_open_rdonly
    module procedure df90_file_open_rdonly
end interface danu_file_open_rdonly  

interface danu_file_open_rdwr
    module procedure df90_file_open_rdwr
end interface danu_file_open_rdwr  

interface danu_file_close
    module procedure df90_file_close
end interface danu_file_close  

interface output_file_is_valid
    module procedure df90_output_file_is_valid_c
    module procedure df90_output_file_is_valid_id
end interface output_file_is_valid

interface output_file_create
    module procedure df90_output_file_create
end interface output_file_create
 
interface output_file_open_rdonly
    module procedure df90_output_file_open_rdonly
end interface output_file_open_rdonly
 
interface output_file_open_rdwr
    module procedure df90_output_file_open_rdwr
end interface output_file_open_rdwr

interface output_file_open_append
    module procedure df90_output_file_open_append
end interface output_file_open_append

interface output_file_close
    module procedure df90_output_file_close
end interface output_file_close

interface attribute_exists
    module procedure df90_attribute_exists
end interface

interface attribute_count
    module procedure df90_attr_count
end interface attribute_count

interface attribute_list
    module procedure df90_attr_names
end interface attribute_list

interface attribute_write
     module procedure df90_attr_write_integer
     module procedure df90_attr_write_real4
     module procedure df90_attr_write_real8
     module procedure df90_attr_write_character0
     module procedure df90_attr_write_character1
end interface attribute_write     
 
interface attribute_read
     module procedure df90_attr_read_integer
     module procedure df90_attr_read_real4
     module procedure df90_attr_read_real8
     module procedure df90_attr_read_character0
end interface attribute_read

interface group_exists
     module procedure df90_group_exists
end interface group_exists
 
interface group_create
     module procedure df90_group_create
end interface group_create
 
interface group_open
     module procedure df90_group_open
end interface group_open
 
interface group_close
     module procedure df90_group_close
end interface group_close
 
interface mesh_count
    module procedure df90_mesh_count
end interface mesh_count

interface mesh_list
    module procedure df90_mesh_list
end interface mesh_list

interface mesh_exists
    module procedure df90_mesh_exists
end interface mesh_exists    

interface mesh_open
    module procedure df90_mesh_open
end interface mesh_open    

interface mesh_add_unstructured
    module procedure df90_mesh_add_unstructured
end interface mesh_add_unstructured  

interface mesh_create_hex_unstruct
    module procedure df90_mesh_create_hex_unstruct
end interface mesh_create_hex_unstruct

interface mesh_write_coordinates
    module procedure df90_mesh_write_coordinates_3d
    module procedure df90_mesh_write_coordinates_2d
    module procedure df90_mesh_write_coordinates_1d
end interface mesh_write_coordinates

interface mesh_read_coordinates
    module procedure df90_mesh_read_coordinates_3d
    module procedure df90_mesh_read_coordinates_2d
    module procedure df90_mesh_read_coordinates_1d
    module procedure df90_mesh_read_coordinates_idx
end interface mesh_read_coordinates

interface mesh_write_connectivity
    module procedure df90_mesh_write_connectivity
end interface mesh_write_connectivity

interface mesh_read_connectivity
    module procedure df90_mesh_read_connectivity
end interface mesh_read_connectivity

interface mesh_get_type
    module procedure df90_mesh_get_type
end interface mesh_get_type  

interface mesh_get_elementtype
    module procedure df90_mesh_get_elementtype
end interface mesh_get_elementtype  

interface mesh_get_dimension
    module procedure df90_mesh_get_dimension
end interface mesh_get_dimension  

interface mesh_get_nnodes
    module procedure df90_mesh_get_nnodes
end interface mesh_get_nnodes  

interface mesh_get_nelem
    module procedure df90_mesh_get_nelem
end interface mesh_get_nelem  

interface mesh_get_elem_order
    module procedure df90_mesh_get_elem_order
end interface mesh_get_elem_order

interface simulation_add
    module procedure df90_simulation_add
end interface simulation_add

interface simulation_open
    module procedure df90_simulation_open
  end interface simulation_open 

interface simulation_count
  module procedure df90_simulation_count  
end interface simulation_count 

interface simulation_exists
  module procedure df90_simulation_exists  
end interface simulation_exists  

interface simulation_list
  module procedure df90_simulation_list  
end interface simulation_list  

interface simulation_link_mesh
  module procedure df90_simulation_link_mesh  
end interface simulation_link_mesh  

interface simulation_open_mesh_link
  module procedure df90_simulation_open_mesh_link  
end interface simulation_open_mesh_link 

interface simulation_mesh_link_exists
  module procedure df90_simulation_mesh_link_exists  
end interface simulation_mesh_link_exists 

interface data_open_dataset
  module procedure df90_data_open_dataset
end interface data_open_dataset

interface data_exists
  module procedure df90_data_exists
end interface data_exists

interface data_count
  module procedure df90_data_count
end interface data_count

interface data_type
  module procedure df90_data_type
end interface data_type

interface data_rank
  module procedure df90_data_rank
end interface data_rank

interface data_dimensions
  module procedure df90_data_dimensions
end interface data_dimensions

interface data_list
  module procedure df90_data_list
end interface data_list

interface data_write
  module procedure df90_data_write_byte_rank0
  module procedure df90_data_write_byte_rank1
  module procedure df90_data_write_byte_rank2
  module procedure df90_data_write_byte_rank3
  module procedure df90_data_write_integer_rank0
  module procedure df90_data_write_integer_rank1
  module procedure df90_data_write_integer_rank2
  module procedure df90_data_write_integer_rank3
  module procedure df90_data_write_real4_rank0
  module procedure df90_data_write_real4_rank1
  module procedure df90_data_write_real4_rank2
  module procedure df90_data_write_real4_rank3
  module procedure df90_data_write_real8_rank0
  module procedure df90_data_write_real8_rank1
  module procedure df90_data_write_real8_rank2
  module procedure df90_data_write_real8_rank3
end interface data_write

interface data_read
  module procedure df90_data_read_byte_rank0
  module procedure df90_data_read_byte_rank1
  module procedure df90_data_read_byte_rank2
  module procedure df90_data_read_byte_rank3
  module procedure df90_data_read_integer_rank0
  module procedure df90_data_read_integer_rank1
  module procedure df90_data_read_integer_rank2
  module procedure df90_data_read_integer_rank3
  module procedure df90_data_read_real4_rank0
  module procedure df90_data_read_real4_rank1
  module procedure df90_data_read_real4_rank2
  module procedure df90_data_read_real4_rank3
  module procedure df90_data_read_real8_rank0
  module procedure df90_data_read_real8_rank1
  module procedure df90_data_read_real8_rank2
  module procedure df90_data_read_real8_rank3
end interface data_read

interface sequence_exists
  module procedure df90_sequence_exists
end interface sequence_exists

interface sequence_count
  module procedure df90_sequence_count
end interface sequence_count

interface sequence_list
  module procedure df90_sequence_list
end interface sequence_list

interface sequence_next_id
  module procedure df90_sequence_get_next_id
end interface

interface sequence_get_id
  module procedure df90_sequence_get_id_c
  module procedure df90_sequence_get_id_i
end interface

interface simulation_open_data
  module procedure df90_simulation_open_data
end interface simulation_open_data

interface simulation_data_exists
  module procedure df90_simulation_data_exists
end interface simulation_data_exists

interface simulation_data_count
  module procedure df90_simulation_data_count
end interface simulation_data_count

interface simulation_data_type
  module procedure df90_simulation_data_type
end interface simulation_data_type

interface simulation_data_rank
  module procedure df90_simulation_data_rank
end interface simulation_data_rank

interface simulation_data_dimensions
  module procedure df90_simulation_data_dimensions
end interface simulation_data_dimensions

interface simulation_data_list
  module procedure df90_simulation_data_list
end interface simulation_data_list

interface simulation_data_write
  module procedure df90_sim_data_write_byte_rank0
  module procedure df90_sim_data_write_byte_rank1
  module procedure df90_sim_data_write_byte_rank2
  module procedure df90_sim_data_write_byte_rank3
  module procedure df90_sim_data_write_integer_rank0
  module procedure df90_sim_data_write_integer_rank1
  module procedure df90_sim_data_write_integer_rank2
  module procedure df90_sim_data_write_integer_rank3
  module procedure df90_sim_data_write_real4_rank0
  module procedure df90_sim_data_write_real4_rank1
  module procedure df90_sim_data_write_real4_rank2
  module procedure df90_sim_data_write_real4_rank3
  module procedure df90_sim_data_write_real8_rank0
  module procedure df90_sim_data_write_real8_rank1
  module procedure df90_sim_data_write_real8_rank2
  module procedure df90_sim_data_write_real8_rank3
end interface simulation_data_write

interface simulation_data_read
  module procedure df90_sim_data_read_byte_rank0
  module procedure df90_sim_data_read_byte_rank1
  module procedure df90_sim_data_read_byte_rank2
  module procedure df90_sim_data_read_byte_rank3
  module procedure df90_sim_data_read_integer_rank0
  module procedure df90_sim_data_read_integer_rank1
  module procedure df90_sim_data_read_integer_rank2
  module procedure df90_sim_data_read_integer_rank3
  module procedure df90_sim_data_read_real4_rank0
  module procedure df90_sim_data_read_real4_rank1
  module procedure df90_sim_data_read_real4_rank2
  module procedure df90_sim_data_read_real4_rank3
  module procedure df90_sim_data_read_real8_rank0
  module procedure df90_sim_data_read_real8_rank1
  module procedure df90_sim_data_read_real8_rank2
  module procedure df90_sim_data_read_real8_rank3
end interface simulation_data_read

interface probe_data_exists
  module procedure df90_probe_data_exists
end interface

interface probe_data_count
  module procedure df90_probe_data_count
end interface

interface probe_data_dimensions
  module procedure df90_probe_data_dimensions
end interface

interface probe_data_list
  module procedure df90_probe_data_list
end interface

interface probe_data_open
  module procedure df90_probe_data_open
end interface probe_data_open

interface probe_create_data
  module procedure df90_probe_create_data_integer0
  module procedure df90_probe_create_data_integer
  module procedure df90_probe_create_data_real40
  module procedure df90_probe_create_data_real4
  module procedure df90_probe_create_data_real80
  module procedure df90_probe_create_data_real8
end interface probe_create_data

interface probe_data_write
  module procedure df90_probe_data_write_integer0
  module procedure df90_probe_data_write_integer
  module procedure df90_probe_data_write_real40
  module procedure df90_probe_data_write_real4
  module procedure df90_probe_data_write_real80
  module procedure df90_probe_data_write_real8
end interface probe_data_write  

interface probe_data_read
  module procedure df90_probe_data_read_integer0
  module procedure df90_probe_data_read_integer
  module procedure df90_probe_data_read_real40
  module procedure df90_probe_data_read_real4
  module procedure df90_probe_data_read_real80
  module procedure df90_probe_data_read_real8
end interface probe_data_read

! ==============================================================================
contains 
! ==============================================================================

! ==============================================================================
! functions
! ==============================================================================

logical function df90_output_file_is_valid_c(filename)

! --- Calling arguments

      character(kind=C_CHAR,len=*), intent(in) :: filename

! --- Local variables
      integer(C_INT)   :: err
      type(C_PTR)         :: fid


      call output_file_open_rdonly(filename,fid,err)

      if ( err .eq. Danu_SUCCESS ) then
          df90_output_file_is_valid_c = df90_output_file_is_valid_id(fid) 
      else
          df90_output_file_is_valid_c = .false.
      end if 

      call output_file_close(fid)

end function df90_output_file_is_valid_c
    
! ------------------------------------------------------------------------------

logical function df90_output_file_is_valid_id(fid)

! --- Calling arguments

      type(C_PTR), intent(in) :: fid

! --- Local variables
      integer(C_INT) :: flag

      call output_file_is_valid_f(fid,flag)

      if ( flag .ne. 0 ) then
           df90_output_file_is_valid_id = .true.
      else
           df90_output_file_is_valid_id = .false.
      end if

end function df90_output_file_is_valid_id
    
! ------------------------------------------------------------------------------

logical function df90_attribute_exists(fid, attr_name)

! --- Calling arguments

      type(C_PTR),                  intent(in) :: fid
      character(kind=C_CHAR,len=*), intent(in) :: attr_name

! --- Local variables
      integer(C_INT) :: flag, name_len, ierr, stat

      name_len = len(attr_name)
      call danu_attr_exists_f(fid, attr_name, name_len, flag,ierr)
      call define_return_status(ierr,stat)

      if ( stat .eq. DANU_SUCCESS ) then
          if ( flag .gt. 0 ) then
              df90_attribute_exists = .true.
          else
              df90_attribute_exists = .false.
          endif
      else
          df90_attribute_exists = .false.
      endif

      
end function df90_attribute_exists

! ------------------------------------------------------------------------------

logical function df90_group_exists(fid, grp_name)

! --- Calling arguments

      type(C_PTR),                  intent(in) :: fid
      character(kind=C_CHAR,len=*), intent(in) :: grp_name

! --- Local variables
      integer(C_INT) :: flag, name_len, ierr, stat

      name_len = len(grp_name)
      call danu_group_exists_f(fid, grp_name, name_len, flag,ierr)
      call define_return_status(ierr,stat)

      if ( stat .eq. DANU_SUCCESS ) then
          if ( flag .gt. 0 ) then
              df90_group_exists = .true.
          else
              df90_group_exists = .false.
          endif
      else
          df90_group_exists = .false.
      endif

      
end function df90_group_exists

! ------------------------------------------------------------------------------

logical function df90_simulation_exists(fid, sim_name)

! --- Calling arguments

      type(C_PTR),                  intent(in) :: fid
      character(kind=C_CHAR,len=*), intent(in) :: sim_name

! --- Local variables
      integer(C_INT) :: flag, name_len, ierr, stat

      name_len = len(sim_name)
      call simulation_exists_f(fid, sim_name, name_len, flag,ierr)
      call define_return_status(ierr,stat)

      if ( stat .eq. DANU_SUCCESS ) then
          if ( flag .gt. 0 ) then
              df90_simulation_exists = .true.
          else
              df90_simulation_exists = .false.
          endif
      else
          df90_simulation_exists = .false.
      endif

      
end function df90_simulation_exists

! ------------------------------------------------------------------------------

logical function df90_data_exists(sid, data_name)

! --- Calling arguments

      type(C_PTR),                  intent(in) :: sid
      character(kind=C_CHAR,len=*), intent(in) :: data_name

! --- Local variables
      integer(C_INT) :: flag, name_len, ierr, stat

      name_len = len(data_name)
      call data_exists_f(sid, data_name, name_len, flag,ierr)
      call define_return_status(ierr,stat)

      if ( stat .eq. DANU_SUCCESS ) then
          if ( flag .gt. 0 ) then
              df90_data_exists = .true.
          else
              df90_data_exists = .false.
          endif
      else
          df90_data_exists = .false.
      endif

      
end function df90_data_exists

! ------------------------------------------------------------------------------

logical function df90_sequence_exists(sid, sname)

! --- Calling arguments

      type(C_PTR),                  intent(in) :: sid
      character(kind=C_CHAR,len=*), intent(in) :: sname

! --- Local variables
      integer(C_INT) :: flag, name_len, ierr, stat

      name_len = len(sname)
      call sequence_exists_f(sid, sname, name_len, flag,ierr)
      call define_return_status(ierr,stat)

      if ( stat .eq. DANU_SUCCESS ) then
          if ( flag .gt. 0 ) then
              df90_sequence_exists = .true.
          else
              df90_sequence_exists = .false.
          endif
      else
          df90_sequence_exists = .false.
      endif

      
end function df90_sequence_exists

! ------------------------------------------------------------------------------

logical function df90_simulation_data_exists(nid, dname)

! --- Calling arguments

      type(C_PTR),                  intent(in) :: nid
      character(kind=C_CHAR,len=*), intent(in) :: dname

! --- Local variables
      integer(C_INT) :: flag, name_len, ierr, stat

      name_len = len(dname)
      call simulation_data_exists_f(nid, dname, name_len, flag, ierr)
      call define_return_status(ierr,stat)

      if ( stat .eq. DANU_SUCCESS ) then
          if ( flag .gt. 0 ) then
              df90_simulation_data_exists = .true.
          else
              df90_simulation_data_exists = .false.
          endif
      else
          df90_simulation_data_exists = .false.
      endif

      
end function df90_simulation_data_exists

! ------------------------------------------------------------------------------

logical function df90_probe_data_exists(sid, pname)

! --- Calling arguments

      type(C_PTR),                  intent(in) :: sid
      character(kind=C_CHAR,len=*), intent(in) :: pname

! --- Local variables
      integer(C_INT) :: flag, name_len, ierr, stat

      name_len = len(pname)
      call probe_data_exists_f(sid, pname, name_len, flag, ierr)
      call define_return_status(ierr,stat)

      if ( stat .eq. DANU_SUCCESS ) then
          if ( flag .gt. 0 ) then
              df90_probe_data_exists = .true.
          else
              df90_probe_data_exists = .false.
          endif
      else
          df90_probe_data_exists = .false.
      endif

      
end function df90_probe_data_exists

! ------------------------------------------------------------------------------
! ==============================================================================
! subroutines
! ==============================================================================

subroutine df90_file_create(filename,fid,stat)

! --- calling arguments
      character(kind=C_CHAR,len=*), intent(in)  :: filename
      type(C_PTR),                     intent(out) :: fid

      integer, intent(out), optional               :: stat

! --- local variables
      integer :: name_len
      integer(C_INT) :: err
    
      name_len = len(filename)
      call danu_file_create_f(filename,name_len,fid,err)
      if (present(stat)) then
          call define_return_status(err,stat)
      end if


end subroutine df90_file_create

! ------------------------------------------------------------------------------

subroutine df90_file_open_rdonly(filename,fid,stat)

! --- calling arguments
      character(kind=C_CHAR,len=*), intent(in)  :: filename
      type(C_PTR),                     intent(out) :: fid

      integer, intent(out), optional               :: stat

! --- local variables
      integer :: name_len
      integer(C_INT) :: err
    
      name_len = len(filename)
      call danu_file_open_rdonly_f(filename,name_len,fid,err)
      if (present(stat)) then
          call define_return_status(err,stat)
      end if


end subroutine df90_file_open_rdonly

! ------------------------------------------------------------------------------

subroutine df90_file_open_rdwr(filename,fid,stat)

! --- calling arguments
      character(kind=C_CHAR,len=*), intent(in)  :: filename
      type(C_PTR),                     intent(out) :: fid

      integer, intent(out), optional               :: stat

! --- local variables
      integer :: name_len
      integer(C_INT) :: err
    
      name_len = len(filename)
      call danu_file_open_rdwr_f(filename,name_len,fid,err)
      if (present(stat)) then
          call define_return_status(err,stat)
      end if


end subroutine df90_file_open_rdwr

! ------------------------------------------------------------------------------

subroutine df90_file_close(fid,stat)

! --- calling arguments
      type(C_PTR),                     intent(inout) :: fid

      integer, intent(out), optional               :: stat

! --- local variables
      integer(C_INT) :: err
    
      call danu_file_close_f(fid,err)
      if( present(stat) ) then
          call define_return_status(err,stat)
      endif    

end subroutine df90_file_close

! ------------------------------------------------------------------------------

subroutine df90_output_file_create(filename,fid,stat)

! --- calling arguments
      character(kind=C_CHAR,len=*), intent(in)  :: filename
      type(C_PTR),                     intent(inout) :: fid

      integer, intent(out), optional               :: stat

! --- local variables
      integer(C_INT) :: err
      integer :: name_len
    
      name_len = len_trim(filename)
      call output_file_create_f(filename,name_len,fid,err)
      if( present(stat) ) then
          call define_return_status(err,stat)
      endif    


end subroutine df90_output_file_create

! ------------------------------------------------------------------------------

subroutine df90_output_file_open_rdonly(filename,fid,stat)

! --- calling arguments
      character(kind=C_CHAR,len=*), intent(in)  :: filename
      type(C_PTR),                     intent(inout) :: fid

      integer, intent(out), optional               :: stat

! --- local variables
      integer(C_INT) :: err
      integer :: name_len
    
      name_len = len_trim(filename)
      call output_file_open_rdonly_f(filename,name_len,fid,err)
      if( present(stat) ) then
          call define_return_status(err,stat)
      endif    


end subroutine df90_output_file_open_rdonly

! ------------------------------------------------------------------------------

subroutine df90_output_file_open_rdwr(filename,fid,stat)

! --- calling arguments
      character(kind=C_CHAR,len=*), intent(in)  :: filename
      type(C_PTR),                     intent(inout) :: fid

      integer, intent(out), optional               :: stat

! --- local variables
      integer(C_INT) :: err
      integer :: name_len
    
      name_len = len_trim(filename)
      call output_file_open_rdwr_f(filename,name_len,fid,err)
      if( present(stat) ) then
          call define_return_status(err,stat)
      endif    

end subroutine df90_output_file_open_rdwr

! ------------------------------------------------------------------------------

subroutine df90_output_file_open_append(filename,fid,stat)

! --- calling arguments
      character(kind=C_CHAR,len=*), intent(in)    :: filename
      type(C_PTR),                  intent(inout) :: fid

      integer, intent(out), optional              :: stat

! --- local variables
      integer(C_INT) :: err
      integer :: name_len
    
      name_len = len_trim(filename)
      call output_file_open_append_f(filename,name_len,fid,err)
      if( present(stat) ) then
          call define_return_status(err,stat)
      endif    

end subroutine df90_output_file_open_append

! ------------------------------------------------------------------------------

subroutine df90_output_file_close(fid,stat)

! --- calling arguments
      type(C_PTR),                     intent(inout) :: fid

      integer, intent(out), optional               :: stat

! --- local variables
      integer(C_INT) :: err
    
      call output_file_close_f(fid,err)
      if( present(stat) ) then
          call define_return_status(err,stat)
      endif 

end subroutine df90_output_file_close

! ------------------------------------------------------------------------------

subroutine df90_attr_count(fid,num_found,stat)

!     
      type(C_PTR),         intent(in)  :: fid
      integer(C_INT),      intent(out) :: num_found

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: err

! --- code
      call danu_attr_count_f(fid,num_found,err)
      if( present(stat) ) then
          call define_return_status(err,stat)
      endif 

end subroutine df90_attr_count

! ------------------------------------------------------------------------------

subroutine df90_attr_names(loc,names,stat)

      type(C_PTR),                                intent(in)  :: loc
      character(kind=C_CHAR,len=*),               intent(out) :: names(:)

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT)                             :: num, length, ierr


! --- code
      length = len(names)
      num    = size(names)

      call danu_attr_names_f(loc,names,length,num,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       


end subroutine df90_attr_names

! ------------------------------------------------------------------------------

subroutine df90_attr_write_integer(id,attr_name,value,stat)

! --- Calling arguments

      type(C_PTR),       intent(in)  :: id
      character(len=*),  intent(in)  :: attr_name
      integer,           intent(in)  :: value

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr

! --- Code

      name_len = len_trim(attr_name)
      call danu_attr_write_int_f(id,attr_name,name_len,value,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_attr_write_integer

! ------------------------------------------------------------------------------

subroutine df90_attr_read_integer(id,attr_name,value,stat)

! --- Calling arguments

      type(C_PTR),       intent(in)   :: id
      character(len=*),  intent(in)   :: attr_name
      integer,           intent(out)  :: value

      integer, intent(out), optional  :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr

! --- Code

      name_len = len_trim(attr_name)
      call danu_attr_read_int_f(id,attr_name,name_len,value,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_attr_read_integer

! ------------------------------------------------------------------------------

subroutine df90_attr_write_real4(id,attr_name,value,stat)

! --- Calling arguments

      type(C_PTR),        intent(in)  :: id
      character(len=*),   intent(in)  :: attr_name
      real(kind=C_FLOAT), intent(in)  :: value

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr

! --- Code

      name_len = len_trim(attr_name)
      call danu_attr_write_real4_f(id,attr_name,name_len,value,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_attr_write_real4

! ------------------------------------------------------------------------------

subroutine df90_attr_read_real4(id,attr_name,value,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)   :: id
      character(len=*),    intent(in)   :: attr_name
      real(kind=C_FLOAT),  intent(out)  :: value

      integer, intent(out), optional    :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr

! --- Code

      name_len = len_trim(attr_name)
      call danu_attr_read_real4_f(id,attr_name,name_len,value,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_attr_read_real4

! ------------------------------------------------------------------------------

subroutine df90_attr_write_real8(id,attr_name,value,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: id
      character(len=*),    intent(in)  :: attr_name
      real(kind=C_DOUBLE), intent(in)  :: value

      integer, intent(out), optional   :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr

! --- Code

      name_len = len_trim(attr_name)
      call danu_attr_write_real8_f(id,attr_name,name_len,value,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_attr_write_real8

! ------------------------------------------------------------------------------

subroutine df90_attr_read_real8(id,attr_name,value,stat)

! --- Calling arguments

      type(C_PTR),          intent(in)   :: id
      character(len=*),     intent(in)   :: attr_name
      real(kind=C_DOUBLE),  intent(out)  :: value

      integer, intent(out), optional     :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr

! --- Code

      name_len = len_trim(attr_name)
      call danu_attr_read_real8_f(id,attr_name,name_len,value,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_attr_read_real8

! ------------------------------------------------------------------------------

subroutine df90_attr_write_character0(id,attr_name,char_data,stat)

! --- Calling arguments

      type(C_PTR),                  intent(in)  :: id
      character(len=*),             intent(in)  :: attr_name
      character(kind=C_CHAR,len=*), intent(in)  :: char_data

      integer, intent(out), optional   :: stat

! --- Local variables

      integer :: name_len, data_len
      integer(C_INT) :: ierr

! --- Code

      name_len = len_trim(attr_name)
      data_len = len(char_data)
      call danu_attr_write_char_f(id,attr_name,name_len,char_data,data_len,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_attr_write_character0

! ------------------------------------------------------------------------------

subroutine df90_attr_read_character0(id,attr_name,char_data,stat)

! --- Calling arguments

      type(C_PTR),                  intent(in)  :: id
      character(kind=C_CHAR,len=*), intent(in)  :: attr_name
      character(kind=C_CHAR,len=*), intent(out)  :: char_data

      integer, intent(out), optional   :: stat

! --- Local variables

      integer :: name_len, data_len
      integer(C_INT) :: ierr

! --- Code

      name_len = len_trim(attr_name)
      data_len = len(char_data)
      call danu_attr_read_char_f(id,attr_name,name_len,char_data,data_len,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       


      
end subroutine df90_attr_read_character0

! ------------------------------------------------------------------------------

subroutine df90_attr_write_character1(id,attr_name,char_data,stat)

! --- Calling arguments

      type(C_PTR),                  intent(in)  :: id
      character(len=*),             intent(in)  :: attr_name
      character(kind=C_CHAR),       intent(in)  :: char_data(:)

      integer, intent(out), optional   :: stat

! --- Local variables

      integer :: name_len, data_len
      integer(C_INT) :: ierr

! --- Code

      name_len = len_trim(attr_name)
      data_len = size(char_data)
      call danu_attr_write_char_f(id,attr_name,name_len,char_data,data_len,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_attr_write_character1

! ------------------------------------------------------------------------------

subroutine df90_group_create(loc,grp_name,gid,stat)

! --- calling arguments
      type(C_PTR),                  intent(in)  :: loc
      character(kind=C_CHAR,len=*), intent(in)  :: grp_name
      type(C_PTR),                  intent(out) :: gid

      integer, intent(out), optional               :: stat

! --- local variables
      integer :: name_len
      integer(C_INT) :: err
    
      name_len = len(grp_name)
      call danu_group_create_f(loc,grp_name,name_len,gid,err)
      if (present(stat)) then
          call define_return_status(err,stat)
      end if


end subroutine df90_group_create

! ------------------------------------------------------------------------------

subroutine df90_group_open(loc,grp_name,gid,stat)

! --- calling arguments
      type(C_PTR),                  intent(in)  :: loc
      character(kind=C_CHAR,len=*), intent(in)  :: grp_name
      type(C_PTR),                  intent(out) :: gid

      integer, intent(out), optional               :: stat

! --- local variables
      integer :: name_len
      integer(C_INT) :: err
    
      name_len = len(grp_name)
      call danu_group_open_f(loc,grp_name,name_len,gid,err)
      if (present(stat)) then
          call define_return_status(err,stat)
      end if


end subroutine df90_group_open

! ------------------------------------------------------------------------------

subroutine df90_group_close(gid,stat)

! --- calling arguments
      type(C_PTR),                  intent(inout) :: gid

      integer, intent(out), optional               :: stat

! --- local variables
      integer(C_INT) :: err
    
      call danu_group_close_f(gid,err)
      if (present(stat)) then
          call define_return_status(err,stat)
      end if


end subroutine df90_group_close

! ------------------------------------------------------------------------------


subroutine df90_mesh_count(fid,mesh_count,stat)

      type(C_PTR),    intent(in)  :: fid
      integer(C_INT), intent(out) :: mesh_count

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: err

! --- code

      call mesh_count_f(fid,mesh_count,err)
      if (present(stat) ) then
          call define_return_status(err,stat)
      endif

end subroutine df90_mesh_count

! ------------------------------------------------------------------------------

subroutine df90_mesh_list(fid,names,stat)

      type(C_PTR),                                intent(in)  :: fid
      character(kind=C_CHAR,len=*), dimension(:), intent(out) :: names

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: name_num, name_len, err

! --- code
      name_len = len(names)
      name_num = size(names)

      call mesh_list_f(fid,names,name_len,name_num,err)
      if (present(stat)) then
          call define_return_status(err,stat)
      endif       


end subroutine df90_mesh_list

! ------------------------------------------------------------------------------

subroutine df90_mesh_exists(fid,mesh_name,exists,stat)
      type(C_PTR),                  intent(in)  :: fid
      character(kind=C_CHAR,len=*), intent(in)  :: mesh_name
      logical,                      intent(out) :: exists

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: name_len, ierr, flag

! --- code

      name_len = len(mesh_name)
      flag = 0
      call mesh_exists_f(fid,mesh_name,name_len,flag,ierr)
      if ( flag .ne. 0 ) then
          exists = .true.
      else
          exists = .false.
      endif
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_mesh_exists

! ------------------------------------------------------------------------------
subroutine df90_mesh_open(fid,mesh_name,mid,stat)

! --- Arguments
      type(C_PTR),                  intent(in)  :: fid
      character(kind=C_CHAR,len=*), intent(in)  :: mesh_name
      type(C_PTR),                  intent(out) :: mid

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: name_len, ierr

! --- code

      name_len = len(mesh_name)
      call mesh_open_f(fid, mesh_name, name_len, mid, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 


end subroutine df90_mesh_open

! ------------------------------------------------------------------------------

subroutine df90_mesh_add_unstructured(fid,mesh_name,elem_order,mesh_dim,mid,stat)
      
      type(C_PTR),                  intent(in)  :: fid
      character(kind=C_CHAR,len=*), intent(in)  :: mesh_name
      integer(kind=C_INT),          intent(in)  :: elem_order
      integer(kind=C_INT),          intent(in)  :: mesh_dim
      type(C_PTR),                  intent(out) :: mid

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: name_len, ierr

! --- code

      name_len = len(mesh_name)
      call mesh_add_unstructured_f(fid, mesh_name, name_len,                   &
                                   elem_order, mesh_dim, mid, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_mesh_add_unstructured

! ------------------------------------------------------------------------------

subroutine df90_mesh_create_hex_unstruct(fid,                                  &
                                         mesh_name,                            &
                                         nelem,                                &
                                         nnodes,                               &
                                         x,y,z,                                &
                                         conn,                                 &
                                         mid,                                  &
                                         stat)
      
      type(C_PTR),                              intent(in)  :: fid
      character(kind=C_CHAR,len=*),             intent(in)  :: mesh_name
      integer(kind=C_INT),                      intent(in)  :: nelem
      integer(kind=C_INT),                      intent(in)  :: nnodes
      real(kind=C_DOUBLE), dimension(nnodes),               intent(in)  :: x, y, z
      integer(kind=C_INT), dimension(HEX_ELEM_ORDER,nelem), intent(in)  :: conn
      type(C_PTR),                                          intent(out) :: mid

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: name_len, ierr

! --- code

      name_len = len(mesh_name)
      call mesh_create_hex_unstruct_f(fid, mesh_name, name_len,                   &
                                      nnodes,x,y,z,                               &
                                      nelem,conn,mid,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_mesh_create_hex_unstruct

! ------------------------------------------------------------------------------

subroutine df90_mesh_write_coordinates_3d(mid,nnodes,x,y,z,stat)

      type(C_PTR),                       intent(in)  :: mid
      integer(kind=C_INT),               intent(in)  :: nnodes
      real(kind=C_DOUBLE), dimension(:), intent(in)  :: x
      real(kind=C_DOUBLE), dimension(:), intent(in)  :: y
      real(kind=C_DOUBLE), dimension(:), intent(in)  :: z

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: ierr

! --- code

      call mesh_write_coordinates_f(mid,nnodes,x,y,z,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 


end subroutine df90_mesh_write_coordinates_3d

! ------------------------------------------------------------------------------

subroutine df90_mesh_write_coordinates_2d(mid,nnodes,x,y,stat)

      type(C_PTR),                       intent(in)  :: mid
      integer(kind=C_INT),               intent(in)  :: nnodes
      real(kind=C_DOUBLE), dimension(:), intent(in)  :: x
      real(kind=C_DOUBLE), dimension(:), intent(in)  :: y

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: ierr

! --- code

      call mesh_write_coordinates_2d_f(mid,nnodes,x,y,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 


end subroutine df90_mesh_write_coordinates_2d

! ------------------------------------------------------------------------------

subroutine df90_mesh_write_coordinates_1d(mid,nnodes,x,stat)

      type(C_PTR),                       intent(in)  :: mid
      integer(kind=C_INT),               intent(in)  :: nnodes
      real(kind=C_DOUBLE), dimension(:), intent(in)  :: x

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: ierr

! --- code

      call mesh_write_coordinates_1d_f(mid,nnodes,x,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 


end subroutine df90_mesh_write_coordinates_1d

! ------------------------------------------------------------------------------

subroutine df90_mesh_read_coordinates_3d(mid,x,y,z,stat)

      type(C_PTR),                       intent(in)  :: mid
      real(kind=C_DOUBLE), dimension(:), intent(out) :: x
      real(kind=C_DOUBLE), dimension(:), intent(out) :: y
      real(kind=C_DOUBLE), dimension(:), intent(out) :: z

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: ierr

! --- code

      call mesh_read_coordinates_f(mid,x,y,z,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 


end subroutine df90_mesh_read_coordinates_3d

! ------------------------------------------------------------------------------

subroutine df90_mesh_read_coordinates_idx(mid,idx,buf,stat)

      type(C_PTR),                       intent(in)  :: mid
      integer(C_INT),                    intent(in)  :: idx
      real(kind=C_DOUBLE), dimension(:), intent(out) :: buf

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: ierr

! --- code

      call mesh_read_coordinates_byindex_f(mid,idx,buf,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_mesh_read_coordinates_idx

! ------------------------------------------------------------------------------

subroutine df90_mesh_read_coordinates_2d(mid,x,y,stat)

      type(C_PTR),                       intent(in)  :: mid
      real(kind=C_DOUBLE), dimension(:), intent(out) :: x
      real(kind=C_DOUBLE), dimension(:), intent(out) :: y

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: ierr

! --- code

      call mesh_read_coordinates_2d_f(mid,x,y,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 


end subroutine df90_mesh_read_coordinates_2d

! ------------------------------------------------------------------------------

subroutine df90_mesh_read_coordinates_1d(mid,x,stat)

      type(C_PTR),                       intent(in)  :: mid
      real(kind=C_DOUBLE), dimension(:), intent(out) :: x

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: ierr

! --- code

      call mesh_read_coordinates_1d_f(mid,x,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 


end subroutine df90_mesh_read_coordinates_1d

! ------------------------------------------------------------------------------

subroutine df90_mesh_write_connectivity(mid,num,idata,stat)

! --- Calling Arguments
 
      type(C_PTR),                         intent(in) :: mid
      integer(kind=C_INT),                 intent(in) :: num
      integer(kind=C_INT), dimension(:,:), intent(in) :: idata              

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

      call mesh_write_connectivity_f(mid,num,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_mesh_write_connectivity

! ------------------------------------------------------------------------------

subroutine df90_mesh_read_connectivity(mid,idata,stat)

! --- Calling Arguments
 
      type(C_PTR),                         intent(in)  :: mid
      integer(kind=C_INT), dimension(:,:), intent(out) :: idata              

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

      call mesh_read_connectivity_f(mid,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_mesh_read_connectivity

! ------------------------------------------------------------------------------

subroutine df90_mesh_get_type(mid,mesh_type,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: mid
      integer(kind=C_INT),               intent(out) :: mesh_type              

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

      call mesh_get_type_f(mid,mesh_type,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_mesh_get_type

! ------------------------------------------------------------------------------

subroutine df90_mesh_get_elementtype(mid,elem_type,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: mid
      integer(kind=C_INT),               intent(out) :: elem_type              

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

      call mesh_get_elementtype_f(mid,elem_type,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_mesh_get_elementtype

! ------------------------------------------------------------------------------

subroutine df90_mesh_get_dimension(mid,mesh_dim,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: mid
      integer(kind=C_INT),               intent(out) :: mesh_dim              

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

      call mesh_get_dimension_f(mid,mesh_dim,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_mesh_get_dimension

! ------------------------------------------------------------------------------

subroutine df90_mesh_get_nnodes(mid,nnodes,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: mid
      integer(kind=C_INT),               intent(out) :: nnodes              

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

      call mesh_get_nnodes_f(mid,nnodes,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_mesh_get_nnodes

! ------------------------------------------------------------------------------

subroutine df90_mesh_get_nelem(mid,nelem,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: mid
      integer(kind=C_INT),               intent(out) :: nelem              

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

      call mesh_get_nelem_f(mid,nelem,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_mesh_get_nelem

! ------------------------------------------------------------------------------

subroutine df90_mesh_get_elem_order(mid,elem_order,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: mid
      integer(kind=C_INT),               intent(out) :: elem_order              

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

      call mesh_get_elem_order_f(mid,elem_order,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_mesh_get_elem_order

! ------------------------------------------------------------------------------

subroutine df90_simulation_add(fid,sim_name,sid,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: fid
      character(kind=C_CHAR,len=*),      intent(in)  :: sim_name
      type(C_PTR),                       intent(out) :: sid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, ierr

! --- code

      name_len = len(sim_name)
      call simulation_add_f(fid, sim_name, name_len, sid, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_simulation_add

! ------------------------------------------------------------------------------

subroutine df90_simulation_open(fid,sim_name,sid,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: fid
      character(kind=C_CHAR,len=*),      intent(in)  :: sim_name
      type(C_PTR),                       intent(out) :: sid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, ierr

! --- code

      name_len = len(sim_name)
      call simulation_open_f(fid, sim_name, name_len, sid, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_simulation_open

! ------------------------------------------------------------------------------

subroutine df90_simulation_count(fid,cnt,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: fid
      integer(kind=C_INT),               intent(out) :: cnt

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

! --- code

      call simulation_count_f(fid, cnt, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_simulation_count

! ------------------------------------------------------------------------------

subroutine df90_simulation_link_mesh(fid,sid,mesh_name,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: fid
      type(C_PTR),                       intent(in)  :: sid
      character(kind=C_CHAR,len=*),      intent(in)  :: mesh_name

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, ierr

! --- code

      name_len = len(mesh_name)
      call simulation_link_mesh_f(fid, sid, mesh_name, name_len, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_simulation_link_mesh

! ------------------------------------------------------------------------------

subroutine df90_simulation_open_mesh_link(sid,mid,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      type(C_PTR),                       intent(out) :: mid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

! --- code

      call simulation_open_mesh_link_f(sid, mid, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_simulation_open_mesh_link

! ------------------------------------------------------------------------------

logical function df90_simulation_mesh_link_exists(sid)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid

! --- Local variables
      integer(C_INT) :: flag, ierr, stat

! --- code

      call simulation_mesh_link_exists_f(sid, flag, ierr)
      call define_return_status(ierr,stat)

      if ( stat .eq. DANU_SUCCESS ) then
        if ( flag .gt. 0 ) then
          df90_simulation_mesh_link_exists = .true.
        else
          df90_simulation_mesh_link_exists = .false.
        endif
      else    
        df90_simulation_mesh_link_exists = .false.
      endif

end function df90_simulation_mesh_link_exists

! ------------------------------------------------------------------------------

subroutine df90_simulation_list(fid,names,stat)

      type(C_PTR),                                intent(in)  :: fid
      character(kind=C_CHAR,len=*), dimension(:), intent(out) :: names

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: name_num, name_len, ierr

! --- code
      name_len = len(names)
      name_num = size(names)

      call simulation_list_f(fid,names,name_len,name_num,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       


end subroutine df90_simulation_list

! ------------------------------------------------------------------------------

subroutine df90_data_open_dataset(sid,data_name,hid,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      character(len=*,kind=C_CHAR),      intent(in)  :: data_name
      type(C_PTR),                       intent(out) :: hid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr, name_len

! --- code
      name_len=len_trim(data_name)
      call data_open_dataset_f(sid,data_name,name_len,hid,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif  

end subroutine df90_data_open_dataset

! ------------------------------------------------------------------------------

subroutine df90_data_count(sid,cnt,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      integer(kind=C_INT),               intent(out) :: cnt

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

! --- code

      call data_count_f(sid, cnt, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_data_count

! ------------------------------------------------------------------------------

subroutine df90_data_type(sid,data_name,typecode,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      character(len=*),                  intent(in)  :: data_name
      integer(kind=C_INT),               intent(out) :: typecode

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr, name_len

! --- code
      name_len=len_trim(data_name)
      call data_type_f(sid,data_name,name_len,typecode,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_data_type

! ------------------------------------------------------------------------------

subroutine df90_data_rank(sid,data_name,rank,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      character(len=*),                  intent(in)  :: data_name
      integer(kind=C_INT),               intent(out) :: rank

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr, name_len

! --- code
      name_len=len_trim(data_name)
      call data_rank_f(sid,data_name,name_len,rank,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_data_rank

! ------------------------------------------------------------------------------

subroutine df90_data_dimensions(sid,data_name,dims,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      character(len=*),                  intent(in)  :: data_name
      integer(kind=C_INT), dimension(:), intent(out) :: dims

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr, name_len, r

! --- code
      r=size(dims,1)
      name_len=len_trim(data_name)
      call data_dimensions_f(sid,data_name,name_len,r,dims,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_data_dimensions

! ------------------------------------------------------------------------------

subroutine df90_data_list(sid,names,stat)

      type(C_PTR),                                intent(in)  :: sid
      character(kind=C_CHAR,len=*), dimension(:), intent(out) :: names

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: name_num, name_len, ierr

! --- code
      name_len = len(names)
      name_num = size(names)

      call data_list_f(sid,names,name_len,name_num,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_list

! ------------------------------------------------------------------------------

subroutine df90_data_write_byte_rank0 (sid, data_name, data, stat)

  type(C_PTR), intent(in) :: sid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(in) :: data
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr

  name_len = len_trim(data_name)
  call data_write_byte0_f (sid, data_name, name_len, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_data_write_byte_rank0

subroutine df90_data_write_byte_rank1 (sid, data_name, data, stat)

  type(C_PTR), intent(in) :: sid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(in) :: data(:)
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr, dims(1)
  integer(C_INT), parameter :: num_dim = 1

  name_len = len_trim(data_name)
  dims = shape(data)
  call data_write_byte_f (sid, data_name, name_len, num_dim, dims, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_data_write_byte_rank1

subroutine df90_data_write_byte_rank2 (sid, data_name, data, stat)

  type(C_PTR), intent(in) :: sid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(in) :: data(:,:)
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr, dims(2)
  integer(C_INT), parameter :: num_dim = 2

  name_len = len_trim(data_name)
  dims = shape(data)
  call data_write_byte_f (sid, data_name, name_len, num_dim, dims, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_data_write_byte_rank2

subroutine df90_data_write_byte_rank3 (sid, data_name, data, stat)

  type(C_PTR), intent(in) :: sid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(in) :: data(:,:,:)
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr, dims(3)
  integer(C_INT), parameter :: num_dim = 3

  name_len = len_trim(data_name)
  dims = shape(data)
  call data_write_byte_f (sid, data_name, name_len, num_dim, dims, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_data_write_byte_rank3

! ------------------------------------------------------------------------------

subroutine df90_data_write_integer_rank0(sid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), intent(in)  :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr


! --- Code

      name_len = len_trim(data_name)
      call data_write_int0_f(sid,data_name,name_len,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_write_integer_rank0

! ------------------------------------------------------------------------------

subroutine df90_data_write_integer_rank1(sid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), dimension(:), intent(in) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:1) :: dims

      integer(C_INT), parameter  :: num_dim = 1

! --- Code

      name_len = len_trim(data_name)
      dims = shape(idata)
      call data_write_int_f(sid,data_name,name_len,num_dim,dims,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_write_integer_rank1

! ------------------------------------------------------------------------------

subroutine df90_data_write_integer_rank2(sid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), dimension(:,:), intent(in) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:2) :: dims

      integer(C_INT), parameter  :: num_dim = 2

! --- Code

      name_len = len_trim(data_name)
      dims = shape(idata)
      call data_write_int_f(sid,data_name,name_len,num_dim,dims,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_write_integer_rank2

! ------------------------------------------------------------------------------

subroutine df90_data_write_integer_rank3(sid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), dimension(:,:,:), intent(in) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:3) :: dims

      integer(C_INT), parameter  :: num_dim = 3

! --- Code

      name_len = len_trim(data_name)
      dims = shape(idata)
      call data_write_int_f(sid,data_name,name_len,num_dim,dims,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_write_integer_rank3

! ------------------------------------------------------------------------------

subroutine df90_data_write_real4_rank0(sid,data_name,r4data,stat)

! --- Calling arguments

      type(C_PTR),        intent(in)  :: sid
      character(len=*),   intent(in)  :: data_name
      real(kind=C_FLOAT), intent(in)  :: r4data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr


! --- Code

      name_len = len_trim(data_name)
      call data_write_float0_f(sid,data_name,name_len,r4data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_write_real4_rank0

! ------------------------------------------------------------------------------


subroutine df90_data_write_real4_rank1(sid,data_name,r4data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_FLOAT), dimension(:), intent(in) :: r4data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:1) :: dims

      integer(C_INT), parameter  :: num_dim = 1

! --- Code

      name_len = len_trim(data_name)
      dims = shape(r4data)
      call data_write_float_f(sid,data_name,name_len,num_dim,dims,r4data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_write_real4_rank1

! ------------------------------------------------------------------------------

subroutine df90_data_write_real4_rank2(sid,data_name,r4data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_FLOAT), dimension(:,:), intent(in) :: r4data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:2) :: dims

      integer(C_INT), parameter  :: num_dim = 2

! --- Code

      name_len = len_trim(data_name)
      dims = shape(r4data)
      call data_write_float_f(sid,data_name,name_len,num_dim,dims,r4data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_write_real4_rank2

! ------------------------------------------------------------------------------

subroutine df90_data_write_real4_rank3(sid,data_name,r4data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_FLOAT), dimension(:,:,:), intent(in) :: r4data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:3) :: dims

      integer(C_INT), parameter  :: num_dim = 3

! --- Code

      name_len = len_trim(data_name)
      dims = shape(r4data)
      call data_write_float_f(sid,data_name,name_len,num_dim,dims,r4data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_write_real4_rank3

! ------------------------------------------------------------------------------

subroutine df90_data_write_real8_rank0(sid,data_name,r8data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_DOUBLE), intent(in)  :: r8data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr


! --- Code

      name_len = len_trim(data_name)
      call data_write_double0_f(sid,data_name,name_len,r8data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_write_real8_rank0

! ------------------------------------------------------------------------------

subroutine df90_data_write_real8_rank1(sid,data_name,r8data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_DOUBLE), dimension(:), intent(in) :: r8data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:1) :: dims

      integer(C_INT), parameter  :: num_dim = 1

! --- Code

      name_len = len_trim(data_name)
      dims = shape(r8data)
      call data_write_double_f(sid,data_name,name_len,num_dim,dims,r8data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_write_real8_rank1

! ------------------------------------------------------------------------------

subroutine df90_data_write_real8_rank2(sid,data_name,r8data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_DOUBLE), dimension(:,:), intent(in) :: r8data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:2) :: dims

      integer(C_INT), parameter  :: num_dim = 2

! --- Code

      name_len = len_trim(data_name)
      dims = shape(r8data)
      call data_write_double_f(sid,data_name,name_len,num_dim,dims,r8data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_write_real8_rank2

! ------------------------------------------------------------------------------

subroutine df90_data_write_real8_rank3(sid,data_name,r8data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_DOUBLE), dimension(:,:,:), intent(in) :: r8data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:3) :: dims

      integer(C_INT), parameter  :: num_dim = 3

! --- Code

      name_len = len_trim(data_name)
      dims = shape(r8data)
      call data_write_double_f(sid,data_name,name_len,num_dim,dims,r8data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_write_real8_rank3

! ------------------------------------------------------------------------------

subroutine df90_data_read_byte_rank0 (sid, data_name, data, stat)

  type(C_PTR), intent(in) :: sid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(out) :: data            
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr

  name_len = len_trim(data_name)
  call data_read_byte0_f (sid, data_name, name_len, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_data_read_byte_rank0

subroutine df90_data_read_byte_rank1 (sid, data_name, data, stat)

  type(C_PTR), intent(in) :: sid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(out) :: data(:)
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr, dims(1)
  integer(C_INT), parameter :: num_dim = 1

  name_len = len_trim(data_name)
  dims = shape(data)
  call data_read_byte_f (sid, data_name, name_len, num_dim, dims, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_data_read_byte_rank1

subroutine df90_data_read_byte_rank2 (sid, data_name, data, stat)

  type(C_PTR), intent(in) :: sid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(out) :: data(:,:)
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr, dims(2)
  integer(C_INT), parameter :: num_dim = 2

  name_len = len_trim(data_name)
  dims = shape(data)
  call data_read_byte_f (sid, data_name, name_len, num_dim, dims, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_data_read_byte_rank2

subroutine df90_data_read_byte_rank3 (sid, data_name, data, stat)

  type(C_PTR), intent(in) :: sid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(out) :: data(:,:,:)
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr, dims(3)
  integer(C_INT), parameter :: num_dim = 3

  name_len = len_trim(data_name)
  dims = shape(data)
  call data_read_byte_f (sid, data_name, name_len, num_dim, dims, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_data_read_byte_rank3

! ------------------------------------------------------------------------------

subroutine df90_data_read_integer_rank0(sid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), intent(out) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr

! --- Code

      name_len = len_trim(data_name)
      call data_read_int0_f(sid,data_name,name_len,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_read_integer_rank0

! ------------------------------------------------------------------------------

subroutine df90_data_read_integer_rank1(sid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), dimension(:), intent(out) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:1) :: dims

      integer(C_INT), parameter  :: num_dim = 1

! --- Code

      name_len = len_trim(data_name)
      dims = shape(idata)
      call data_read_int_f(sid,data_name,name_len,num_dim,dims,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_read_integer_rank1

! ------------------------------------------------------------------------------

subroutine df90_data_read_integer_rank2(sid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), dimension(:,:), intent(out) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:2) :: dims

      integer(C_INT), parameter  :: num_dim = 2

! --- Code

      name_len = len_trim(data_name)
      dims = shape(idata)
      call data_read_int_f(sid,data_name,name_len,num_dim,dims,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_read_integer_rank2

! ------------------------------------------------------------------------------

subroutine df90_data_read_integer_rank3(sid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), dimension(:,:,:), intent(out) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:3) :: dims

      integer(C_INT), parameter  :: num_dim = 3

! --- Code

      name_len = len_trim(data_name)
      dims = shape(idata)
      call data_read_int_f(sid,data_name,name_len,num_dim,dims,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_read_integer_rank3

! ------------------------------------------------------------------------------

subroutine df90_data_read_real4_rank0(sid,data_name,r4data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_FLOAT),  intent(out) :: r4data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr

! --- Code

      name_len = len_trim(data_name)
      call data_read_float0_f(sid,data_name,name_len,r4data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_read_real4_rank0

! ------------------------------------------------------------------------------


subroutine df90_data_read_real4_rank1(sid,data_name,r4data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_FLOAT), dimension(:), intent(out) :: r4data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:1) :: dims

      integer(C_INT), parameter  :: num_dim = 1

! --- Code

      name_len = len_trim(data_name)
      dims = shape(r4data)
      call data_read_float_f(sid,data_name,name_len,num_dim,dims,r4data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_read_real4_rank1

! ------------------------------------------------------------------------------

subroutine df90_data_read_real4_rank2(sid,data_name,r4data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_FLOAT), dimension(:,:), intent(out) :: r4data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:2) :: dims

      integer(C_INT), parameter  :: num_dim = 2

! --- Code

      name_len = len_trim(data_name)
      dims = shape(r4data)
      call data_read_float_f(sid,data_name,name_len,num_dim,dims,r4data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_read_real4_rank2

! ------------------------------------------------------------------------------

subroutine df90_data_read_real4_rank3(sid,data_name,r4data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_FLOAT), dimension(:,:,:), intent(out) :: r4data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:3) :: dims

      integer(C_INT), parameter  :: num_dim = 3

! --- Code

      name_len = len_trim(data_name)
      dims = shape(r4data)
      call data_read_float_f(sid,data_name,name_len,num_dim,dims,r4data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_read_real4_rank3

! ------------------------------------------------------------------------------

subroutine df90_data_read_real8_rank0(sid,data_name,r8data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_DOUBLE), intent(out) :: r8data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr

! --- Code

      name_len = len_trim(data_name)
      call data_read_double0_f(sid,data_name,name_len,r8data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_read_real8_rank0

! ------------------------------------------------------------------------------

subroutine df90_data_read_real8_rank1(sid,data_name,r8data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_DOUBLE), dimension(:), intent(out) :: r8data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:1) :: dims

      integer(C_INT), parameter  :: num_dim = 1

! --- Code

      name_len = len_trim(data_name)
      dims = shape(r8data)
      call data_read_double_f(sid,data_name,name_len,num_dim,dims,r8data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_read_real8_rank1

! ------------------------------------------------------------------------------

subroutine df90_data_read_real8_rank2(sid,data_name,r8data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_DOUBLE), dimension(:,:), intent(out) :: r8data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:2) :: dims

      integer(C_INT), parameter  :: num_dim = 2

! --- Code

      name_len = len_trim(data_name)
      dims = shape(r8data)
      call data_read_double_f(sid,data_name,name_len,num_dim,dims,r8data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_read_real8_rank2

! ------------------------------------------------------------------------------

subroutine df90_data_read_real8_rank3(sid,data_name,r8data,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: sid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_DOUBLE), dimension(:,:,:), intent(out) :: r8data            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:3) :: dims

      integer(C_INT), parameter  :: num_dim = 3

! --- Code

      name_len = len_trim(data_name)
      dims = shape(r8data)
      call data_read_double_f(sid,data_name,name_len,num_dim,dims,r8data,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_data_read_real8_rank3

! ------------------------------------------------------------------------------

subroutine df90_sequence_count(sid,cnt,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      integer(kind=C_INT),               intent(out) :: cnt

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

! --- code

      call sequence_count_f(sid, cnt, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_sequence_count

! ------------------------------------------------------------------------------

subroutine df90_sequence_list(sid,names,stat)

      type(C_PTR),                                intent(in)  :: sid
      character(kind=C_CHAR,len=*), dimension(:), intent(out) :: names

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: name_num, name_len, ierr

! --- code
      name_len = len(names)
      name_num = size(names)

      call sequence_list_f(sid,names,name_len,name_num,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sequence_list

! ------------------------------------------------------------------------------

subroutine df90_sequence_get_next_id(sid,cyc,time,nid,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      integer(kind=C_INT),               intent(in)  :: cyc 
      real(kind=C_DOUBLE),               intent(in)  :: time
      type(C_PTR),                       intent(out) :: nid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

! --- code

      call sequence_get_next_id_f(sid,cyc,time, nid, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_sequence_get_next_id

! ------------------------------------------------------------------------------

subroutine df90_sequence_get_id_c(sid,sname,nid,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      character(kind=C_CHAR,len=*),      intent(in)  :: sname
      type(C_PTR),                       intent(out) :: nid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len
      integer(C_INT) :: ierr

! --- code

      name_len = len(sname)
      call sequence_get_id_f(sid,sname,name_len,nid,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_sequence_get_id_c

! ------------------------------------------------------------------------------

subroutine df90_sequence_get_id_i(sid,num,nid,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      integer(kind=C_INT),               intent(in)  :: num
      type(C_PTR),                       intent(out) :: nid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

! --- code

      call sequence_get_id_bynum_f(sid,num,nid,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_sequence_get_id_i

! ------------------------------------------------------------------------------

subroutine df90_simulation_open_data(nid,dataname,id,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: nid
      character(len=*,kind=C_CHAR),      intent(in)  :: dataname
      type(C_PTR),                       intent(out) :: id

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr, flen

! --- code

      flen = len_trim(dataname)
      call simulation_open_data_f(nid,dataname,flen,id,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_simulation_open_data

! ------------------------------------------------------------------------------

subroutine df90_simulation_data_count(nid,cnt,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: nid
      integer(kind=C_INT),               intent(out) :: cnt

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

! --- code

      call simulation_data_count_f(nid, cnt, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_simulation_data_count

! ------------------------------------------------------------------------------

subroutine df90_simulation_data_rank(nid,data_name,rank,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: nid
      character(len=*),                  intent(in)  :: data_name
      integer(kind=C_INT),               intent(out) :: rank

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, ierr

! --- code
      name_len=len(data_name)
      call simulation_data_rank_f(nid, data_name,name_len,rank, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_simulation_data_rank

! ------------------------------------------------------------------------------

subroutine df90_simulation_data_type(nid,data_name,code,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: nid
      character(len=*),                  intent(in)  :: data_name
      integer(kind=C_INT),               intent(out) :: code

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, ierr

! --- code
      name_len=len(data_name)
      call simulation_data_type_f(nid, data_name,name_len,code,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_simulation_data_type

! ------------------------------------------------------------------------------

subroutine df90_simulation_data_dimensions(nid,data_name,dimensions,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: nid
      character(len=*),                  intent(in)  :: data_name
      integer(kind=C_INT), dimension(:), intent(out) :: dimensions

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, ndims, ierr

! --- code
      name_len=len_trim(data_name)
      ndims=size(dimensions,1)
      call simulation_data_dimensions_f(nid, data_name,name_len,ndims,dimensions,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_simulation_data_dimensions

! ------------------------------------------------------------------------------

subroutine df90_simulation_data_list(nid,names,stat)

      type(C_PTR),                                intent(in)  :: nid
      character(kind=C_CHAR,len=*), dimension(:), intent(out) :: names

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: name_num, name_len, ierr

! --- code
      name_len = len(names)
      name_num = size(names)

      call simulation_data_list_f(nid,names,name_len,name_num,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_simulation_data_list

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_byte_rank0 (nid, data_name, data, stat)

  type(C_PTR), intent(in) :: nid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(in) :: data            
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr

  name_len = len_trim(data_name)
  call simulation_data_write_byte0_f (nid, data_name, name_len, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_sim_data_write_byte_rank0

subroutine df90_sim_data_write_byte_rank1 (nid, data_name, data, stat)

  type(C_PTR), intent(in) :: nid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(in) :: data(:)          
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr, dims(1)
  integer(C_INT), parameter :: num_dim = 1

  name_len = len_trim(data_name)
  dims = shape(data)
  call simulation_data_write_byte_f (nid, data_name, name_len, num_dim, dims, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_sim_data_write_byte_rank1

subroutine df90_sim_data_write_byte_rank2 (nid, data_name, data, stat)

  type(C_PTR), intent(in) :: nid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(in) :: data(:,:)          
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr, dims(2)
  integer(C_INT), parameter :: num_dim = 2

  name_len = len_trim(data_name)
  dims = shape(data)
  call simulation_data_write_byte_f (nid, data_name, name_len, num_dim, dims, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_sim_data_write_byte_rank2

subroutine df90_sim_data_write_byte_rank3 (nid, data_name, data, stat)

  type(C_PTR), intent(in) :: nid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(in) :: data(:,:,:)        
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr, dims(3)
  integer(C_INT), parameter :: num_dim = 3

  name_len = len_trim(data_name)
  dims = shape(data)
  call simulation_data_write_byte_f (nid, data_name, name_len, num_dim, dims, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_sim_data_write_byte_rank3

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_integer_rank0(nid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: nid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), intent(in)  :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr


! --- Code

      name_len = len_trim(data_name)
      call simulation_data_write_int0_f(nid,data_name,name_len,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_write_integer_rank0

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_integer_rank1(nid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: nid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), dimension(:), intent(in) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:1) :: dims

      integer(C_INT), parameter  :: num_dim = 1

! --- Code

      name_len = len_trim(data_name)
      dims = shape(idata)
      call simulation_data_write_int_f(nid,data_name,name_len,num_dim,dims,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_write_integer_rank1

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_integer_rank2(nid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: nid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), dimension(:,:), intent(in) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:2) :: dims

      integer(C_INT), parameter  :: num_dim = 2

! --- Code

      name_len = len_trim(data_name)
      dims = shape(idata)
      call simulation_data_write_int_f(nid,data_name,name_len,num_dim,dims,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_write_integer_rank2

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_integer_rank3(nid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: nid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), dimension(:,:,:), intent(in) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:3) :: dims

      integer(C_INT), parameter  :: num_dim = 3

! --- Code

      name_len = len_trim(data_name)
      dims = shape(idata)
      call simulation_data_write_int_f(nid,data_name,name_len,num_dim,dims,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_write_integer_rank3

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_real4_rank0(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),        intent(in)  :: nid
      character(len=*),   intent(in)  :: data_name
      real(kind=C_FLOAT), intent(in)  :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr


! --- Code

      name_len = len_trim(data_name)
      call simulation_data_write_float0_f(nid,data_name,name_len,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_write_real4_rank0

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_real4_rank1(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),                       intent(in)  :: nid
      character(len=*),                  intent(in)  :: data_name
      real(kind=C_FLOAT),  dimension(:), intent(in) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:1) :: dims

      integer(C_INT), parameter  :: num_dim = 1

! --- Code

      name_len = len_trim(data_name)
      dims = shape(rdata)
      call simulation_data_write_float_f(nid,data_name,name_len,num_dim,dims,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_write_real4_rank1

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_real4_rank2(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),                         intent(in)  :: nid
      character(len=*),                    intent(in)  :: data_name
      real(kind=C_FLOAT),  dimension(:,:), intent(in) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:2) :: dims

      integer(C_INT), parameter  :: num_dim = 2

! --- Code

      name_len = len_trim(data_name)
      dims = shape(rdata)
      call simulation_data_write_float_f(nid,data_name,name_len,num_dim,dims,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_write_real4_rank2

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_real4_rank3(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),                           intent(in) :: nid
      character(len=*),                      intent(in) :: data_name
      real(kind=C_FLOAT),  dimension(:,:,:), intent(in) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:3) :: dims

      integer(C_INT), parameter  :: num_dim = 3

! --- Code

      name_len = len_trim(data_name)
      dims = shape(rdata)
      call simulation_data_write_float_f(nid,data_name,name_len,num_dim,dims,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_write_real4_rank3

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_real8_rank0(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: nid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_DOUBLE), intent(in)  :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr


! --- Code

      name_len = len_trim(data_name)
      call simulation_data_write_double0_f(nid,data_name,name_len,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_write_real8_rank0

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_real8_rank1(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),                        intent(in)  :: nid
      character(len=*),                   intent(in)  :: data_name
      real(kind=C_DOUBLE),  dimension(:), intent(in) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:1) :: dims

      integer(C_INT), parameter  :: num_dim = 1

! --- Code

      name_len = len_trim(data_name)
      dims = shape(rdata)
      call simulation_data_write_double_f(nid,data_name,name_len,num_dim,dims,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_write_real8_rank1

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_real8_rank2(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),                         intent(in)  :: nid
      character(len=*),                    intent(in)  :: data_name
      real(kind=C_DOUBLE), dimension(:,:), intent(in) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:2) :: dims

      integer(C_INT), parameter  :: num_dim = 2

! --- Code

      name_len = len_trim(data_name)
      dims = shape(rdata)
      call simulation_data_write_double_f(nid,data_name,name_len,num_dim,dims,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_write_real8_rank2

! ------------------------------------------------------------------------------

subroutine df90_sim_data_write_real8_rank3(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),                           intent(in) :: nid
      character(len=*),                      intent(in) :: data_name
      real(kind=C_DOUBLE), dimension(:,:,:), intent(in) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:3) :: dims

      integer(C_INT), parameter  :: num_dim = 3

! --- Code

      name_len = len_trim(data_name)
      dims = shape(rdata)
      call simulation_data_write_double_f(nid,data_name,name_len,num_dim,dims,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_write_real8_rank3

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_byte_rank0 (nid,data_name,data,stat)

  type(C_PTR), intent(in) :: nid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(out) :: data            
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr

  name_len = len_trim(data_name)
  call simulation_data_read_byte0_f (nid,data_name,name_len,data,ierr)
  if (present(stat)) call define_return_status (ierr,stat)

end subroutine df90_sim_data_read_byte_rank0

subroutine df90_sim_data_read_byte_rank1 (nid, data_name, data, stat)

  type(C_PTR), intent(in) :: nid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(out) :: data(:)
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr, dims(1)
  integer(C_INT), parameter :: num_dim = 1

  name_len = len_trim(data_name)
  dims = shape(data)
  call simulation_data_read_byte_f (nid, data_name, name_len, num_dim, dims, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_sim_data_read_byte_rank1

subroutine df90_sim_data_read_byte_rank2 (nid, data_name, data, stat)

  type(C_PTR), intent(in) :: nid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(out) :: data(:,:)
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr, dims(2)
  integer(C_INT), parameter :: num_dim = 2

  name_len = len_trim(data_name)
  dims = shape(data)
  call simulation_data_read_byte_f (nid, data_name, name_len, num_dim, dims, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_sim_data_read_byte_rank2

subroutine df90_sim_data_read_byte_rank3 (nid, data_name, data, stat)

  type(C_PTR), intent(in) :: nid
  character(*), intent(in) :: data_name
  integer(C_INT8_T), intent(out) :: data(:,:,:)
  integer, intent(out), optional :: stat

  integer :: name_len
  integer(C_INT) :: ierr, dims(3)
  integer(C_INT), parameter :: num_dim = 3

  name_len = len_trim(data_name)
  dims = shape(data)
  call simulation_data_read_byte_f (nid, data_name, name_len, num_dim, dims, data, ierr)
  if (present(stat)) call define_return_status (ierr, stat)

end subroutine df90_sim_data_read_byte_rank3

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_integer_rank0(nid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: nid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), intent(out)  :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr


! --- Code

      name_len = len_trim(data_name)
      call simulation_data_read_int0_f(nid,data_name,name_len,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_read_integer_rank0

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_integer_rank1(nid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: nid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), dimension(:), intent(out) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:1) :: dims

      integer(C_INT), parameter  :: num_dim = 1

! --- Code

      name_len = len_trim(data_name)
      dims = shape(idata)
      call simulation_data_read_int_f(nid,data_name,name_len,num_dim,dims,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_read_integer_rank1

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_integer_rank2(nid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: nid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), dimension(:,:), intent(out) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:2) :: dims

      integer(C_INT), parameter  :: num_dim = 2

! --- Code

      name_len = len_trim(data_name)
      dims = shape(idata)
      call simulation_data_read_int_f(nid,data_name,name_len,num_dim,dims,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_read_integer_rank2

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_integer_rank3(nid,data_name,idata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: nid
      character(len=*),    intent(in)  :: data_name
      integer(kind=C_INT), dimension(:,:,:), intent(out) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:3) :: dims

      integer(C_INT), parameter  :: num_dim = 3

! --- Code

      name_len = len_trim(data_name)
      dims = shape(idata)
      call simulation_data_read_int_f(nid,data_name,name_len,num_dim,dims,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_read_integer_rank3

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_real4_rank0(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),        intent(in)  :: nid
      character(len=*),   intent(in)  :: data_name
      real(kind=C_FLOAT), intent(out) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr


! --- Code

      name_len = len_trim(data_name)
      call simulation_data_read_float0_f(nid,data_name,name_len,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_read_real4_rank0

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_real4_rank1(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),                      intent(in)  :: nid
      character(len=*),                 intent(in)  :: data_name
      real(kind=C_FLOAT), dimension(:), intent(out) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:1) :: dims

      integer(C_INT), parameter  :: num_dim = 1

! --- Code

      name_len = len_trim(data_name)
      dims = shape(rdata)
      call simulation_data_read_float_f(nid,data_name,name_len,num_dim,dims,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_read_real4_rank1

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_real4_rank2(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),                        intent(in)  :: nid
      character(len=*),                   intent(in)  :: data_name
      real(kind=C_FLOAT), dimension(:,:), intent(out) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:2) :: dims

      integer(C_INT), parameter  :: num_dim = 2

! --- Code

      name_len = len_trim(data_name)
      dims = shape(rdata)
      call simulation_data_read_float_f(nid,data_name,name_len,num_dim,dims,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_read_real4_rank2

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_real4_rank3(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),                          intent(in)  :: nid
      character(len=*),                     intent(in)  :: data_name
      real(kind=C_FLOAT), dimension(:,:,:), intent(out) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:3) :: dims

      integer(C_INT), parameter  :: num_dim = 3

! --- Code

      name_len = len_trim(data_name)
      dims = shape(rdata)
      call simulation_data_read_float_f(nid,data_name,name_len,num_dim,dims,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_read_real4_rank3

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_real8_rank0(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),         intent(in)  :: nid
      character(len=*),    intent(in)  :: data_name
      real(kind=C_DOUBLE), intent(out) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr


! --- Code

      name_len = len_trim(data_name)
      call simulation_data_read_double0_f(nid,data_name,name_len,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_read_real8_rank0

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_real8_rank1(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),                       intent(in)  :: nid
      character(len=*),                  intent(in)  :: data_name
      real(kind=C_DOUBLE), dimension(:), intent(out) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:1) :: dims

      integer(C_INT), parameter  :: num_dim = 1

! --- Code

      name_len = len_trim(data_name)
      dims = shape(rdata)
      call simulation_data_read_double_f(nid,data_name,name_len,num_dim,dims,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_read_real8_rank1

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_real8_rank2(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),                         intent(in)  :: nid
      character(len=*),                    intent(in)  :: data_name
      real(kind=C_DOUBLE), dimension(:,:), intent(out) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:2) :: dims

      integer(C_INT), parameter  :: num_dim = 2

! --- Code

      name_len = len_trim(data_name)
      dims = shape(rdata)
      call simulation_data_read_double_f(nid,data_name,name_len,num_dim,dims,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_read_real8_rank2

! ------------------------------------------------------------------------------

subroutine df90_sim_data_read_real8_rank3(nid,data_name,rdata,stat)

! --- Calling arguments

      type(C_PTR),                           intent(in)  :: nid
      character(len=*),                      intent(in)  :: data_name
      real(kind=C_DOUBLE), dimension(:,:,:), intent(out) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: name_len
      integer(C_INT) :: ierr
      integer(C_INT), dimension(1:3) :: dims

      integer(C_INT), parameter  :: num_dim = 3

! --- Code

      name_len = len_trim(data_name)
      dims = shape(rdata)
      call simulation_data_read_double_f(nid,data_name,name_len,num_dim,dims,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_sim_data_read_real8_rank3

! ------------------------------------------------------------------------------

subroutine df90_probe_data_count(sid,cnt,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      integer(kind=C_INT),               intent(out) :: cnt

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr

! --- code

      call probe_data_count_f(sid, cnt, ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_probe_data_count

! ------------------------------------------------------------------------------

subroutine df90_probe_data_dimensions(sid,pname,dlen,num,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      character(kind=C_CHAR,len=*),      intent(in)  :: pname
      integer(kind=C_INT),               intent(out) :: dlen
      integer(kind=C_INT),               intent(out) :: num

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: ierr, flen

! --- code

      flen=len(pname)
      call probe_data_dimensions_f(sid, pname,flen,dlen,num,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      end if 

end subroutine df90_probe_data_dimensions

! ------------------------------------------------------------------------------

subroutine df90_probe_data_list(sid,names,stat)

      type(C_PTR),                                intent(in)  :: sid
      character(kind=C_CHAR,len=*), dimension(:), intent(out) :: names

      integer, intent(out), optional :: stat

! --- local variables
      integer(C_INT) :: name_num, name_len, ierr

! --- code
      name_len = len(names)
      name_num = size(names)

      call probe_data_list_f(sid,names,name_len,name_num,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_list

! ------------------------------------------------------------------------------

subroutine df90_probe_data_open(sid,probe_name,pid,stat)

! --- Calling Arguments
 
      type(C_PTR),                       intent(in)  :: sid
      character(kind=C_CHAR,len=*),      intent(in)  :: probe_name
      type(C_PTR),                       intent(out) :: pid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, ierr

! --- code

      name_len = len(probe_name)
      call probe_data_open_f(sid, probe_name, name_len, pid, ierr)
      if (present(stat)) then
        call define_return_status(ierr,stat)
      end if 

end subroutine df90_probe_data_open

! ------------------------------------------------------------------------------

subroutine df90_probe_create_data_integer(sid,probe_name,idata,pid,stat)

! --- Calling Arguments
 
      type(C_PTR),                        intent(in)  :: sid
      character(kind=C_CHAR,len=*),       intent(in)  :: probe_name
      integer(kind=C_INT),dimension(:,:), intent(in)  :: idata
      type(C_PTR),                        intent(out) :: pid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, num, dlen, ierr

! --- code

      name_len = len(probe_name)
      dlen = size(idata,1)
      num = size(idata,2)
      call probe_create_data_int_f(sid, probe_name,name_len,dlen,num,idata,pid,ierr)
      if (present(stat)) then
        call define_return_status(ierr,stat)
      end if 


end subroutine df90_probe_create_data_integer

! ------------------------------------------------------------------------------

subroutine df90_probe_create_data_integer0(sid,probe_name,idata,pid,stat)

! --- Calling Arguments
 
      type(C_PTR),                        intent(in)  :: sid
      character(kind=C_CHAR,len=*),       intent(in)  :: probe_name
      integer(kind=C_INT),dimension(:),   intent(in)  :: idata
      type(C_PTR),                        intent(out) :: pid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, num, dlen, ierr

! --- code

      name_len = len(probe_name)
      dlen = 1
      num = size(idata,1)
      call probe_create_data_int_f(sid, probe_name,name_len,dlen,num,idata,pid,ierr)
      if (present(stat)) then
        call define_return_status(ierr,stat)
      end if 


end subroutine df90_probe_create_data_integer0

! ------------------------------------------------------------------------------

subroutine df90_probe_create_data_real4(sid,probe_name,rdata,pid,stat)

! --- Calling Arguments
 
      type(C_PTR),                        intent(in)  :: sid
      character(kind=C_CHAR,len=*),       intent(in)  :: probe_name
      real(kind=C_FLOAT), dimension(:,:), intent(in)  :: rdata
      type(C_PTR),                        intent(out) :: pid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, num, dlen, ierr

! --- code

      name_len = len(probe_name)
      dlen = size(rdata,1)
      num = size(rdata,2)
      call probe_create_data_float_f(sid, probe_name,name_len,dlen,num,rdata,pid,ierr)
      if (present(stat)) then
        call define_return_status(ierr,stat)
      end if 


end subroutine df90_probe_create_data_real4

! ------------------------------------------------------------------------------

subroutine df90_probe_create_data_real40(sid,probe_name,rdata,pid,stat)

! --- Calling Arguments
 
      type(C_PTR),                        intent(in)  :: sid
      character(kind=C_CHAR,len=*),       intent(in)  :: probe_name
      real(kind=C_FLOAT), dimension(:),   intent(in)  :: rdata
      type(C_PTR),                        intent(out) :: pid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, num, dlen, ierr

! --- code

      name_len = len(probe_name)
      dlen = 1
      num = size(rdata,1)
      call probe_create_data_float_f(sid, probe_name,name_len,dlen,num,rdata,pid,ierr)
      if (present(stat)) then
        call define_return_status(ierr,stat)
      end if 


end subroutine df90_probe_create_data_real40

! ------------------------------------------------------------------------------

subroutine df90_probe_create_data_real8(sid,probe_name,rdata,pid,stat)

! --- Calling Arguments
 
      type(C_PTR),                         intent(in)  :: sid
      character(kind=C_CHAR,len=*),        intent(in)  :: probe_name
      real(kind=C_DOUBLE), dimension(:,:), intent(in)  :: rdata
      type(C_PTR),                         intent(out) :: pid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, num, dlen, ierr

! --- code

      name_len = len(probe_name)
      dlen = size(rdata,1)
      num = size(rdata,2)
      call probe_create_data_double_f(sid, probe_name,name_len,dlen,num,rdata,pid,ierr)
      if (present(stat)) then
        call define_return_status(ierr,stat)
      end if 


end subroutine df90_probe_create_data_real8

! ------------------------------------------------------------------------------

subroutine df90_probe_create_data_real80(sid,probe_name,rdata,pid,stat)

! --- Calling Arguments
 
      type(C_PTR),                         intent(in)  :: sid
      character(kind=C_CHAR,len=*),        intent(in)  :: probe_name
      real(kind=C_DOUBLE), dimension(:),   intent(in)  :: rdata
      type(C_PTR),                         intent(out) :: pid

      integer, intent(out), optional :: stat

! --- Local variables
      integer(C_INT) :: name_len, num, dlen, ierr

! --- code

      name_len = len(probe_name)
      dlen = 1
      num = size(rdata,1)
      call probe_create_data_double_f(sid, probe_name,name_len,dlen,num,rdata,pid,ierr)
      if (present(stat)) then
        call define_return_status(ierr,stat)
      end if 


end subroutine df90_probe_create_data_real80

! ------------------------------------------------------------------------------

subroutine df90_probe_data_write_integer(pid,idata,stat)

! --- Calling arguments

      type(C_PTR),                         intent(in)  :: pid
      integer(kind=C_INT), dimension(:,:), intent(in)  :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: num, dlen
      integer(C_INT) :: ierr


! --- Code

      dlen=size(idata,1)
      num=size(idata,2)
      call probe_data_write_int_f(pid,num,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_write_integer

! ------------------------------------------------------------------------------

subroutine df90_probe_data_write_integer0(pid,idata,stat)

! --- Calling arguments

      type(C_PTR),                         intent(in)  :: pid
      integer(kind=C_INT), dimension(:),   intent(in)  :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: num, dlen
      integer(C_INT) :: ierr


! --- Code

      dlen=1
      num=size(idata,1)
      call probe_data_write_int_f(pid,num,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_write_integer0

! ------------------------------------------------------------------------------

subroutine df90_probe_data_write_real4(pid,rdata,stat)

! --- Calling arguments

      type(C_PTR),                         intent(in)  :: pid
      real(kind=C_FLOAT), dimension(:,:),  intent(in)  :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: num, dlen
      integer(C_INT) :: ierr


! --- Code

      dlen=size(rdata,1)
      num=size(rdata,2)
      call probe_data_write_float_f(pid,num,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_write_real4

! ------------------------------------------------------------------------------

subroutine df90_probe_data_write_real40(pid,rdata,stat)

! --- Calling arguments

      type(C_PTR),                         intent(in)  :: pid
      real(kind=C_FLOAT), dimension(:),    intent(in)  :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: num, dlen
      integer(C_INT) :: ierr


! --- Code

      dlen=1
      num=size(rdata,1)
      call probe_data_write_float_f(pid,num,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_write_real40

! ------------------------------------------------------------------------------

subroutine df90_probe_data_write_real8(pid,rdata,stat)

! --- Calling arguments

      type(C_PTR),                          intent(in)  :: pid
      real(kind=C_DOUBLE), dimension(:,:),  intent(in)  :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: num, dlen
      integer(C_INT) :: ierr


! --- Code

      dlen=size(rdata,1)
      num=size(rdata,2)
      call probe_data_write_double_f(pid,num,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_write_real8

! ------------------------------------------------------------------------------

subroutine df90_probe_data_write_real80(pid,rdata,stat)

! --- Calling arguments

      type(C_PTR),                          intent(in)  :: pid
      real(kind=C_DOUBLE), dimension(:),    intent(in)  :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: num, dlen
      integer(C_INT) :: ierr


! --- Code

      dlen=1
      num=size(rdata,1)
      call probe_data_write_double_f(pid,num,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_write_real80

! ------------------------------------------------------------------------------

subroutine df90_probe_data_read_integer0(pid,idata,stat)

! --- Calling arguments

      type(C_PTR),                       intent(in)  :: pid
      integer(kind=C_INT), dimension(:), intent(out) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer(C_INT) :: ierr

! --- Code
  
      call probe_data_read_int_f(pid,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_read_integer0

! ------------------------------------------------------------------------------

subroutine df90_probe_data_read_integer(pid,idata,stat)

! --- Calling arguments

      type(C_PTR),                         intent(in)  :: pid
      integer(kind=C_INT), dimension(:,:), intent(out) :: idata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: ierr

! --- Code

      call probe_data_read_int_f(pid,idata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_read_integer

! ------------------------------------------------------------------------------

subroutine df90_probe_data_read_real40(pid,rdata,stat)

! --- Calling arguments

      type(C_PTR),                       intent(in)  :: pid
      real(kind=C_FLOAT), dimension(:),  intent(out) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer(C_INT) :: ierr


! --- Code

      call probe_data_read_float_f(pid,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_read_real40

! ------------------------------------------------------------------------------

subroutine df90_probe_data_read_real4(pid,rdata,stat)

! --- Calling arguments

      type(C_PTR),                         intent(in)  :: pid
      real(kind=C_FLOAT), dimension(:,:),  intent(out) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: ierr

! --- Code
    
      call probe_data_read_float_f(pid,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_read_real4

! ------------------------------------------------------------------------------

subroutine df90_probe_data_read_real80(pid,rdata,stat)

! --- Calling arguments

      type(C_PTR),                        intent(in)  :: pid
      real(kind=C_DOUBLE), dimension(:),  intent(out) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer(C_INT) :: ierr

! --- Code

      call probe_data_read_double_f(pid,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_read_real80

! ------------------------------------------------------------------------------

subroutine df90_probe_data_read_real8(pid,rdata,stat)

! --- Calling arguments

      type(C_PTR),                          intent(in)  :: pid
      real(kind=C_DOUBLE), dimension(:,:),  intent(out) :: rdata            

      integer, intent(out), optional :: stat

! --- Local variables

      integer :: ierr

! --- Code

      call probe_data_read_double_f(pid,rdata,ierr)
      if (present(stat)) then
          call define_return_status(ierr,stat)
      endif       

end subroutine df90_probe_data_read_real8

! ------------------------------------------------------------------------------


!! =============================================================================
!! PRIVATE
!! =============================================================================
!
subroutine define_return_status(err,stat)

! --- calling arguments
      integer(C_INT), intent(in)  :: err
      integer,           intent(out) :: stat


      if ( err .ne. 0 ) then
          stat = Danu_FAILURE
      else
          stat = Danu_SUCCESS
      end if


end subroutine define_return_status      
! ==============================================================================
end module danu_module
! ==============================================================================

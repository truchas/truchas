! ==============================================================================
!                                                                              !
!                                                                              !
!                             Copyright  (C) 20xx,                             !
!                      Los Alamos National Security, LLC                       !
!                                                                              !
!                             LA-CC-xxxxxx                                     !
!                                                                              !
! ==============================================================================
!
!  DANU
!     Interface module 
!     Defines interface blocks for C functions and subroutines
!     See danu_fort_*.h for the C prototypes

!
module danu_iface

! ---
     implicit none

! ==============================================================================
! Enumerated Types
! ==============================================================================

    
! ==============================================================================
! Interfaces
! ==============================================================================
    
! --- File Control
interface 
    subroutine danu_file_create_f(filename,length,ptr,ierr) bind(c)
    use, intrinsic :: iso_c_binding
    character(C_CHAR)  :: filename(*)
    integer(C_INT)     :: length
    type(C_PTR)        :: ptr
    integer(C_INT)     :: ierr
    end subroutine danu_file_create_f
end interface

interface 
    subroutine danu_file_open_rdonly_f(filename,length,ptr,ierr) bind(c)
    use, intrinsic :: iso_c_binding
    character(C_CHAR)  :: filename(*)
    integer(C_INT)     :: length
    type(C_PTR)        :: ptr
    integer(C_INT)     :: ierr
    end subroutine danu_file_open_rdonly_f
end interface

interface 
    subroutine danu_file_open_rdwr_f(filename,length,ptr,ierr) bind(c)
    use, intrinsic :: iso_c_binding
    character(C_CHAR)  :: filename(*)
    integer(C_INT)     :: length
    type(C_PTR)        :: ptr
    integer(C_INT)     :: ierr
    end subroutine danu_file_open_rdwr_f
end interface

interface 
    subroutine danu_file_close_f(ptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                     intent(inout) :: ptr
    integer(C_INT),                  intent(out)   :: ierr
    end subroutine danu_file_close_f
end interface 

interface 
    subroutine output_file_is_valid_f(ptr,flag) bind(c)
    use iso_c_binding
    type(C_PTR),                     intent(in)    :: ptr
    integer(C_INT),                  intent(out)   :: flag
    end subroutine output_file_is_valid_f
end interface

interface 
    subroutine output_file_create_f(filename,length,ptr,ierr) bind(c)
    use iso_c_binding
    character(C_CHAR),               intent(in)  :: filename(*)
    integer(C_INT),                  intent(in)  :: length
    type(C_PTR),                     intent(out) :: ptr
    integer(C_INT),                  intent(out) :: ierr
    end subroutine output_file_create_f
end interface

interface 
    subroutine output_file_open_rdonly_f(filename,length,ptr,ierr) bind(c)
    use iso_c_binding
    character(C_CHAR),               intent(in)  :: filename(*)
    integer(C_INT),                  intent(in)  :: length
    type(C_PTR),                     intent(out) :: ptr
    integer(C_INT),                  intent(out) :: ierr
    end subroutine output_file_open_rdonly_f
end interface

interface 
    subroutine output_file_open_rdwr_f(filename,length,ptr,ierr) bind(c)
    use iso_c_binding
    character(C_CHAR),               intent(in)  :: filename(*)
    integer(C_INT),                  intent(in)  :: length
    type(C_PTR),                     intent(out) :: ptr
    integer(C_INT),                  intent(out) :: ierr
    end subroutine output_file_open_rdwr_f
end interface

interface 
    subroutine output_file_open_append_f(filename,length,ptr,ierr) bind(c)
    use iso_c_binding
    character(C_CHAR),               intent(in)  :: filename(*)
    integer(C_INT),                  intent(in)  :: length
    type(C_PTR),                     intent(out) :: ptr
    integer(C_INT),                  intent(out) :: ierr
    end subroutine output_file_open_append_f
end interface

interface 
    subroutine output_file_close_f(ptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                     intent(inout) :: ptr
    integer(C_INT),                  intent(out)   :: ierr
    end subroutine output_file_close_f
end interface

! --- Attribute Control
interface
    subroutine danu_attr_write_int_f(ptr,attr_name,name_len,value,flag) bind(c)
    use iso_c_binding
    type(C_PTR),                     intent(in)    :: ptr
    character(kind=C_CHAR),          intent(in)    :: attr_name(*)
    integer(kind=C_INT),             intent(in)    :: name_len
    integer(kind=C_INT),             intent(in)    :: value
    integer(kind=C_INT),             intent(out)   :: flag
    end subroutine danu_attr_write_int_f
end interface

interface
    subroutine danu_attr_read_int_f(ptr,attr_name,name_len,buffer,flag) bind(c)
    use iso_c_binding
    type(C_PTR),                     intent(in)    :: ptr
    character(kind=C_CHAR),          intent(in)    :: attr_name(*)
    integer(kind=C_INT),             intent(in)    :: name_len
    integer(kind=C_INT),             intent(out)   :: buffer
    integer(kind=C_INT),             intent(out)   :: flag
    end subroutine danu_attr_read_int_f
end interface

interface
    subroutine danu_attr_write_real4_f(ptr,attr_name,name_len,value,flag) bind(c)
    use iso_c_binding
    type(C_PTR),                     intent(in)    :: ptr
    character(kind=C_CHAR),          intent(in)    :: attr_name(*)
    integer(kind=C_INT),             intent(in)    :: name_len
    real(kind=C_FLOAT),              intent(in)    :: value
    integer(kind=C_INT),             intent(out)   :: flag
    end subroutine danu_attr_write_real4_f
end interface

interface
    subroutine danu_attr_read_real4_f(ptr,attr_name,name_len,buffer,flag) &
                                 bind(c)
    use iso_c_binding
    type(C_PTR),                     intent(in)    :: ptr
    character(kind=C_CHAR),          intent(in)    :: attr_name(*)
    integer(kind=C_INT),             intent(in)    :: name_len
    real(kind=C_FLOAT),              intent(out)   :: buffer
    integer(kind=C_INT),             intent(out)   :: flag
    end subroutine danu_attr_read_real4_f
end interface

interface
    subroutine danu_attr_write_real8_f(ptr,attr_name,name_len,value,flag) &
                                       bind(c)
    use iso_c_binding
    type(C_PTR),                     intent(in)    :: ptr
    character(kind=C_CHAR),          intent(in)    :: attr_name(*)
    integer(kind=C_INT),             intent(in)    :: name_len
    real(kind=C_DOUBLE),             intent(in)    :: value
    integer(kind=C_INT),             intent(out)   :: flag
    end subroutine danu_attr_write_real8_f
end interface

interface
    subroutine danu_attr_read_real8_f(ptr,attr_name,name_len,buffer,flag) &
                                 bind(c)
    use iso_c_binding
    type(C_PTR),                     intent(in)    :: ptr
    character(kind=C_CHAR),          intent(in)    :: attr_name(*)
    integer(kind=C_INT),             intent(in)    :: name_len
    real(kind=C_DOUBLE),             intent(out)   :: buffer
    integer(kind=C_INT),             intent(out)   :: flag
    end subroutine danu_attr_read_real8_f
end interface

interface
    subroutine danu_attr_write_char_f(ptr,          &
                                      attr_name,    &
                                      name_len,     &
                                      char_data,    &
                                      length,       &
                                      flag)         &
                                      bind(c)
    use iso_c_binding
    type(C_PTR),                     intent(in)    :: ptr
    character(kind=C_CHAR),          intent(in)    :: attr_name(*)
    integer(kind=C_INT),             intent(in)    :: name_len
    character(kind=C_CHAR),          intent(in)    :: char_data(*)
    integer(kind=C_INT),             intent(in)    :: length
    integer(kind=C_INT),             intent(out)   :: flag
    end subroutine danu_attr_write_char_f
end interface

interface
    subroutine danu_attr_read_char_f(ptr,          &
                                      attr_name,    &
                                      name_len,     &
                                      char_data,    &
                                      length,       &
                                      flag)         &
                                      bind(c)
    use iso_c_binding
    type(C_PTR),                     intent(in)    :: ptr
    character(kind=C_CHAR),          intent(in)    :: attr_name(*)
    integer(kind=C_INT),             intent(in)    :: name_len
    character(kind=C_CHAR),          intent(out)   :: char_data(*)
    integer(kind=C_INT),             intent(in)    :: length
    integer(kind=C_INT),             intent(out)   :: flag
    end subroutine danu_attr_read_char_f
end interface

interface 
    subroutine danu_attr_exists_f(ptr,attr_name,name_len,exists,flag) &
                                  bind(c)
    use iso_c_binding
    type(C_PTR),                 intent(in)  :: ptr
    character(kind=C_CHAR),      intent(in)  :: attr_name(*)
    integer(kind=C_INT),         intent(in)  :: name_len
    integer(kind=C_INT),         intent(out) :: exists
    integer(kind=C_INT),         intent(out) :: flag
    end subroutine danu_attr_exists_f
end interface    
 
interface
    subroutine danu_attr_count_f(ptr,attr_count,flag) bind(c)
    use iso_c_binding
    type(C_PTR),         intent(in)  :: ptr
    integer(kind=C_INT), intent(out) :: attr_count
    integer(kind=C_INT), intent(out) :: flag
    end subroutine danu_attr_count_f
end interface    

interface
    subroutine danu_attr_names_f(ptr,names,name_len,name_num,flag) &
                                 bind(c)
    use iso_c_binding
    type(C_PTR),                               intent(in)  :: ptr
    character(kind=C_CHAR),                    intent(out) :: names(*)
    integer(kind=C_INT),                       intent(in)  :: name_len
    integer(kind=C_INT),                       intent(in)  :: name_num
    integer(kind=C_INT),                       intent(out) :: flag
    end subroutine danu_attr_names_f
end interface    

! --- Group control

interface
    subroutine danu_group_exists_f(ptr,grp_name,nlen, exists, ierr) bind(c)
    use iso_c_binding
    type(C_PTR),            intent(in) :: ptr
    character(kind=C_CHAR), intent(in) :: grp_name
    integer(kind=C_INT),    intent(in) :: nlen
    integer(kind=C_INT),    intent(out) :: exists
    integer(kind=C_INT),    intent(out) :: ierr
    end subroutine danu_group_exists_f
end interface

interface
    subroutine danu_group_create_f(ptr,grp_name,nlen,group,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),            intent(in)  :: ptr
    character(kind=C_CHAR), intent(in)  :: grp_name(*)
    integer(kind=C_INT),    intent(in)  :: nlen
    type(C_PTR),            intent(out) :: group
    integer(kind=C_INT),    intent(out) :: ierr
    end subroutine danu_group_create_f
end interface

interface
    subroutine danu_group_open_f(ptr,grp_name,nlen,group,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),            intent(in)  :: ptr
    character(kind=C_CHAR), intent(in)  :: grp_name(*)
    integer(kind=C_INT),    intent(in)  :: nlen
    type(C_PTR),            intent(out) :: group
    integer(kind=C_INT),    intent(out) :: ierr
    end subroutine danu_group_open_f
end interface

interface
    subroutine danu_group_close_f(group,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),            intent(inout) :: group
    integer(kind=C_INT),    intent(out)   :: ierr
    end subroutine danu_group_close_f
end interface

! --- Meshes

interface
    subroutine mesh_count_f(ptr,mesh_count,flag) bind(c)
    use iso_c_binding
    type(C_PTR),         intent(in)  :: ptr
    integer(kind=C_INT), intent(out) :: mesh_count
    integer(kind=C_INT), intent(out) :: flag
    end subroutine mesh_count_f
end interface   

interface
    subroutine mesh_list_f(ptr,names,name_len,name_size,flag) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(out) :: names(*)
    integer(kind=C_INT),                        intent(in)  :: name_len
    integer(kind=C_INT),                        intent(in)  :: name_size
    integer(kind=C_INT),                        intent(out) :: flag
    end subroutine mesh_list_f
end interface    

interface
    subroutine mesh_exists_f(ptr,mesh_name,name_len,flag,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: ptr
    character(kind=C_CHAR),       intent(in)  :: mesh_name(*)
    integer(kind=C_INT),          intent(in)  :: name_len
    integer(kind=C_INT),          intent(out) :: flag
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_exists_f
end interface    

interface
    subroutine mesh_open_f(ptr,mesh_name,name_len,mid,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: ptr
    character(kind=C_CHAR),       intent(in)  :: mesh_name(*)
    integer(kind=C_INT),          intent(in)  :: name_len
    type(C_PTR),                  intent(out) :: mid
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_open_f
end interface    

interface
    subroutine mesh_add_unstructured_f(ptr,                                   &
                                      mesh_name,                              &
                                      name_len,                               &
                                      elemorder,                              &
                                      mesh_dim,                               &
                                      mid,                                    &
                                      ierr)                                   &
                                      bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: ptr
    character(kind=C_CHAR),       intent(in)  :: mesh_name(*)
    integer(kind=C_INT),          intent(in)  :: name_len
    integer(kind=C_INT),          intent(in)  :: elemorder
    integer(kind=C_INT),          intent(in)  :: mesh_dim
    type(C_PTR),                  intent(out) :: mid
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_add_unstructured_f
end interface    

interface
    subroutine mesh_write_coordinates_f(mptr,                                 &
                                        nnodes,                               &
                                        x,                                    &
                                        y,                                    &
                                        z,                                    &
                                        ierr)                                 &
                                        bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    integer(kind=C_INT),          intent(in)  :: nnodes
    real(kind=C_DOUBLE),          intent(in)  :: x(*)
    real(kind=C_DOUBLE),          intent(in)  :: y(*)
    real(kind=C_DOUBLE),          intent(in)  :: z(*)
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_write_coordinates_f
end interface 

interface
    subroutine mesh_write_coordinates_1d_f(mptr,                                 &
                                           nnodes,                               &
                                           x,                                    &
                                           ierr)                                 &
                                           bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    integer(kind=C_INT),          intent(in)  :: nnodes
    real(kind=C_DOUBLE),          intent(in)  :: x(*)
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_write_coordinates_1d_f
end interface 

interface
    subroutine mesh_write_coordinates_2d_f(mptr,                                 &
                                           nnodes,                               &
                                           x,                                    &
                                           y,                                    &
                                           ierr)                                 &
                                           bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    integer(kind=C_INT),          intent(in)  :: nnodes
    real(kind=C_DOUBLE),          intent(in)  :: x(*)
    real(kind=C_DOUBLE),          intent(in)  :: y(*)
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_write_coordinates_2d_f
end interface 

interface
    subroutine mesh_read_coordinates_f(mptr,                                 &
                                        x,                                   &
                                        y,                                   &
                                        z,                                   &
                                        ierr)                                &
                                        bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    real(kind=C_DOUBLE),          intent(out)  :: x(*)
    real(kind=C_DOUBLE),          intent(out)  :: y(*)
    real(kind=C_DOUBLE),          intent(out)  :: z(*)
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_read_coordinates_f
end interface 

interface
    subroutine mesh_read_coordinates_byindex_f(mptr,                          &
                                               idx,                           &
                                               buf,                           &
                                               ierr)                          &
                                               bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    integer(kind=C_INT),          intent(in)  :: idx
    real(kind=C_DOUBLE),          intent(out) :: buf(*)
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_read_coordinates_byindex_f
end interface

interface
    subroutine mesh_read_coordinates_1d_f(mptr,                                 &
                                          x,                                    &
                                          ierr)                                 &
                                          bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    real(kind=C_DOUBLE),          intent(out)  :: x(*)
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_read_coordinates_1d_f
end interface 

interface
    subroutine mesh_read_coordinates_2d_f(mptr,                                 &
                                          x,                                    &
                                          y,                                    &
                                          ierr)                                 &
                                          bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    real(kind=C_DOUBLE),          intent(out) :: x(*)
    real(kind=C_DOUBLE),          intent(out) :: y(*)
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_read_coordinates_2d_f
end interface 

interface
    subroutine mesh_write_connectivity_f(mptr,                            & 
                                         num,                             &
                                         idata,                           &
                                         ierr)                            &
                                         bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    integer(kind=C_INT),          intent(in)  :: num
    integer(kind=C_INT),          intent(in)  :: idata(*)
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_write_connectivity_f
end interface 

interface
    subroutine mesh_read_connectivity_f(mptr,                             & 
                                        idata,                            &
                                        ierr)                             &
                                        bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    integer(kind=C_INT),          intent(out) :: idata(*)
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_read_connectivity_f
end interface 

interface
  subroutine mesh_connectivity_size_f(mptr,dataname,flen,isize,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    character(kind=C_CHAR),       intent(in)  :: dataname
    integer(kind=C_INT),          intent(in)  :: flen
    integer(kind=C_INT),          intent(out) :: isize
    integer(kind=C_INT),          intent(out) :: ierr
  end subroutine mesh_connectivity_size_f
end interface

interface
    subroutine mesh_get_type_f(mptr,mesh_type,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    integer(kind=C_INT),          intent(out) :: mesh_type
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_get_type_f
end interface

interface
    subroutine mesh_get_elementtype_f(mptr,elem_type,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    integer(kind=C_INT),          intent(out) :: elem_type
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_get_elementtype_f
end interface

interface
    subroutine mesh_get_dimension_f(mptr,d,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    integer(kind=C_INT),          intent(out) :: d
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_get_dimension_f
end interface

interface
    subroutine mesh_get_nnodes_f(mptr,nnodes,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    integer(kind=C_INT),          intent(out) :: nnodes
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_get_nnodes_f
end interface

interface
    subroutine mesh_get_nelem_f(mptr,nelem,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    integer(kind=C_INT),          intent(out) :: nelem
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_get_nelem_f
end interface

interface
    subroutine mesh_get_elem_order_f(mptr,elem_order,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: mptr
    integer(kind=C_INT),          intent(out) :: elem_order
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine mesh_get_elem_order_f
end interface

interface
  subroutine mesh_create_hex_unstruct_f(ptr,                           &
                                        mname,                         &
                                        flen,                          &
                                        nnodes,x,y,z,                  &
                                        nelem,conn,                    &
                                        mptr,                          &
                                        ierr)                          &
                                        bind(c)
  use iso_c_binding
  type(C_PTR),                  intent(in)  :: ptr
  character(kind=C_CHAR),       intent(in)  :: mname(*)
  integer(kind=C_INT),          intent(in)  :: flen
  integer(kind=C_INT),          intent(in)  :: nnodes
  real(kind=C_DOUBLE),          intent(in)  :: x(*)
  real(kind=C_DOUBLE),          intent(in)  :: y(*)
  real(kind=C_DOUBLE),          intent(in)  :: z(*)
  integer(kind=C_INT),          intent(in)  :: nelem
  integer(kind=C_INT),          intent(in)  :: conn(*)
  type(C_PTR),                  intent(out) :: mptr
  integer(kind=C_INT),          intent(out) :: ierr
  end subroutine mesh_create_hex_unstruct_f
end interface 

! --- Simulations

interface
    subroutine simulation_exists_f(fptr,sim_name,flen,flag,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: fptr
    character(kind=C_CHAR),       intent(in)  :: sim_name
    integer(kind=C_INT),          intent(in)  :: flen
    integer(kind=C_INT),          intent(out) :: flag
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine simulation_exists_f
end interface 

interface
    subroutine simulation_count_f(fptr,cnt,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: fptr
    integer(kind=C_INT),          intent(out) :: cnt
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine simulation_count_f
end interface 

interface
    subroutine simulation_list_f(ptr,names,name_len,name_size,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(out) :: names(*)
    integer(kind=C_INT),                        intent(in)  :: name_len
    integer(kind=C_INT),                        intent(in)  :: name_size
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_list_f
end interface  

interface
    subroutine simulation_add_f(fptr,sim_name,flen,sptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: fptr
    character(kind=C_CHAR),       intent(in)  :: sim_name(*)
    integer(kind=C_INT),          intent(in)  :: flen
    type(C_PTR),                  intent(out) :: sptr
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine simulation_add_f
end interface 

interface
    subroutine simulation_open_f(fptr,sim_name,flen,sptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: fptr
    character(kind=C_CHAR),       intent(in)  :: sim_name(*)
    integer(kind=C_INT),          intent(in)  :: flen
    type(C_PTR),                  intent(out) :: sptr
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine simulation_open_f
end interface

interface
    subroutine simulation_link_mesh_f(fptr,sptr,mesh_name,flen,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: fptr
    type(C_PTR),                  intent(in)  :: sptr
    character(kind=C_CHAR),       intent(in)  :: mesh_name
    integer(kind=C_INT),          intent(in)  :: flen
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine simulation_link_mesh_f
end interface

! --- Non-series Datasets

interface
    subroutine data_open_dataset_f(fptr,data_name,flen,hptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: fptr
    character(kind=C_CHAR),       intent(in)  :: data_name
    integer(kind=C_INT),          intent(in)  :: flen
    type(C_PTR),                  intent(out) :: hptr
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine data_open_dataset_f
end interface 

interface
    subroutine data_exists_f(fptr,data_name,flen,flag,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: fptr
    character(kind=C_CHAR),       intent(in)  :: data_name
    integer(kind=C_INT),          intent(in)  :: flen
    integer(kind=C_INT),          intent(out) :: flag
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine data_exists_f
end interface 

interface
    subroutine data_count_f(fptr,cnt,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: fptr
    integer(kind=C_INT),          intent(out) :: cnt
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine data_count_f
end interface

interface
  subroutine data_type_f(fptr,data_name,flen,typecode,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: fptr
    character(kind=C_CHAR),       intent(in)  :: data_name
    integer(kind=C_INT),          intent(in)  :: flen
    integer(kind=C_INT),          intent(out) :: typecode
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine data_type_f
end interface

interface
    subroutine data_rank_f(fptr,data_name,flen,rank,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: fptr
    character(kind=C_CHAR),       intent(in)  :: data_name
    integer(kind=C_INT),          intent(in)  :: flen
    integer(kind=C_INT),          intent(out) :: rank
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine data_rank_f
end interface

interface
    subroutine data_dimensions_f(fptr,data_name,flen,rank,dimensions,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                  intent(in)  :: fptr
    character(kind=C_CHAR),       intent(in)  :: data_name
    integer(kind=C_INT),          intent(in)  :: flen
    integer(kind=C_INT),          intent(in)  :: rank
    integer(kind=C_INT),          intent(out) :: dimensions(*)
    integer(kind=C_INT),          intent(out) :: ierr
    end subroutine data_dimensions_f
end interface

interface
    subroutine data_list_f(ptr,names,name_len,name_size,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(out) :: names(*)
    integer(kind=C_INT),                        intent(in)  :: name_len
    integer(kind=C_INT),                        intent(in)  :: name_size
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_list_f
end interface  

interface
    subroutine data_write_int0_f(ptr,data_name,flen,idata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: idata
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_write_int0_f
end interface

interface
    subroutine data_write_int_f(ptr,data_name,flen,ndim,dims,idata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(in)  :: dims(*)
    integer(kind=C_INT),                        intent(in)  :: idata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_write_int_f
end interface

interface
    subroutine data_write_float0_f(ptr,data_name,flen,r4data,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    real(kind=C_FLOAT),                         intent(in)  :: r4data
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_write_float0_f
end interface

interface
    subroutine data_write_float_f(ptr,data_name,flen,ndim,dims,r4data,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(in)  :: dims(*)
    real(kind=C_FLOAT),                         intent(in)  :: r4data(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_write_float_f
end interface

interface
    subroutine data_write_double0_f(ptr,data_name,flen,r8data,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    real(kind=C_DOUBLE),                        intent(in)  :: r8data
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_write_double0_f
end interface

interface
    subroutine data_write_double_f(ptr,data_name,flen,ndim,dims,r8data,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(in)  :: dims(*)
    real(kind=C_DOUBLE),                        intent(in)  :: r8data(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_write_double_f
end interface

interface
    subroutine data_read_int0_f(ptr,data_name,flen,idata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(out) :: idata
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_read_int0_f
end interface

interface
    subroutine data_read_int_f(ptr,data_name,flen,ndim,dims,idata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(in)  :: dims(*)
    integer(kind=C_INT),                        intent(out) :: idata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_read_int_f
end interface

interface
    subroutine data_read_float0_f(ptr,data_name,flen,r4data,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    real(kind=C_FLOAT),                         intent(out) :: r4data
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_read_float0_f
end interface

interface
    subroutine data_read_float_f(ptr,data_name,flen,ndim,dims,r4data,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(in)  :: dims(*)
    real(kind=C_FLOAT),                         intent(out) :: r4data(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_read_float_f
end interface

interface
    subroutine data_read_double0_f(ptr,data_name,flen,r8data,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    real(kind=C_DOUBLE),                        intent(out) :: r8data
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_read_double0_f
end interface

interface
    subroutine data_read_double_f(ptr,data_name,flen,ndim,dims,r8data,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(in)  :: dims(*)
    real(kind=C_DOUBLE),                        intent(out) :: r8data(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine data_read_double_f
end interface

! --- Series Group access and control

interface
    subroutine sequence_exists_f(sptr, sname, slen, flag, ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: sptr
    character(kind=C_CHAR),                     intent(in)  :: sname(*)
    integer(kind=C_INT),                        intent(in)  :: slen
    integer(kind=C_INT),                        intent(out) :: flag
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine sequence_exists_f
end interface    

interface
    subroutine sequence_count_f(sptr, num, ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: sptr
    integer(kind=C_INT),                        intent(out) :: num
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine sequence_count_f
end interface 

interface 
    subroutine sequence_get_next_id_f(sptr,cyc,time,nptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                 intent(in)  :: sptr
    integer(kind=C_INT),                         intent(in)  :: cyc
    real(kind=C_DOUBLE),                         intent(in)  :: time
    type(C_PTR),                                 intent(out) :: nptr
    integer(kind=C_INT),                         intent(out) :: ierr
    end subroutine sequence_get_next_id_f
end interface 

interface
    subroutine sequence_list_f(ptr,names,name_len,name_size,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(out) :: names(*)
    integer(kind=C_INT),                        intent(in)  :: name_len
    integer(kind=C_INT),                        intent(in)  :: name_size
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine sequence_list_f
end interface  

interface
    subroutine sequence_get_id_f(sptr,sname,slen,nptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: sptr
    character(kind=C_CHAR),                     intent(in)  :: sname(*)
    integer(kind=C_INT),                        intent(in)  :: slen
    type(C_PTR),                                intent(out) :: nptr
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine sequence_get_id_f
end interface

interface
    subroutine sequence_get_id_bynum_f(sptr,num,nptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: sptr
    integer(kind=C_INT),                        intent(in)  :: num
    type(C_PTR),                                intent(out) :: nptr
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine sequence_get_id_bynum_f
end interface

! --- Series Data access and control

interface
    subroutine simulation_open_data_f(nptr,dname,dlen,dptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: nptr
    character(kind=C_CHAR),                     intent(in)  :: dname(*)
    integer(kind=C_INT),                        intent(in)  :: dlen
    type(C_PTR),                                intent(out) :: dptr
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_open_data_f
end interface    

interface
    subroutine simulation_data_exists_f(nptr,dname,dlen,flag,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: nptr
    character(kind=C_CHAR),                     intent(in)  :: dname(*)
    integer(kind=C_INT),                        intent(in)  :: dlen
    integer(kind=C_INT),                        intent(out) :: flag
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_exists_f
end interface    

interface
    subroutine simulation_data_count_f(nptr,cnt,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: nptr
    integer(kind=C_INT),                        intent(out) :: cnt
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_count_f
end interface

interface
    subroutine simulation_data_rank_f(nptr,dname,flen,rank,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: nptr
    character(kind=C_CHAR),                     intent(in)  :: dname(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(out) :: rank 
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_rank_f
end interface

interface
    subroutine simulation_data_type_f(nptr,dname,flen,typecode,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: nptr
    character(kind=C_CHAR),                     intent(in)  :: dname(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(out) :: typecode
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_type_f
end interface

interface
    subroutine simulation_data_dimensions_f(nptr,dname,flen,ndim,dims,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: nptr
    character(kind=C_CHAR),                     intent(in)  :: dname(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(out) :: dims(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_dimensions_f
end interface

interface
    subroutine simulation_data_list_f(ptr,names,name_len,name_size,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(out) :: names(*)
    integer(kind=C_INT),                        intent(in)  :: name_len
    integer(kind=C_INT),                        intent(in)  :: name_size
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_list_f
end interface  


interface
    subroutine simulation_data_write_int0_f(ptr,data_name,flen,idata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: idata
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_write_int0_f
end interface

interface
    subroutine simulation_data_write_int_f(ptr,data_name,flen,ndim,dims,idata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(in)  :: dims(*)
    integer(kind=C_INT),                        intent(in)  :: idata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_write_int_f
end interface

interface
    subroutine simulation_data_write_float0_f(ptr,data_name,flen,rdata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    real(kind=C_FLOAT),                         intent(in)  :: rdata
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_write_float0_f
end interface

interface
    subroutine simulation_data_write_float_f(ptr,data_name,flen,ndim,dims,rdata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(in)  :: dims(*)
    real(kind=C_FLOAT),                         intent(in)  :: rdata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_write_float_f
end interface

interface
    subroutine simulation_data_write_double0_f(ptr,data_name,flen,rdata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    real(kind=C_DOUBLE),                        intent(in)  :: rdata
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_write_double0_f
end interface

interface
    subroutine simulation_data_write_double_f(ptr,data_name,flen,ndim,dims,rdata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(in)  :: dims(*)
    real(kind=C_DOUBLE),                        intent(in)  :: rdata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_write_double_f
end interface


interface
    subroutine simulation_data_read_int0_f(ptr,data_name,flen,idata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(out) :: idata
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_read_int0_f
end interface

interface
    subroutine simulation_data_read_int_f(ptr,data_name,flen,ndim,dims,idata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(in)  :: dims(*)
    integer(kind=C_INT),                        intent(out) :: idata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_read_int_f
end interface

interface
    subroutine simulation_data_read_float0_f(ptr,data_name,flen,rdata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    real(kind=C_FLOAT),                         intent(out) :: rdata
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_read_float0_f
end interface

interface
    subroutine simulation_data_read_float_f(ptr,data_name,flen,ndim,dims,rdata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(in)  :: dims(*)
    real(kind=C_FLOAT),                         intent(out) :: rdata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_read_float_f
end interface

interface
    subroutine simulation_data_read_double0_f(ptr,data_name,flen,rdata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    real(kind=C_DOUBLE),                        intent(out) :: rdata
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_read_double0_f
end interface

interface
    subroutine simulation_data_read_double_f(ptr,data_name,flen,ndim,dims,rdata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: ptr
    character(kind=C_CHAR),                     intent(in)  :: data_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: ndim
    integer(kind=C_INT),                        intent(in)  :: dims(*)
    real(kind=C_DOUBLE),                        intent(out) :: rdata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine simulation_data_read_double_f
end interface

! --- Probes

interface
    subroutine probe_data_exists_f(sptr,probe_name,flen,flag,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: sptr
    character(kind=C_CHAR),                     intent(in)  :: probe_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(out) :: flag
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_data_exists_f
end interface

interface
    subroutine probe_data_count_f(sptr,cnt,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: sptr
    integer(kind=C_INT),                        intent(out) :: cnt
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_data_count_f
end interface

interface
    subroutine probe_data_dimensions_f(sptr,pname,flen,dlen,num,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: sptr
    character(kind=C_CHAR),                     intent(in)  :: pname(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(out) :: dlen
    integer(kind=C_INT),                        intent(out) :: num
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_data_dimensions_f
end interface

interface
    subroutine probe_data_dimensions2_f(pptr,dlen,num,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: pptr
    integer(kind=C_INT),                        intent(out) :: dlen
    integer(kind=C_INT),                        intent(out) :: num
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_data_dimensions2_f
end interface

interface
    subroutine probe_data_list_f(sptr,names,name_len,name_size,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: sptr
    character(kind=C_CHAR),                     intent(out) :: names(*)
    integer(kind=C_INT),                        intent(in)  :: name_len
    integer(kind=C_INT),                        intent(in)  :: name_size
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_data_list_f
end interface  

interface
    subroutine probe_data_open_f(sptr,probe_name,flen,pptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: sptr
    character(kind=C_CHAR),                     intent(in)  :: probe_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    type(C_PTR),                                intent(out) :: pptr
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_data_open_f
end interface

interface
    subroutine probe_create_data_int_f(sptr,probe_name,                       &
                                       flen,dlen,num,idata,pptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: sptr
    character(kind=C_CHAR),                     intent(in)  :: probe_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: dlen
    integer(kind=C_INT),                        intent(in)  :: num
    integer(kind=C_INT),                        intent(in)  :: idata(*)
    type(C_PTR),                                intent(out) :: pptr
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_create_data_int_f
end interface

interface
    subroutine probe_create_data_float_f(sptr,probe_name,                       &
                                       flen,dlen,num,rdata,pptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: sptr
    character(kind=C_CHAR),                     intent(in)  :: probe_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: dlen
    integer(kind=C_INT),                        intent(in)  :: num
    real(kind=C_FLOAT),                         intent(in)  :: rdata(*)
    type(C_PTR),                                intent(out) :: pptr
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_create_data_float_f
end interface

interface
    subroutine probe_create_data_double_f(sptr,probe_name,                       &
                                       flen,dlen,num,rdata,pptr,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: sptr
    character(kind=C_CHAR),                     intent(in)  :: probe_name(*)
    integer(kind=C_INT),                        intent(in)  :: flen
    integer(kind=C_INT),                        intent(in)  :: dlen
    integer(kind=C_INT),                        intent(in)  :: num
    real(kind=C_DOUBLE),                        intent(in)  :: rdata(*)
    type(C_PTR),                                intent(out) :: pptr
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_create_data_double_f
end interface

interface
    subroutine probe_data_write_int_f(pptr,num,idata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: pptr
    integer(kind=C_INT),                        intent(in)  :: num
    integer(kind=C_INT),                        intent(in)  :: idata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_data_write_int_f
end interface

interface
    subroutine probe_data_write_float_f(pptr,num,rdata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: pptr
    integer(kind=C_INT),                        intent(in)  :: num
    real(kind=C_FLOAT),                         intent(in)  :: rdata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_data_write_float_f
end interface

interface
    subroutine probe_data_write_double_f(pptr,num,rdata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: pptr
    integer(kind=C_INT),                        intent(in)  :: num
    real(kind=C_DOUBLE),                        intent(in)  :: rdata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_data_write_double_f
end interface

interface
    subroutine probe_data_read_int_f(pptr,idata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: pptr
    integer(kind=C_INT),                        intent(out) :: idata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_data_read_int_f
end interface

interface
    subroutine probe_data_read_float_f(pptr,rdata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: pptr
    real(kind=C_FLOAT),                         intent(out) :: rdata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_data_read_float_f
end interface

interface
    subroutine probe_data_read_double_f(pptr,rdata,ierr) bind(c)
    use iso_c_binding
    type(C_PTR),                                intent(in)  :: pptr
    real(kind=C_DOUBLE),                        intent(out) :: rdata(*)
    integer(kind=C_INT),                        intent(out) :: ierr
    end subroutine probe_data_read_double_f
end interface


! =============================================================================
! =============================================================================
end module danu_iface
! =============================================================================


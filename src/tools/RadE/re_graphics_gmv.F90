#include "f90_assert.fpp"

module re_graphics_gmv

  use kinds, only: r8
  use re_encl_type
  implicit none
  private
  
  public :: gmv_open, gmv_close, gmv_write_encl
  public :: gmv_begin_variables, gmv_end_variables, gmv_write_face_var

  integer, parameter :: CELLDATA = 0, NODEDATA = 1

  !! Interfaces to the external Fortran-compatible C-procedures in f_gmvwrite.c
  interface
    subroutine gmvwrite_openfile_ir_f (filenam, isize, rsize)
      character(len=*) :: filenam
      integer :: isize, rsize
      end subroutine
    subroutine gmvwrite_closefile_f ()
      end subroutine
    subroutine gmvwrite_node_data_f (nndes, x, y, z)
      integer :: nndes
      real :: x(*), y(*), z(*)
      end subroutine
    subroutine gmvwrite_cell_header_f (ncells)
      integer :: ncells
      end subroutine
    subroutine gmvwrite_cell_type_f (cell_type, nverts, nodes)
      character(len=*) :: cell_type
      integer :: nverts, nodes(*)
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
      real :: ptime
      end subroutine gmvwrite_probtime_f
    subroutine gmvwrite_cycleno_f (cyclenum)
      integer :: cyclenum
      end subroutine gmvwrite_cycleno_f
    subroutine gmvwrite_variable_header_f ()
      end subroutine gmvwrite_variable_header_f
    subroutine gmvwrite_variable_name_data_f (data_type, varname, vids)
      integer :: data_type
      character(len=*) :: varname
      real :: vids(*)
      end subroutine
    subroutine gmvwrite_variable_endvars_f ()
      end subroutine
  end interface

contains

  subroutine gmv_open (file)
    character(len=*) :: file
    !! 4-byte integer data and 8-byte real data.
    call gmvwrite_openfile_ir_f (file, 4, 4)
  end subroutine gmv_open

  subroutine gmv_close ()
    call gmvwrite_closefile_f ()
  end subroutine gmv_close
  
  subroutine gmv_write_encl (e, full)
  
    type(encl), intent(in) :: e
    logical, intent(in), optional :: full
    
    logical :: call_sym
    
    !! Don't write the fully-developed enclosure surface unless the caller
    !! requests it and the enclosure actually includes symmetry attributes.
    call_sym = .false.
    if (present(full)) call_sym = full
    call_sym = call_sym .and. (any(e%mirror).or.(e%rot_axis>0))
    
    if (call_sym) then
      call gmv_write_encl_sym (e)
    else
      call gmv_write_encl_nosym (e)
    end if
    
  end subroutine gmv_write_encl
  
  
  subroutine gmv_write_encl_nosym (e)
  
    use string_utilities, only: i_to_c

    type(encl), intent(in) :: e
    
    integer :: j
    integer, pointer :: list(:)
    real, dimension(e%nnode) :: x, y, z
    
    !! Write the node coordinate data.
    x = real(e%x(1,:))
    y = real(e%x(2,:))
    z = real(e%x(3,:))
    call gmvwrite_node_data_f (e%nnode, x, y, z)
    
    !! Write the cell data.
    call gmvwrite_cell_header_f (e%nface)
    do j = 1, e%nface
      list => e%fnode(e%xface(j):e%xface(j+1)-1)
      select case (size(list))
      case (3)
        call gmvwrite_cell_type_f ('tri', 3, list)
      case (4)
        call gmvwrite_cell_type_f ('quad', 4, list)
      case default
        INSIST( .false. )
      end select
    end do
    
    !! Write the group IDs as  the cell meterial.
    call gmvwrite_material_header_f (size(e%group_id_list), CELLDATA)
    do j = 1, size(e%group_id_list)
      call gmvwrite_material_name_f ('Group'//i_to_c(e%group_id_list(j)))
    end do
    call gmvwrite_material_ids_f (e%gnum, CELLDATA)
    
  end subroutine gmv_write_encl_nosym
  
  subroutine gmv_write_encl_sym (e)
  
    use string_utilities, only: i_to_c

    type(encl), intent(in) :: e
    
    integer :: ncopy, i, j, k, n, nnode, nface
    integer, pointer :: gnum(:), gen(:), copy(:), fnode(:,:)
    character(len=8), pointer :: gen_name(:), copy_name(:)
    real, pointer :: x(:), y(:), z(:), xc(:), yc(:), zc(:)
    real(r8) :: theta, ct, st
    real(r8), parameter :: TWOPI= 6.2831853071795864769_r8
    
    !! Allocate coordinate arrays
    ncopy = 2**count(e%mirror)
    if (e%rot_axis > 0) ncopy = ncopy * e%num_rot
    n = ncopy * e%nnode
    allocate(x(n), y(n), z(n))
    
    !! Generating coordinates go in first.
    x(:e%nnode) = real(e%x(1,:))
    y(:e%nnode) = real(e%x(2,:))
    z(:e%nnode) = real(e%x(3,:))
    
    !! Generate the mirror symmetry copies.
    nnode = e%nnode
    do k = 1, 3
      if (e%mirror(k)) then
        select case (k)
        case (1)
          do i = 1, nnode
            x(nnode+i) = -x(i)
            y(nnode+i) =  y(i)
            z(nnode+i) =  z(i)
          end do
        case (2)
          do i = 1, nnode
            x(nnode+i) =  x(i)
            y(nnode+i) = -y(i)
            z(nnode+i) =  z(i)
          end do
        case (3)
          do i = 1, nnode
            x(nnode+i) =  x(i)
            y(nnode+i) =  y(i)
            z(nnode+i) = -z(i)
          end do
        end select
        nnode = 2*nnode
      end if
    end do
    
    !! Generate the rotation symmetry copies.
    if (e%rot_axis /= 0) then
      theta = TWOPI/e%num_rot
      do k = 1, e%num_rot-1
        ct = cos(k*theta)
        st = sin(k*theta)
        select case (e%rot_axis)
        case (1)
          do i = 1, nnode
            x(k*nnode+i) = x(i)
            y(k*nnode+i) = ct * y(i) - st * z(i)
            z(k*nnode+i) = st * y(i) + ct * z(i)
          end do
        case (2)
          do i = 1, nnode
            y(k*nnode+i) = y(i)
            z(k*nnode+i) = ct * z(i) - st * x(i)
            y(k*nnode+i) = st * z(i) + ct * x(i)
          end do
        case (3)
          do i = 1, nnode
            z(k*nnode+i) = z(i)
            x(k*nnode+i) = ct * x(i) - st * y(i)
            y(k*nnode+i) = st * x(i) + ct * y(i)
          end do
        end select
      end do
    end if
    
    call gmvwrite_node_data_f (size(x), x, y, z)
    
    call generated_faces (fnode, gnum, gen, gen_name, copy, copy_name)
    n = size(fnode,dim=2)
    allocate(xc(n), yc(n), zc(n))
    call gmvwrite_cell_header_f (size(fnode,dim=2))
    do j = 1, size(fnode,dim=2)
      if (fnode(4,j) == 0) then
        call gmvwrite_cell_type_f ('tri', 3, fnode(:,j))
        xc(j) = sum(x(fnode(:3,j))) / 3
        yc(j) = sum(y(fnode(:3,j))) / 3
        zc(j) = sum(z(fnode(:3,j))) / 3
      else
        call gmvwrite_cell_type_f ('quad', 4, fnode(:,j))
        xc(j) = sum(x(fnode(:4,j))) / 4
        yc(j) = sum(y(fnode(:4,j))) / 4
        zc(j) = sum(z(fnode(:4,j))) / 4
      end if
    end do
    deallocate(x, y, z)
    deallocate(fnode)
    
    !! Write the face centroid coordinates as cell data.
    call gmvwrite_variable_header_f ()
    call gmvwrite_variable_name_data_f (CELLDATA, 'Xc', xc)
    call gmvwrite_variable_name_data_f (CELLDATA, 'Yc', yc)
    call gmvwrite_variable_name_data_f (CELLDATA, 'Zc', zc)
    call gmvwrite_variable_endvars_f ()
    deallocate(xc, yc, zc)
    
    !! Write the group IDs as  the cell meterial.
    call gmvwrite_material_header_f (size(e%group_id_list), CELLDATA)
    do j = 1, size(e%group_id_list)
      call gmvwrite_material_name_f ('Group'//i_to_c(e%group_id_list(j)))
    end do
    call gmvwrite_material_ids_f (gnum, CELLDATA)
    deallocate(gnum)
    
    !! Face generation flag
    call gmvwrite_flag_header_f ()
    call gmvwrite_flag_name_f ('Generate', size(gen_name), CELLDATA)
    do j = 1, size(gen_name)
      call gmvwrite_flag_subname_f(gen_name(j))
    end do
    call gmvwrite_flag_data_f (CELLDATA, gen)
    deallocate(gen_name, gen)
    
    !! Face copy flag
    call gmvwrite_flag_name_f ('Copy', size(copy_name), CELLDATA)
    do j = 1, size(copy_name)
      call gmvwrite_flag_subname_f(copy_name(j))
    end do
    call gmvwrite_flag_data_f (CELLDATA, copy)
    deallocate(copy_name, copy)
    call gmvwrite_flag_endflag_f ()
    
  contains
  
    subroutine generated_faces (fnode, gnum, gen, gen_name, copy, copy_name)
    
      integer, pointer :: fnode(:,:), gnum(:), gen(:), copy(:)
      character(len=*), pointer :: gen_name(:), copy_name(:)
      
      integer :: n, nnode, i, j, k, ncopy
      integer, pointer :: list(:)
      
      !! The number of mesh copies from applying the symmetries.
      ncopy = 2**count(e%mirror)
      if (e%rot_axis > 0) ncopy = ncopy * e%num_rot

      !! The faces in the generated mesh
      allocate(fnode(4,ncopy*e%nface), gnum(ncopy*e%nface))
      allocate(gen_name(ncopy), gen(ncopy*e%nface))
      
      !! Copy the generating surface faces into the initial part of FNODE.
      gen_name(1) = 'Original'
      do j = 1, e%nface
        list => e%fnode(e%xface(j):e%xface(j+1)-1)
        select case (size(list))
        case (3)
          fnode(:3,j) = list
          fnode(4,j)  = 0
        case (4)
          fnode(:,j) = list
        case default
          INSIST(.false.)
        end select
        gen(j)  = 1
        gnum(j) = e%gnum(j)
      end do
      
      n = 1
      nface = e%nface
      nnode = e%nnode
      do k = 1, 3
        if (e%mirror(k)) then
          n = n + 1
          gen_name(n) = 'Mirror' // 'XYZ'(k:k)
          do i = 1, nface
            if (fnode(4,i) == 0) then ! tri cell
              fnode(1,nface+i) = fnode(1,i) + nnode
              fnode(2,nface+i) = fnode(3,i) + nnode
              fnode(3,nface+i) = fnode(2,i) + nnode
              fnode(4,nface+i) = 0
            else  ! quad cell
              fnode(1,nface+i) = fnode(1,i) + nnode
              fnode(2,nface+i) = fnode(4,i) + nnode
              fnode(3,nface+i) = fnode(3,i) + nnode
              fnode(4,nface+i) = fnode(2,i) + nnode
            end if
            gen(nface+i)  = n
            gnum(nface+i) = gnum(i)
          end do
          nface = 2*nface
          nnode = 2*nnode
        end if
      end do

      if (e%rot_axis > 0) then
        do k = 1, e%num_rot-1
          n = n + 1
          gen_name(n) = 'Rot' // 'XYZ'(e%rot_axis:e%rot_axis) // i_to_c(k)
          do i = 1, nface
            if (fnode(4,i) == 0) then ! tri cell
              fnode(:3,k*nface+i) = fnode(:3,i) + k*nnode
              fnode(4,k*nface+i) = 0
            else  ! quad cell
              fnode(:,k*nface+i) = fnode(:,i) + k*nnode
            end if
            gen(k*nface+i)  = n
            gnum(k*nface+i) = gnum(i)
          end do
        end do
      end if
      
      !! Copy selection flag data.
      allocate(copy_name(ncopy), copy(ncopy*e%nface))
      n = 0
      do i = 1, ncopy
        copy_name(i) = 'Copy' // i_to_c(i)
        do j = 1, e%nface
          n = n + 1
          copy(n) = i
        end do
      end do
      
    end subroutine generated_faces

  end subroutine gmv_write_encl_sym
  
  subroutine gmv_begin_variables (time, seq)
    real, intent(in), optional :: time
    integer, intent(in), optional :: seq
    if (present(time)) call gmvwrite_probtime_f (time)
    if (present(seq))  call gmvwrite_cycleno_f (seq)
    call gmvwrite_variable_header_f()
  end subroutine gmv_begin_variables

  subroutine gmv_end_variables ()
    call gmvwrite_variable_endvars_f()
  end subroutine gmv_end_variables

  subroutine gmv_write_face_var (e, u, name, full)
  
    type(encl),  intent(in) :: e
    real,             intent(in) :: u(:)
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: full
    
    logical :: call_sym
    
    ASSERT( size(u) == e%nface )
    
    call_sym = .false.
    if (present(full)) call_sym = full
    call_sym = call_sym .and. (any(e%mirror).or.(e%rot_axis>0))
    
    if (call_sym) then
      call gmv_write_face_var_sym (e, u, name)
    else
      call gmv_write_face_var_nosym (e, u, name)
    end if
    
  end subroutine gmv_write_face_var
    
  subroutine gmv_write_face_var_nosym (e, u, name)
  
    type(encl),  intent(in) :: e
    real,             intent(in) :: u(:)
    character(len=*), intent(in) :: name
    
    ASSERT( size(u) == e%nface )
    
    call gmvwrite_variable_name_data_f (CELLDATA, name, u)
    
  end subroutine gmv_write_face_var_nosym

  subroutine gmv_write_face_var_sym (e, u, name)
  
    type(encl),  intent(in) :: e
    real,             intent(in) :: u(:)
    character(len=*), intent(in) :: name
    
    integer :: i, k, ncopy, nface
    real, allocatable :: ufull(:)
    
    ASSERT( size(u) == e%nface )
    
    !! The number of mesh copies from applying the symmetries.
    ncopy = 2**count(e%mirror)
    if (e%rot_axis > 0) ncopy = ncopy * e%num_rot
    
    allocate(ufull(ncopy*e%nface))
    
    ufull(:e%nface) = u

    nface = e%nface
    do k = 1, 3
      if (e%mirror(k)) then
        do i = 1, nface
          ufull(nface+i) = ufull(i)
        end do
        nface = 2*nface
      end if
    end do

    if (e%rot_axis > 0) then
      do k = 1, e%num_rot-1
        do i = 1, nface
          ufull(k*nface+i) = ufull(i)
        end do
      end do
    end if
      
    call gmvwrite_variable_name_data_f (CELLDATA, name, ufull)
    
    deallocate(ufull)
    
  end subroutine gmv_write_face_var_sym

end module re_graphics_gmv

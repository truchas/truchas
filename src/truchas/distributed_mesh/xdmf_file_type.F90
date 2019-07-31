!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module xdmf_file_type

  use kinds, only: r8
  implicit none
  private

  type, public :: xdmf_file
    private
    integer :: nnode, nnode_tot, node_offset
    integer :: ncell, ncell_tot, cell_offset
    integer :: xmf, bin ! logical units
    integer :: pos = 0  ! current file position after parallel writes complete
    integer :: stage = 0
    character(:), allocatable :: binfile
  contains
    procedure :: open => open_file
    procedure :: close => close_file
    procedure :: write_mesh
    procedure :: begin_variables
    procedure :: end_variables
    generic   :: write_cell_var => write_cell_var_scalar, write_cell_var_vector
    generic   :: write_node_var => write_node_var_scalar, write_node_var_vector
    procedure, private :: write_cell_var_scalar
    procedure, private :: write_cell_var_vector
    procedure, private :: write_node_var_scalar
    procedure, private :: write_node_var_vector
  end type xdmf_file

contains

  subroutine open_file(this, filename)

    use parallel_communication, only: is_IOP
    use pgslib_module, only: pgslib_barrier

    class(xdmf_file), intent(out) :: this
    character(*), intent(in) :: filename

    !! The IO process handles the XDMF header file
    if (is_IOP) then
      open(newunit=this%xmf,file=filename//'.xmf',action='write',status='replace')
      !TODO: need iostat error handling
      write(this%xmf,'(a)') '<?xml version="1.0" ?>'
      write(this%xmf,'(a)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(this%xmf,'(a)') '<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
      write(this%xmf,'(t3,a)') '<Domain>'
    end if

    !! Everybody opens the binary data file.
    this%binfile = filename // '.bin'
    if (is_IOP) then  ! create an empty file
      open(newunit=this%bin,file=this%binfile,access='stream',action='write',status='replace')
      call pgslib_barrier
    else  ! wait for the IOP to create the file and then open it
      call pgslib_barrier
      open(newunit=this%bin,file=this%binfile,access='stream',action='write',status='old')
    end if
    !TODO: need iostat error handling

    !! Sanity check
    inquire(this%bin,pos=this%pos)
    INSIST(this%pos == 1)

    this%stage = 1  ! initialized, ready for mesh

  end subroutine open_file

  subroutine close_file(this)
    use parallel_communication, only: is_IOP
    class(xdmf_file), intent(inout) :: this
    INSIST(this%stage == 2 .or. this%stage == 3)
    if (is_IOP) then
      if (this%stage == 2) then ! no time snapshots were written
        write(this%xmf,'(t5,a)') '<Grid GridType="Uniform">'
        write(this%xmf,'(t7,a)') '<Geometry Reference="/Xdmf/Domain/Geometry"/>'
        write(this%xmf,'(t7,a)') '<Topology Reference="/Xdmf/Domain/Topology"/>'
      end if
      write(this%xmf,'(t5,a)') '</Grid>'
      write(this%xmf,'(t3,a)') '</Domain>'
      write(this%xmf,'(a)')    '</Xdmf>'
      close(this%xmf)
    end if
    close(this%bin)
    this%stage = 0
  end subroutine close_file

  subroutine write_mesh(this, mesh)

    use unstr_2d_mesh_type
    use parallel_communication, only: is_IOP, global_sum
    use,intrinsic :: iso_fortran_env, only: file_storage_size

    class(xdmf_file), intent(inout) :: this
    type(unstr_2d_mesh), intent(in) :: mesh

    integer :: j, n, nvrtx, nvrtx_tot, vrtx_offset, iolength, my_pos
    integer, allocatable :: connect(:)

    INSIST(this%stage == 1)

    this%nnode = mesh%nnode_onP               ! number of nodes I own
    this%nnode_tot = global_sum(this%nnode)   ! total number of nodes
    this%node_offset = my_offset(this%nnode)  ! my offset for writing node-based data

    this%ncell = mesh%ncell_onP               ! number of cells I own
    this%ncell_tot = global_sum(this%ncell)   ! total number of cells
    this%cell_offset = my_offset(this%ncell)  ! my offset for writing cell-based data

    nvrtx = mesh%cstart(mesh%ncell_onP+1) - 1 ! number of vertices I own
    nvrtx_tot = global_sum(nvrtx)             ! total number of vertices
    vrtx_offset = my_offset(nvrtx)            ! my offset for writing vertex-based data

    !! Write the XDMF node coordinate element
    if (is_IOP) then
      write(this%xmf,'(t5,a)') '<Geometry GeometryType="XY">'
      write(this%xmf,'(t7,a,i0,a,i0,a)') '<DataItem DataType="Float" Precision="8" Rank="2" &
          &Dimensions="', this%nnode_tot, ' 2" Format="Binary" Seek="', this%pos-1, '">'
      write(this%xmf,'(t9,a)') this%binfile
      write(this%xmf,'(t7,a)') '</DataItem>'
      write(this%xmf,'(t5,a)') '</Geometry>'
    end if

    !! Everybody writes their slab of the coordinate array.
    iolength = storage_size(mesh%x) / file_storage_size
    my_pos = this%pos + 2*this%node_offset*iolength
    write(this%bin,pos=my_pos) mesh%x(:,1:this%nnode)
    this%pos = this%pos + 2*this%nnode_tot*iolength

    !! Write the XDMF mesh connectivity element
    if (is_IOP) then
      write(this%xmf,'(t5,a,i0,a)') '<Topology TopologyType="Mixed" NumberOfElements="', this%ncell_tot, '">'
      write(this%xmf,'(t7,a,2(i0,a))') '<DataItem DataType="Int" Precision="4" Dimensions="', &
          this%ncell_tot + nvrtx_tot, '" Format="Binary" Seek="', this%pos-1, '">'
      write(this%xmf,'(t9,a)') this%binfile
      write(this%xmf,'(t7,a)') '</DataItem>'
      write(this%xmf,'(t5,a)') '</Topology>'
    end if

    !! Generate the XDMF-format connectivity array
    allocate(connect(this%ncell+nvrtx))
    n = 1
    do j = 1, this%ncell
      associate (list => mesh%cnode(mesh%cstart(j):mesh%cstart(j+1)-1))
        select case (size(list))
        case (3)
          connect(n) = 4  ! XDMF tri cell identifier
        case (4)
          connect(n) = 5  ! XDMF quad cell identifier
        case default
          INSIST(.false.)
        end select
        connect(n+1:n+size(list)) = mesh%node_ip%global_index(list) - 1  ! XDMF uses 0-based indexing
        n = n + size(list) + 1
      end associate
    end do
    ASSERT(n == size(connect)+1)

    !! Everybody writes their slab of the connectivity array.
    iolength = storage_size(connect)/file_storage_size
    my_pos = this%pos + (this%cell_offset + vrtx_offset)*iolength
    write(this%bin,pos=my_pos) connect
    this%pos = this%pos + (this%ncell_tot + nvrtx_tot)*iolength

    this%stage = 2  ! mesh written

  end subroutine write_mesh

  subroutine begin_variables(this, time)
    use parallel_communication, only: is_IOP
    class(xdmf_file), intent(inout) :: this
    real(r8), intent(in) :: time
    INSIST(this%stage == 2 .or. this%stage == 3)
    if (is_IOP) then
      if (this%stage == 2) &
          write(this%xmf,'(t5,a)') '<Grid GridType="Collection" CollectionType="Temporal">'
      write(this%xmf,'(t7,a)') '<Grid GridType="Uniform">'
      write(this%xmf,'(t9,a)') '<Geometry Reference="/Xdmf/Domain/Geometry"/>'
      write(this%xmf,'(t9,a)') '<Topology Reference="/Xdmf/Domain/Topology"/>'
      write(this%xmf,'(t9,a,es23.15,a)') '<Time Value="', time, '"/>'
    end if
    this%stage = 4 ! ! ready for writing snapshot variable
  end subroutine begin_variables

  subroutine end_variables(this)
    use parallel_communication, only: is_IOP
    class(xdmf_file), intent(inout) :: this
    INSIST(this%stage == 4)
    if (is_IOP) write(this%xmf,'(t7,a)') '</Grid>' ! matching tag written by begin_variables
    this%stage = 3  ! ready for next time snapshot
  end subroutine end_variables

  subroutine write_cell_var_scalar(this, u, name)

    use,intrinsic :: iso_fortran_env, only: file_storage_size
    use parallel_communication, only: is_IOP

    class(xdmf_file), intent(inout) :: this
    real(r8), intent(in) :: u(:)
    character(*), intent(in) :: name

    integer :: iolength, my_pos

    INSIST(this%stage == 4)
    ASSERT(size(u) >= this%ncell)

    !! Write the XDMF header for the variable.
    if (is_IOP) then
      write(this%xmf,'(t9,3a)') '<Attribute Name="', name, '" Center="Cell">'
      write(this%xmf,'(t11,a,i0,2(a,i0))') '<DataItem DataType="Float" Precision="8" Dimensions="', &
          this%ncell_tot, '" Format="Binary" Seek="', this%pos-1, '">'
      write(this%xmf,'(t13,a)') this%binfile
      write(this%xmf,'(t11,a)') '</DataItem>'
      write(this%xmf,'(t9,a)') '</Attribute>'
    end if

    !! Everybody writes their slab of the array.
    iolength = storage_size(u)/file_storage_size
    my_pos = this%pos + this%cell_offset*iolength
    write(this%bin,pos=my_pos) u(:this%ncell)
    this%pos = this%pos + this%ncell_tot*iolength

  end subroutine write_cell_var_scalar

  subroutine write_cell_var_vector(this, u, name)

    use,intrinsic :: iso_fortran_env, only: file_storage_size
    use parallel_communication, only: is_IOP

    class(xdmf_file), intent(inout) :: this
    real(r8), intent(in) :: u(:,:)
    character(*), intent(in) :: name

    integer :: iolength, my_pos, j

    INSIST(this%stage == 4)
    ASSERT(size(u,2) >= this%ncell)

    !! XDMF vector attributes are required to have 3 values per grid location
    !! even for a 2D mesh. We hack around it here by supplying a dummy 0 value
    !! for the z-component of the vector. Ideally the number of vector
    !! components would not be fixed at all.
    !! See https://gitlab.kitware.com/xdmf/xdmf/issues/17
    INSIST(size(u,1) == 2)

    !! Write the XDMF header for the variable.
    if (is_IOP) then
      write(this%xmf,'(t7,3a)') '<Attribute Name="', name, '" Center="Cell" AttributeType="Vector">'
      write(this%xmf,'(t9,a,2(i0,a))') '<DataItem DataType="Float" Precision="8" Rank="2" &
          &Dimensions="', this%ncell_tot, ' 3" Format="Binary" Seek="', this%pos-1, '">'
      write(this%xmf,'(t11,a)') this%binfile
      write(this%xmf,'(t9,a)') '</DataItem>'
      write(this%xmf,'(t7,a)') '</Attribute>'
    end if

    !! Everybody writes their slab of the array.
    iolength = storage_size(u)/file_storage_size
    my_pos = this%pos + 3*this%cell_offset*iolength
    write(this%bin,pos=my_pos) (u(:,j), 0.0_r8, j = 1, this%ncell)
    this%pos = this%pos + 3*this%ncell_tot*iolength

  end subroutine write_cell_var_vector

  subroutine write_node_var_scalar(this, u, name)

    use,intrinsic :: iso_fortran_env, only: file_storage_size
    use parallel_communication, only: is_IOP

    class(xdmf_file), intent(inout) :: this
    real(r8), intent(in) :: u(:)
    character(*), intent(in) :: name

    integer :: iolength, my_pos

    INSIST(this%stage == 4)
    ASSERT(size(u) >= this%nnode)

    !! Write the XDMF header for the variable.
    if (is_IOP) then
      write(this%xmf,'(t7,3a)') '<Attribute Name="', name, '" Center="Node">'
      write(this%xmf,'(t9,a,i0,2(a,i0))') '<DataItem DataType="Float" Precision="8" Dimensions="', &
          this%nnode_tot, '" Format="Binary" Seek="', this%pos-1, '">'
      write(this%xmf,'(t11,a)') this%binfile
      write(this%xmf,'(t9,a)') '</DataItem>'
      write(this%xmf,'(t7,a)') '</Attribute>'
    end if

    !! Everybody writes their slab of the array.
    iolength = storage_size(u)/file_storage_size
    my_pos = this%pos + this%node_offset*iolength
    write(this%bin,pos=my_pos) u(:this%nnode)
    this%pos = this%pos + this%nnode_tot*iolength

  end subroutine write_node_var_scalar

  subroutine write_node_var_vector(this, u, name)

    use,intrinsic :: iso_fortran_env, only: file_storage_size
    use parallel_communication, only: is_IOP

    class(xdmf_file), intent(inout) :: this
    real(r8), intent(in) :: u(:,:)
    character(*), intent(in) :: name

    integer :: iolength, my_pos, j

    INSIST(this%stage == 4)
    ASSERT(size(u,2) >= this%nnode)

    !! XDMF vector attributes are required to have 3 values per grid location
    !! even for a 2D mesh. We hack around it here by supplying a dummy 0 value
    !! for the z-component of the vector. Ideally the number of vector
    !! components would not be fixed at all.
    !! See https://gitlab.kitware.com/xdmf/xdmf/issues/17
    INSIST(size(u,1) == 2)

    !! Write the XDMF header for the variable.
    if (is_IOP) then
      write(this%xmf,'(t7,3a)') '<Attribute Name="', name, '" Center="Node" AttributeType="Vector">'
      write(this%xmf,'(t9,a,2(i0,a))') '<DataItem DataType="Float" Precision="8" Rank="2" &
          &Dimensions="', this%nnode_tot, ' 3" Format="Binary" Seek="', this%pos-1, '">'
      write(this%xmf,'(t11,a)') this%binfile
      write(this%xmf,'(t9,a)') '</DataItem>'
      write(this%xmf,'(t7,a)') '</Attribute>'
    end if

    !! Everybody writes their slab of the array.
    iolength = storage_size(u)/file_storage_size
    my_pos = this%pos + 3*this%node_offset*iolength
    write(this%bin,pos=my_pos) (u(:,j), 0.0_r8, j = 1, this%nnode)
    this%pos = this%pos + 3*this%nnode_tot*iolength

  end subroutine write_node_var_vector

  !! Quick-n-dirty exclusive prefix sum (done inefficiently!)
  integer function my_offset(n)
    use parallel_communication
    integer, intent(in) :: n
    integer :: buffer(nPE)
    call collate(buffer, n)
    call broadcast(buffer)
    my_offset = sum(buffer(1:this_PE-1))
  end function my_offset

end module xdmf_file_type

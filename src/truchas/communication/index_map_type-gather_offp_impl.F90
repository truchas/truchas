!! Implementation of INDEX_MAP GATHER_OFFP Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

submodule(index_map_type) gather_offp_impl
implicit none
contains

!!!! I4 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath1_i4_1(this, local_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: local_data(:)
    call gath2_i4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine gath2_i4_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(in) :: onp_data(:)
    integer(i4), intent(inout) :: offp_data(:)
    integer :: ierr
    integer(i4), allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    onp_buf = onp_data(this%onp_index)
    call MPI_Neighbor_alltoallv(onp_buf, this%onp_counts, this%onp_displs, MPI_INTEGER4, &
        offp_data, this%offp_counts, this%offp_displs, MPI_INTEGER4, this%gather_comm, ierr)
  end subroutine

  module subroutine gath1_i4_2(this, local_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: local_data(:,:)
    call gath2_i4_2(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_i4_2(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(in) :: onp_data(:,:)
    integer(i4), intent(inout) :: offp_data(:,:)
    integer :: ierr
    integer(i4), allocatable :: onp_buf(:,:)
    integer :: block_type
    if (.not.allocated(this%offp_index)) return
    onp_buf = onp_data(:,this%onp_index)
    call MPI_Type_contiguous(size(onp_data,dim=1), MPI_INTEGER4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Neighbor_alltoallv(onp_buf, this%onp_counts, this%onp_displs, block_type, &
        offp_data, this%offp_counts, this%offp_displs, block_type, this%gather_comm, ierr)
    call MPI_Type_free(block_type, ierr)
  end subroutine

  module subroutine gath1_i4_3(this, local_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: local_data(:,:,:)
    call gath2_i4_3(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_i4_3(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(in) :: onp_data(:,:,:)
    integer(i4), intent(inout) :: offp_data(:,:,:)
    integer :: ierr
    integer(i4), allocatable :: onp_buf(:,:,:)
    integer :: block_type
    if (.not.allocated(this%offp_index)) return
    onp_buf = onp_data(:,:,this%onp_index)
    call MPI_Type_contiguous(size(onp_data(:,:,1)), MPI_INTEGER4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Neighbor_alltoallv(onp_buf, this%onp_counts, this%onp_displs, block_type, &
        offp_data, this%offp_counts, this%offp_displs, block_type, this%gather_comm, ierr)
    call MPI_Type_free(block_type, ierr)
  end subroutine

!!!! R4 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath1_r4_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: local_data(:)
    call gath2_r4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine gath2_r4_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(r4), intent(in) :: onp_data(:)
    real(r4), intent(inout) :: offp_data(:)
    integer :: ierr
    real(r4), allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    onp_buf = onp_data(this%onp_index)
    call MPI_Neighbor_alltoallv(onp_buf, this%onp_counts, this%onp_displs, MPI_REAL4, &
        offp_data, this%offp_counts, this%offp_displs, MPI_REAL4, this%gather_comm, ierr)
  end subroutine

  module subroutine gath1_r4_2(this, local_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: local_data(:,:)
    call gath2_r4_2(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_r4_2(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(r4), intent(in) :: onp_data(:,:)
    real(r4), intent(inout) :: offp_data(:,:)
    integer :: ierr
    real(r4), allocatable :: onp_buf(:,:)
    integer :: block_type
    if (.not.allocated(this%offp_index)) return
    onp_buf = onp_data(:,this%onp_index)
    call MPI_Type_contiguous(size(onp_data,dim=1), MPI_REAL4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Neighbor_alltoallv(onp_buf, this%onp_counts, this%onp_displs, block_type, &
        offp_data, this%offp_counts, this%offp_displs, block_type, this%gather_comm, ierr)
    call MPI_Type_free(block_type, ierr)
  end subroutine

  module subroutine gath1_r4_3(this, local_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: local_data(:,:,:)
    call gath2_r4_3(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_r4_3(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(r4), intent(in) :: onp_data(:,:,:)
    real(r4), intent(inout) :: offp_data(:,:,:)
    integer :: ierr
    real(r4), allocatable :: onp_buf(:,:,:)
    integer :: block_type
    if (.not.allocated(this%offp_index)) return
    onp_buf = onp_data(:,:,this%onp_index)
    call MPI_Type_contiguous(size(onp_data(:,:,1)), MPI_REAL4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Neighbor_alltoallv(onp_buf, this%onp_counts, this%onp_displs, block_type, &
        offp_data, this%offp_counts, this%offp_displs, block_type, this%gather_comm, ierr)
    call MPI_Type_free(block_type, ierr)
  end subroutine

!!!! R8 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath1_r8_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: local_data(:)
    call gath2_r8_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine gath2_r8_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(r8), intent(in) :: onp_data(:)
    real(r8), intent(inout) :: offp_data(:)
    integer :: ierr
    real(r8), allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    onp_buf = onp_data(this%onp_index)
    call MPI_Neighbor_alltoallv(onp_buf, this%onp_counts, this%onp_displs, MPI_REAL8, &
        offp_data, this%offp_counts, this%offp_displs, MPI_REAL8, this%gather_comm, ierr)
  end subroutine

  module subroutine gath1_r8_2(this, local_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: local_data(:,:)
    call gath2_r8_2(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_r8_2(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(r8), intent(in) :: onp_data(:,:)
    real(r8), intent(inout) :: offp_data(:,:)
    integer :: ierr
    real(r8), allocatable :: onp_buf(:,:)
    integer :: block_type
    if (.not.allocated(this%offp_index)) return
    onp_buf = onp_data(:,this%onp_index)
    call MPI_Type_contiguous(size(onp_data,dim=1), MPI_REAL8, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Neighbor_alltoallv(onp_buf, this%onp_counts, this%onp_displs, block_type, &
        offp_data, this%offp_counts, this%offp_displs, block_type, this%gather_comm, ierr)
    call MPI_Type_free(block_type, ierr)
  end subroutine

  module subroutine gath1_r8_3(this, local_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: local_data(:,:,:)
    call gath2_r8_3(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_r8_3(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(r8), intent(in) :: onp_data(:,:,:)
    real(r8), intent(inout) :: offp_data(:,:,:)
    integer :: ierr
    real(r8), allocatable :: onp_buf(:,:,:)
    integer :: block_type
    if (.not.allocated(this%offp_index)) return
    onp_buf = onp_data(:,:,this%onp_index)
    call MPI_Type_contiguous(size(onp_data(:,:,1)), MPI_REAL8, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Neighbor_alltoallv(onp_buf, this%onp_counts, this%onp_displs, block_type, &
        offp_data, this%offp_counts, this%offp_displs, block_type, this%gather_comm, ierr)
    call MPI_Type_free(block_type, ierr)
  end subroutine

!!!! LOGICAL DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath1_dl_1(this, local_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: local_data(:)
    call gath2_dl_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine gath2_dl_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    logical, intent(in) :: onp_data(:)
    logical, intent(inout) :: offp_data(:)
    integer :: ierr
    logical, allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    onp_buf = onp_data(this%onp_index)
    call MPI_Neighbor_alltoallv(onp_buf, this%onp_counts, this%onp_displs, MPI_LOGICAL, &
        offp_data, this%offp_counts, this%offp_displs, MPI_LOGICAL, this%gather_comm, ierr)
  end subroutine

  module subroutine gath1_dl_2(this, local_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: local_data(:,:)
    call gath2_dl_2(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_dl_2(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    logical, intent(in) :: onp_data(:,:)
    logical, intent(inout) :: offp_data(:,:)
    integer :: ierr
    logical, allocatable :: onp_buf(:,:)
    integer :: block_type
    if (.not.allocated(this%offp_index)) return
    onp_buf = onp_data(:,this%onp_index)
    call MPI_Type_contiguous(size(onp_data,dim=1), MPI_LOGICAL, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Neighbor_alltoallv(onp_buf, this%onp_counts, this%onp_displs, block_type, &
        offp_data, this%offp_counts, this%offp_displs, block_type, this%gather_comm, ierr)
    call MPI_Type_free(block_type, ierr)
  end subroutine

  module subroutine gath1_dl_3(this, local_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: local_data(:,:,:)
    call gath2_dl_3(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_dl_3(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    logical, intent(in) :: onp_data(:,:,:)
    logical, intent(inout) :: offp_data(:,:,:)
    integer :: ierr
    logical, allocatable :: onp_buf(:,:,:)
    integer :: block_type
    if (.not.allocated(this%offp_index)) return
    onp_buf = onp_data(:,:,this%onp_index)
    call MPI_Type_contiguous(size(onp_data(:,:,1)), MPI_LOGICAL, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Neighbor_alltoallv(onp_buf, this%onp_counts, this%onp_displs, block_type, &
        offp_data, this%offp_counts, this%offp_displs, block_type, this%gather_comm, ierr)
    call MPI_Type_free(block_type, ierr)
  end subroutine

end submodule gather_offp_impl

!! Implementation of INDEX_MAP SCATTER_OFFP Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

submodule(index_map_type) scatter_offp_impl
implicit none
contains

  module subroutine scat1_sum_i4_1(this, local_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: local_data(:)
    call scat2_sum_i4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_sum_i4_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: onp_data(:)
    integer(i4), intent(in) :: offp_data(:)
    integer :: j, ierr
    integer(i4), allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    allocate(onp_buf(size(this%onp_index)))
    call MPI_Neighbor_alltoallv(offp_data, this%offp_counts, this%offp_displs, MPI_INTEGER4, &
        onp_buf, this%onp_counts, this%onp_displs, MPI_INTEGER4, this%scatter_comm, ierr)
    do j = 1, size(onp_buf) !NB: must be sequential; ONP_INDEX is many-to-one
      onp_data(this%onp_index(j)) = onp_data(this%onp_index(j)) + onp_buf(j)
    end do
  end subroutine

  module subroutine scat1_sum_r4_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: local_data(:)
    call scat2_sum_r4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_sum_r4_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: onp_data(:)
    real(r4), intent(in) :: offp_data(:)
    integer :: j, ierr
    real(r4), allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    allocate(onp_buf(size(this%onp_index)))
    call MPI_Neighbor_alltoallv(offp_data, this%offp_counts, this%offp_displs, MPI_REAL4, &
        onp_buf, this%onp_counts, this%onp_displs, MPI_REAL4, this%scatter_comm, ierr)
    do j = 1, size(onp_buf) !NB: must be sequential; ONP_INDEX is many-to-one
      onp_data(this%onp_index(j)) = onp_data(this%onp_index(j)) + onp_buf(j)
    end do
  end subroutine

  module subroutine scat1_sum_r8_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: local_data(:)
    call scat2_sum_r8_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_sum_r8_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: onp_data(:)
    real(r8), intent(in) :: offp_data(:)
    integer :: j, ierr
    real(r8), allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    allocate(onp_buf(size(this%onp_index)))
    call MPI_Neighbor_alltoallv(offp_data, this%offp_counts, this%offp_displs, MPI_REAL8, &
        onp_buf, this%onp_counts, this%onp_displs, MPI_REAL8, this%scatter_comm, ierr)
    do j = 1, size(onp_buf) !NB: must be sequential; ONP_INDEX is many-to-one
      onp_data(this%onp_index(j)) = onp_data(this%onp_index(j)) + onp_buf(j)
    end do
  end subroutine

  module subroutine scat1_min_i4_1(this, local_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: local_data(:)
    call scat2_min_i4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_min_i4_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: onp_data(:)
    integer(i4), intent(in) :: offp_data(:)
    integer :: j, ierr
    integer(i4), allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    allocate(onp_buf(size(this%onp_index)))
    call MPI_Neighbor_alltoallv(offp_data, this%offp_counts, this%offp_displs, MPI_INTEGER4, &
        onp_buf, this%onp_counts, this%onp_displs, MPI_INTEGER4, this%scatter_comm, ierr)
    do j = 1, size(onp_buf) !NB: must be sequential; ONP_INDEX is many-to-one
      onp_data(this%onp_index(j)) = min(onp_data(this%onp_index(j)), onp_buf(j))
    end do
  end subroutine

  module subroutine scat1_min_r4_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: local_data(:)
    call scat2_min_r4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_min_r4_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: onp_data(:)
    real(r4), intent(in) :: offp_data(:)
    integer :: j, ierr
    real(r4), allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    allocate(onp_buf(size(this%onp_index)))
    call MPI_Neighbor_alltoallv(offp_data, this%offp_counts, this%offp_displs, MPI_REAL4, &
        onp_buf, this%onp_counts, this%onp_displs, MPI_REAL4, this%scatter_comm, ierr)
    do j = 1, size(onp_buf) !NB: must be sequential; ONP_INDEX is many-to-one
      onp_data(this%onp_index(j)) = min(onp_data(this%onp_index(j)), onp_buf(j))
    end do
  end subroutine

  module subroutine scat1_min_r8_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: local_data(:)
    call scat2_min_r8_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_min_r8_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: onp_data(:)
    real(r8), intent(in) :: offp_data(:)
    integer :: j, ierr
    real(r8), allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    allocate(onp_buf(size(this%onp_index)))
    call MPI_Neighbor_alltoallv(offp_data, this%offp_counts, this%offp_displs, MPI_REAL8, &
        onp_buf, this%onp_counts, this%onp_displs, MPI_REAL8, this%scatter_comm, ierr)
    do j = 1, size(onp_buf) !NB: must be sequential; ONP_INDEX is many-to-one
      onp_data(this%onp_index(j)) = min(onp_data(this%onp_index(j)), onp_buf(j))
    end do
  end subroutine

  module subroutine scat1_max_i4_1(this, local_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: local_data(:)
    call scat2_max_i4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_max_i4_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: onp_data(:)
    integer(i4), intent(in) :: offp_data(:)
    integer :: j, ierr
    integer(i4), allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    allocate(onp_buf(size(this%onp_index)))
    call MPI_Neighbor_alltoallv(offp_data, this%offp_counts, this%offp_displs, MPI_INTEGER4, &
        onp_buf, this%onp_counts, this%onp_displs, MPI_INTEGER4, this%scatter_comm, ierr)
    do j = 1, size(onp_buf) !NB: must be sequential; ONP_INDEX is many-to-one
      onp_data(this%onp_index(j)) = max(onp_data(this%onp_index(j)), onp_buf(j))
    end do
  end subroutine

  module subroutine scat1_max_r4_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: local_data(:)
    call scat2_max_r4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_max_r4_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: onp_data(:)
    real(r4), intent(in) :: offp_data(:)
    integer :: j, ierr
    real(r4), allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    allocate(onp_buf(size(this%onp_index)))
    call MPI_Neighbor_alltoallv(offp_data, this%offp_counts, this%offp_displs, MPI_REAL4, &
        onp_buf, this%onp_counts, this%onp_displs, MPI_REAL4, this%scatter_comm, ierr)
    do j = 1, size(onp_buf) !NB: must be sequential; ONP_INDEX is many-to-one
      onp_data(this%onp_index(j)) = max(onp_data(this%onp_index(j)), onp_buf(j))
    end do
  end subroutine

  module subroutine scat1_max_r8_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: local_data(:)
    call scat2_max_r8_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_max_r8_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: onp_data(:)
    real(r8), intent(in) :: offp_data(:)
    integer :: j, ierr
    real(r8), allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    allocate(onp_buf(size(this%onp_index)))
    call MPI_Neighbor_alltoallv(offp_data, this%offp_counts, this%offp_displs, MPI_REAL8, &
        onp_buf, this%onp_counts, this%onp_displs, MPI_REAL8, this%scatter_comm, ierr)
    do j = 1, size(onp_buf) !NB: must be sequential; ONP_INDEX is many-to-one
      onp_data(this%onp_index(j)) = max(onp_data(this%onp_index(j)), onp_buf(j))
    end do
  end subroutine

  module subroutine scat1_or_dl_1(this, local_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: local_data(:)
    call scat2_or_dl_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_or_dl_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: onp_data(:)
    logical, intent(in) :: offp_data(:)
    integer :: j, ierr
    logical, allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    allocate(onp_buf(size(this%onp_index)))
    call MPI_Neighbor_alltoallv(offp_data, this%offp_counts, this%offp_displs, MPI_LOGICAL, &
        onp_buf, this%onp_counts, this%onp_displs, MPI_LOGICAL, this%scatter_comm, ierr)
    do j = 1, size(onp_buf) !NB: must be sequential; ONP_INDEX is many-to-one
      onp_data(this%onp_index(j)) = onp_data(this%onp_index(j)) .or. onp_buf(j)
    end do
  end subroutine

  module subroutine scat1_and_dl_1(this, local_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: local_data(:)
    call scat2_and_dl_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_and_dl_1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: onp_data(:)
    logical, intent(in) :: offp_data(:)
    integer :: j, ierr
    logical, allocatable :: onp_buf(:)
    if (.not.allocated(this%offp_index)) return
    allocate(onp_buf(size(this%onp_index)))
    call MPI_Neighbor_alltoallv(offp_data, this%offp_counts, this%offp_displs, MPI_LOGICAL, &
        onp_buf, this%onp_counts, this%onp_displs, MPI_LOGICAL, this%scatter_comm, ierr)
    do j = 1, size(onp_buf) !NB: must be sequential; ONP_INDEX is many-to-one
      onp_data(this%onp_index(j)) = onp_data(this%onp_index(j)) .and. onp_buf(j)
    end do
  end subroutine

end submodule scatter_offp_impl

#include "f90_assert.fpp"

module parallel_csr_matrix

  use kinds, only: r8
  use index_partitioning
  use GraphModule, only: NGraphType
  implicit none
  private
  
  public :: pcsr_graph_create, pcsr_graph_delete
  public :: pcsr_graph_add_clique, pcsr_graph_add_edge, pcsr_graph_fill_complete
  public :: pcsr_graph_filled
  
  public :: pcsr_matrix_create, pcsr_matrix_delete
  public :: pcsr_matrix_set_all_values, pcsr_matrix_set_value, pcsr_matrix_add_to_value
  public :: pcsr_matrix_project_out_index
  public :: pcsr_matrix_get_row_view, pcsr_matrix_get_diagonal_copy
  public :: pcsr_matrix_graph

  type, public :: pcsr_graph
    type(ip_desc), pointer :: row_ip => null()
    integer, allocatable :: xadj(:), adjncy(:)
    type(NGraphType), pointer :: g => null() ! temporary used during construction
  end type pcsr_graph
  
  type, public :: pcsr_matrix
    real(r8), pointer :: data(:) => null()
    type(pcsr_graph), pointer :: graph => null()
    logical :: graph_dealloc = .false.
    integer :: nrow = 0, nrow_onP = 0
  end type pcsr_matrix
  
  interface pcsr_matrix_create
    module procedure pcsr_matrix_create, pcsr_matrix_create_mold
  end interface
  
  interface pcsr_graph_add_edge
    module procedure pcsr_graph_add_edge, pcsr_graph_add_edge1
  end interface
  
contains

  subroutine pcsr_graph_create (this, row_ip)
    use GraphModule, only: CreateGraph
    type(pcsr_graph), intent(out) :: this
    type(ip_desc), pointer :: row_ip
    ASSERT(associated(row_ip))
    this%row_ip => row_ip
    allocate(this%g)
    this%g = CreateGraph(local_size(row_ip), self_edge=.true.)
  end subroutine pcsr_graph_create
  
  subroutine pcsr_graph_add_clique (this, indices)
    use GraphModule, only: AddClique
    type(pcsr_graph), intent(inout) :: this
    integer, intent(in) :: indices(:)
    ASSERT(associated(this%g))
    call AddClique (this%g, indices)
  end subroutine pcsr_graph_add_clique
  
  subroutine pcsr_graph_add_edge (this, n, nbr)
    use GraphModule, only: AddEdge
    type(pcsr_graph), intent(inout) :: this
    integer, intent(in) :: n, nbr
    ASSERT(associated(this%g))
    call AddEdge (this%g, n, nbr)
  end subroutine pcsr_graph_add_edge
  
  subroutine pcsr_graph_add_edge1 (this, n, nbr)
    use GraphModule, only: AddEdge
    type(pcsr_graph), intent(inout) :: this
    integer, intent(in) :: n, nbr(:)
    ASSERT(associated(this%g))
    call AddEdge (this%g, n, nbr)
  end subroutine pcsr_graph_add_edge1
  
  subroutine pcsr_graph_fill_complete (this)
    use GraphModule, only: GetNeighborStructure, DeleteGraph
    type(pcsr_graph), intent(inout) :: this
    ASSERT(associated(this%g))
    call GetNeighborStructure (this%g, this%xadj, this%adjncy)
    call DeleteGraph (this%g)
    deallocate(this%g)
  end subroutine pcsr_graph_fill_complete
  
  subroutine pcsr_graph_delete (this)
    use GraphModule, only: DeleteGraph
    type(pcsr_graph), intent(inout) :: this
    this%row_ip => null()
    if (associated(this%g)) then
      call DeleteGraph (this%g)
      deallocate(this%g)
    end if
  end subroutine pcsr_graph_delete
  
  logical function pcsr_graph_filled (this)
    type(pcsr_graph), intent(in) :: this
    pcsr_graph_filled = allocated(this%xadj)
  end function pcsr_graph_filled

  subroutine pcsr_matrix_create (this, graph, take_graph)
    type(pcsr_matrix), intent(out) :: this
    type(pcsr_graph), pointer :: graph
    logical, intent(in) :: take_graph
    ASSERT(associated(graph))
    ASSERT(pcsr_graph_filled(graph))
    this%graph => graph
    this%graph_dealloc = take_graph
    allocate(this%data(size(graph%adjncy)))
    this%nrow = local_size(graph%row_ip)
    this%nrow_onP = onP_size(graph%row_ip)
  end subroutine pcsr_matrix_create
  
  subroutine pcsr_matrix_create_mold (this, mold)
    type(pcsr_matrix), intent(out) :: this
    type(pcsr_matrix), intent(in)  :: mold
    call pcsr_matrix_create (this, mold%graph, take_graph=.false.)
  end subroutine pcsr_matrix_create_mold
  
  subroutine pcsr_matrix_delete (this)
    type(pcsr_matrix), intent(inout) :: this
    type(pcsr_matrix) :: default
    if (associated(this%data)) deallocate(this%data)
    if (this%graph_dealloc) then
      if (associated(this%graph)) then
        call pcsr_graph_delete (this%graph)
        deallocate(this%graph)
      end if
    end if
    this = default  ! assign default initialization values
  end subroutine pcsr_matrix_delete
  
  function pcsr_matrix_graph (this) result (graph)
    type(pcsr_matrix), intent(in) :: this
    type(pcsr_graph), pointer :: graph
    graph => this%graph
  end function pcsr_matrix_graph
  
  subroutine pcsr_matrix_set_all_values (this, value)
    type(pcsr_matrix), intent(inout) :: this
    real(r8), intent(in) :: value
    ASSERT(associated(this%data))
    this%data = value
  end subroutine pcsr_matrix_set_all_values
  
  subroutine pcsr_matrix_add_to_value (this, row, col, value)
    type(pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value
    integer :: n
    ASSERT(row >= 1 .and. row <= this%nrow)
    n = array_index(this%graph, row, col)
    ASSERT(n /= 0)
    this%data(n) = this%data(n) + value
  end subroutine pcsr_matrix_add_to_value
  
  subroutine pcsr_matrix_set_value (this, row, col, value)
    type(pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: row, col
    real(r8), intent(in) :: value
    integer :: n
    ASSERT(row >= 1 .and. row <= this%nrow)
    n = array_index(this%graph, row, col)
    ASSERT(n /= 0)
    this%data(n) = value
  end subroutine pcsr_matrix_set_value

  pure integer function array_index (graph, row, col) result (k)
    type(pcsr_graph), intent(in) :: graph
    integer, intent(in) :: row, col
    do k = graph%xadj(row), graph%xadj(row+1)-1
      if (graph%adjncy(k) == col) return
    end do
    k = 0
  end function array_index
    
  subroutine pcsr_matrix_project_out_index (this, index)
    type(pcsr_matrix), intent(inout) :: this
    integer, intent(in) :: index
    integer :: m, n, lmn, lnm
    m = index
    do lmn = this%graph%xadj(m), this%graph%xadj(m+1)-1
      this%data(lmn) = 0.0_r8
      n = this%graph%adjncy(lmn)
      if (n == m) cycle ! diagonal element
      do lnm = this%graph%xadj(n), this%graph%xadj(n+1)-1
        if (this%graph%adjncy(lnm) == m) exit
      end do
      ASSERT(lnm < this%graph%xadj(n+1))
      this%data(lnm) = 0.0_r8
    end do
  end subroutine pcsr_matrix_project_out_index
  
  subroutine pcsr_matrix_get_row_view (this, row, values, indices)
    type(pcsr_matrix), intent(in) :: this
    integer, intent(in) :: row
    real(r8), pointer :: values(:)
    integer, pointer :: indices(:)
    ASSERT(row >= 1 .and. row <= this%nrow)
    values => this%data(this%graph%xadj(row):this%graph%xadj(row+1)-1)
    indices => this%graph%adjncy(this%graph%xadj(row):this%graph%xadj(row+1)-1)
  end subroutine pcsr_matrix_get_row_view
  
  subroutine pcsr_matrix_get_diagonal_copy (this, d)
    type(pcsr_matrix), intent(in) :: this
    real(r8), intent(inout) :: d(:)
    integer :: i, k
    ASSERT(size(d) >= this%nrow_onP)
    do i = 1, this%nrow_onP
      do k = this%graph%xadj(i), this%graph%xadj(i+1)-1
        if (this%graph%adjncy(k) == i) then
          d(i) = this%data(k)
          exit
        end if
      end do
    end do
  end subroutine pcsr_matrix_get_diagonal_copy

end module parallel_csr_matrix

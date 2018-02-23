module flow_projection_type

  use kinds, only: r8
  use truchas_logging_services
  use truchas_timers
  use flow_mesh_type
  use unstr_mesh_type
  use index_partitioning
  use hypre_hybrid
  use pcsr_matrix_type
  use parameter_list_type
  implicit none
  private

  public :: flow_projection

  type :: flow_projection
    type(flow_mesh), pointer :: mesh ! unowned reference
    type(hypre_hybrid) :: solver
    type(parameter_list) :: p
    real(r8), allocatable :: rhs(:)
  contains
    procedure :: read_params
    procedure :: init
    procedure :: setup
    procedure :: solve
  end type flow_projection

contains

  subroutine read_params(this, p)
    class(flow_projection), intent(inout) :: this
    type(parameter_list), intent(inout) :: p

    this%p = p
  end subroutine read_params

  subroutine init(this, m)
    class(flow_projection), intent(inout) :: this
    type(flow_mesh), pointer, intent(in) :: m
    !-
    integer :: j, i
    type(pcsr_graph), pointer :: g
    type(pcsr_matrix), pointer :: A
    type(ip_desc), pointer :: row_ip

    this%mesh => m

    allocate(this%rhs(m%ncell))

    !! Create a CSR matrix graph for the pressure poisson system.
    allocate(g)
    row_ip => m%cell_ip
    call g%init (row_ip)
    do j = 1, m%ncell_onP
      call g%add_edge(j,j)
      associate (cn => m%cnhbr(m%xcnhbr(j):m%xcnhbr(j+1)-1))
        do i = 1, size(cn)
          if (cn(i) > 0) call g%add_edge(j, i)
        end do
      end associate
    end do
    call g%add_complete

    allocate(A)
    call A%init(g, take_graph=.true.)
    call this%solver%init(A, this%p)

  end subroutine init


  subroutine setup(this, props)
    class(flow_projection), intent(inout) :: this
    type(flow_props), intent(in) :: props
    !-
    type(pcsr_matrix), pointer :: A
    type(unstr_mesh), pointer :: m
    integer :: i, j, fi, ni
    real(r8) :: coeff, dx(3), length2, n

    A = this%solver%matrix()
    call A%set_all(0.0_r8)

    m => this%mesh%mesh

    do j = 1, m%ncell_onP
      associate (cn => m%cnhbr(m%xcnhbr(j):m%xcnhbr(j+1)-1), &
          fn => m%xcface(j):m%xcface(j+1)-1)

        ASSERT(size(cn) == size(fn))

        do i = 1, size(cn)
          fi = fn(i) ! face index
          ni = cn(i) ! neighbor index

          !! FIXME: COME BACK AND HANDLE SOLID FACES, AND SOLID/VOID CELLS
          if (ni > 0) then
            dx = this%mesh%cell_centroid(:,ni) - this%mesh%cell_centroid(:,j)
            length2 = sum(dx**2)
            ! note that normal is already weighted by face area
            if (btest(m%cfpar(j), pos=i)) then
              n = dot_product(dx, m%normal(:,fi))
            else
              n = -dot_product(dx, m%normal(:,fi))
            end if
            coeff = n/(props%rho_fc(fi)*length2)
            ! DOUBLE CHECK THE SIGN OF THESE COEFFICIENTS
            call A%add_to(j, j, coeff)
            call A%add_to(j, ni, -coeff)
          else
            dx = this%mesh%face_centroid(:,fi) - this%mesh%cell_centroid(:,j)
            length2 = sum(dx**2)
            ! note that normal is already weighted by face area
            if (btest(m%cfpar(j), pos=i)) then
              n = dot_product(dx, m%normal(:,fi))
            else
              n = -dot_product(dx, m%normal(:,fi))
            end if
            coeff = n/(props%rho_fc(fi)*length2)

            call A%add_to(j, j, coeff)
          end if
        end do
      end associate
    end do

    call this%solver%setup()

    this%rhs = 0.0_r8


  end subroutine setup


  subroutine solve(this, solution)
    class(flow_projection), intent(inout) :: this
    real(r8), intent(inout) :: solution(:)
    !-
    integer :: ierr

    call this%solver%solve(this%rhs, solution, ierr)
    if (ierr /= 0) call tls_error("projection solve unsuccessful")

  end subroutine solve

end module flow_type

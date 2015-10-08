MODULE PRECONDITIONERS
  !=============================================================================
  ! Purpose:
  !
  !    Perform preconditioning of linear systems of equations by
  !    solving a related system of linear equations.  In the overall
  !    algorithmic notation, this related system of equations is given
  !    the notation Mw = r.  Here M is a matrix resulting from
  !    approximating the discrete operators embodied in the real
  !    matrix A with simple, orthogonal, regular mesh-connectivity
  !    assumptions.  The routines in this module are general linear
  !    algebra routines, and the notation for the equations is of the
  !    form Ax = b.
  !
  !   Public Interface:
  !
  !     * call PRECONDITION (B, X, status)
  !
  !     This is a general driver and switch, that selects which specific
  !     preconditioner to apply to the system of equations.  Here b is the
  !     right hand side, x is the solution vector (which is overwritten),
  !     and status returns a success/failure code.
  !
  ! Contains: PRECONDITION
  !
  !           BLI
  !           GET_A_COEFFICIENTS_FULL
  !           GET_A_COEFFICIENTS_SPARSE
  !           GET_B_COEFFICIENTS
  !           ILU0_FACTOR_A
  !           ILU0_FACTOR_B
  !           BUBBLE_PERMUTE
  !           ILU0_PRECONDITION
  !           ILU0_SOLVE
  !           INVERSE_DISTANCE_PROLONGATION
  !           JACOBI
  !           LU_FACTOR_A_RF
  !           LU_FACTOR_A_CF
  !           LU_FACTOR_B
  !           LU_PRECONDITION
  !           LU_SOLVE
  !           PIECEWISE_LINEAR_PROLONGATION
  !           REPORT_NORMS
  !           RESIDUAL
  !           SSOR
  !           TM_DIAG
  !           TM_SSOR
  !
  ! The matrix A (or M, if you prefer) is stored as the transpose of
  ! what you would normally expect in ITPACK-ELLPACK form.  In
  ! addition, the diagonal elements are stored in column 0, and there
  ! is no map pointer (in Mesh%Ngbr_Cell) pointing at them.  The column
  ! data is stored in Mesh%Ngbr_Cell, but you generally will collect
  ! them with an EE_Gather call, rather than doing the offset
  ! arithmetic directly using Mesh%Ngbr_Cell.
  !
  ! Mesh%Ngbr_Cell (MNC) contains the map that translates the matrix
  ! A.  if an entry in MNC is 0, then there is no valid coefficient in
  ! that entry of A.  if an MNC entry is > 0, then the entry
  ! represents the index into the X%values array where the
  ! on-processor data that the A coefficient multiplies is.  if an MNC
  ! entry is < 0, then the entry represents the index into the X%aux1
  ! array where the off-processor data that the A coefficient
  ! multiplies is.
  !
  ! The solid mechanics package uses separate routines and preconditioning matrices.
  ! TM_SSOR and TM_DIAG are for the solid mechanics(thermo-mechanics) displacement
  ! solutions.  TM_P and TM_P_Map are the preconditioning matix and connectivity map,
  ! and are ragged arrays (var_vector).
  !
  ! Authors: Douglas B. Kothe, LANL (dbk@lanl.gov)
  !          John Turner, LANL (turner@lanl.gov)
  !          Bryan Lally, LANL (lally@lanl.gov)
  !          Dave Korzekwa, LANL (dak@lanl.gov) (solid mechanics only)
  !
  !=============================================================================
  use kinds, only: r8
  use var_vector_module
  use UbikSolve_module
  use timing_tree
  use truchas_logging_services
  implicit none
  private

  public :: PRECONDITION

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Preconditioning matrix P. The calling physics routine MUST point this to
  ! a matrix of dimension (0:nfc,nunk) which approximates the actual physics,
  ! where nunk = number of unknowns (an integer multiple of nnodes or ncells)

  real(r8), pointer, public, save, dimension(:,:) :: P => null()

  ! Simple diagonal scaling using the diagonal of the preconditioner
  real(r8), pointer, public, save, dimension(:) :: DIAG_P => null()

  ! Preconditioning matrix and connectivity map for solid mechanics (thermo-mechanics)
  ! preconditioners.  Since the number of nodes adjacent to a given node can be highly
  ! variable, we use ragged arrays (var_vector) for these matrices.
  type(real_var_vector), pointer, public, save, dimension(:) :: TM_P => null()
  type(int_var_vector), pointer, public, save, dimension(:)  :: TM_P_Map => null()

CONTAINS

  ! <><><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><><><>

  SUBROUTINE PRECONDITION (B, X_vec, status)
    !===========================================================================
    ! Purpose:
    !
    !   Preconditioner driver
    !===========================================================================
    use linear_solution, only: Ubik_solver, &
                               PRECOND_DIAGONAL, PRECOND_JACOBI, &
                               PRECOND_SSOR, PRECOND_ILU0, PRECOND_LU, &
                               PRECOND_NONE, &
                               PRECOND_TM_SSOR, PRECOND_TM_DIAG
 
    ! arguments
    real(r8), dimension(:), target, intent(IN) :: B
    type(Ubik_vector_type), intent(INOUT) :: X_vec
    integer, intent(OUT) :: status

    ! local variables
    integer :: n1, n2, n3
    real(r8), dimension(:), pointer :: X =>NULL()
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the timer.
    call start_timer ("Solver TMP1")

    ! Alias for solution vector values.
    X => Ubik_values_ptr(X_vec)

    ! Check to make sure P has been associated, but only if we are preconditioning
    ! This is not currently smart enough to tell which preconditioner size and format 
    ! we are supposed to be using (temperature for ncells vs. displacements for (ndim*nnodes),
    ! and may fail to detect a problem.
    if (Ubik_solver%precond /= PRECOND_NONE) then
       if ( (.not. ASSOCIATED(P))      .and. &
            (.not. ASSOCIATED(DIAG_P)) .and. &
            (.not. ASSOCIATED(TM_P))   ) then
          call TLS_panic ('PRECONDITION: Preconditioning matrix not initialized')
       else 
          n1 = 0; n2 = 0; n3 = 0
          if (ASSOCIATED(P)) then
             n1 = SIZE(P,1)
             n2 = SIZE(P,2)
          end if
          if (ASSOCIATED(DIAG_P)) then
             n3 = SIZE(DIAG_P)
          end if
          if (ASSOCIATED(TM_P)) then
             n3 = SIZE(TM_P,1)
          end if
          if ((n2 /= SIZE(B)).and.(n3 /= SIZE(B))) call TLS_panic ('PRECONDITION: Preconditioning matrix size not correct')
       end if
    end if

    select case (Ubik_solver%precond) 
    case (PRECOND_NONE)                 ! no preconditioning
       X = B
    case (PRECOND_DIAGONAL)              ! Simple diagonal scaling
       call DIAGONAL(DIAG_P, X_vec, B, Ubik_steps(Ubik_solver%control), &
            Ubik_omega(Ubik_solver%control), Ubik_solver%precond_scope)
       Ubik_solver%precond_iter = &
            Ubik_solver%precond_iter + Ubik_steps(Ubik_solver%control)
    case (PRECOND_JACOBI)               ! multistep weighted-Jacobi
       call JACOBI (P, X_vec, B, Ubik_steps(Ubik_solver%control), &
            Ubik_omega(Ubik_solver%control), Ubik_solver%precond_scope)
       Ubik_solver%precond_iter = &
            Ubik_solver%precond_iter + Ubik_steps(Ubik_solver%control)
    case (PRECOND_SSOR)                 ! multistep weighted Gauss-Seidel
       call SSOR (P, X_vec, B, Ubik_steps(Ubik_solver%control), &
            Ubik_omega(Ubik_solver%control), Ubik_solver%precond_scope)
       Ubik_solver%precond_iter = &
            Ubik_solver%precond_iter + Ubik_steps(Ubik_solver%control)
    case (PRECOND_TM_SSOR)              ! SSOR for thermo-mechanics only
       call TM_SSOR (TM_P, TM_P_Map, X_vec, B, Ubik_steps(Ubik_solver%control), &
            Ubik_omega(Ubik_solver%control), Ubik_solver%precond_scope)
       Ubik_solver%precond_iter = &
            Ubik_solver%precond_iter + Ubik_steps(Ubik_solver%control)
    case (PRECOND_TM_DIAG)              ! Diagonal scaling for thermo-mechanics only
       call TM_DIAG (TM_P, X_vec, B, Ubik_steps(Ubik_solver%control), &
            Ubik_omega(Ubik_solver%control), Ubik_solver%precond_scope)
       Ubik_solver%precond_iter = &
            Ubik_solver%precond_iter + Ubik_steps(Ubik_solver%control)
    case (PRECOND_ILU0)                 ! ILU(0)
       call ILU0_PRECONDITION (P, X, B, Ubik_solver)
    case (PRECOND_LU)                   ! full LU
       call LU_PRECONDITION (P, X, B, Ubik_solver)
    case default                        ! error
       call TLS_panic ('PRECONDITION: internal error - invalid preconditioner')
    end select

    ! Stop the timer.
    call stop_timer ("Solver TMP1")
    
    status = 0

  END SUBROUTINE PRECONDITION

  ! <><><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><><><>

  SUBROUTINE Y_EQ_AX_ORTHOGONAL (X_vec, Y, status)
    !=======================================================================
    ! Purpose:
    !
    !   Compute y = Ax where A is a matrix stored as array P(nfc,:)
    !   in ELL format.
    !=======================================================================
    use gs_module,        only: EE_GATHER
    use parameter_module, only: nfc, ncells

    ! Arguments
    type(Ubik_vector_type), target, intent(INOUT) :: X_vec
    real(r8), dimension(:), target, intent(INOUT) :: Y
    integer, intent(OUT)   :: status

    ! Local Variables
    real(r8), dimension(nfc,ncells) :: X_Neighbors
    integer :: f
    real(Ubik_real_type), dimension(:), pointer :: X

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the timer.
    call start_timer ("Solver TMP2")

    X => Ubik_values_ptr(X_vec)

    ! Gather element face neighbors
    call EE_GATHER (X_Neighbors, X, BOUNDARY=X_vec%aux1)
!    aux1_ptr => Ubik_aux1_ptr(X_vec)
!    call EE_GATHER (X_Neighbors, X, BOUNDARY=aux1_ptr)
!    call Ubik_set_aux1_ptr(X_vec, aux1_ptr)

    ! Perform the matrix-vector multiply
    Y = P(0,:) * X
    do f = 1,nfc
       Y = Y + P(f,:) * X_Neighbors(f,:)
    end do

    ! Stop the timer.
    call stop_timer ("Solver TMP2")

    status = 0

  END SUBROUTINE Y_EQ_AX_ORTHOGONAL

  SUBROUTINE DIAGONAL(A, X_vec, B, steps, omega, scope)
    !===========================================================================
    ! Purpose:
    !
    !   Diagonal scaling using diagonal elements of the preconditioning matrix
    !===========================================================================
    use linear_solution,  only: PRECOND_SCOPE_LOCAL, PRECOND_SCOPE_GLOBAL

    ! arguments
    real(r8), dimension(:), intent(IN) :: A
    type(Ubik_vector_type), target, intent(INOUT) :: X_vec
    real(r8), dimension(:), intent(IN) :: B
    integer,  intent(IN) :: steps
    real(r8), intent(IN) :: omega
    integer,  intent(IN) :: scope

    ! local variables
    integer :: i,j
    real(r8), dimension(UBOUND(A,1)) :: Dot
    real(r8), dimension(:), pointer  :: X

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    X => Ubik_values_ptr(X_vec)

    do i = 1, steps

       select case (scope)
       case (PRECOND_SCOPE_GLOBAL)
          call Ubik_destroy (X_vec, overlap_only = .true.)
       case (PRECOND_SCOPE_LOCAL)
       case Default
          call TLS_panic ('DIAGONAL: unknown preconditioner scope')
       end select

       ! Diagonal scaling
       Dot = 0.0_r8
       !
       ! update iterate
       do j = 1,UBOUND(A,1)
          Dot(j) = (b(j) - Dot(j))/A(j)
       end do
       X = X + omega*(Dot - X)

    end do

  END SUBROUTINE DIAGONAL

  SUBROUTINE JACOBI (A, X_vec, B, steps, omega, scope)
    !===========================================================================
    ! Purpose:
    !
    !   perform Jacobi iterations
    !===========================================================================
    use linear_solution,  only: PRECOND_SCOPE_LOCAL, PRECOND_SCOPE_GLOBAL
    use mesh_module,      only: DEGENERATE_FACE, Mesh
    use parameter_module, only: nfc, ncells
    use gs_module,        only: EE_GATHER

    ! arguments
    real(r8), dimension(0:,:), intent(IN) :: A
    type(Ubik_vector_type), target, intent(INOUT) :: X_vec
    real(r8), dimension(:), intent(IN) :: B
    integer,  intent(IN) :: steps
    real(r8), intent(IN) :: omega
    integer,  intent(IN) :: scope

    ! local variables
    integer :: f
    integer :: i
    integer :: c
    real(r8), dimension(nfc,ncells) :: X_Neighbors
    integer :: neq
    real(r8), dimension(:), pointer :: Dot
    real(r8), dimension(:), pointer :: X

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    X => Ubik_values_ptr(X_vec)

    neq = size(X)
    allocate(Dot(neq))

    do i = 1, steps

       select case (scope)
       case (PRECOND_SCOPE_GLOBAL)
          call Ubik_destroy (X_vec, overlap_only = .true.)
       case (PRECOND_SCOPE_LOCAL)
       case Default
          call TLS_panic ('JACOBI: unknown preconditioner scope')
       end select

       ! gather element face neighbors, communicating if necessary
       call EE_GATHER (X_Neighbors, X, BOUNDARY=X_vec%aux1)
!       aux1_ptr => Ubik_aux1_ptr(X_vec)
!       call EE_GATHER (X_Neighbors, X, BOUNDARY=aux1_ptr)
!       call Ubik_set_aux1_ptr(X_vec, aux1_ptr)

       ! perform a matrix-vector multiply using coefficients from the orthogonal
       ! operator, except for the diagonal
       Dot = 0.0_r8
       do f = 1, nfc
          do c = 1,ncells
             if (Mesh(c)%Ngbr_Cell(f) == DEGENERATE_FACE) then
                cycle
             end if
             Dot(c) = Dot(c) + A(f,c) * X_Neighbors(f,c)
          end do
       end do

       ! update iterate
       Dot = (b - Dot) / A(0,:)
       X = X + omega*(Dot - X)

    end do

    ! Deallocate -- nicer if this were allocated once
    deallocate(Dot)

  END SUBROUTINE JACOBI

  SUBROUTINE xSSORx (A, X_vec, B, steps, omega, scope)
    !===========================================================================
    ! Purpose:
    !
    !   perform symmetric successive-over-relaxation (SSOR) iterations
    !===========================================================================
    use gs_module,        only: Gather_BoundaryData
    use linear_solution,  only: PRECOND_SCOPE_LOCAL, PRECOND_SCOPE_GLOBAL
    use mesh_module,      only: Mesh, DEGENERATE_FACE
    use parameter_module, only: nfc, ncells

    ! arguments
    real(r8), dimension(0:,:), intent(IN) :: A
    type(Ubik_vector_type), target, intent(INOUT) :: X_vec
    real(r8), dimension(:), intent(IN) :: B
    integer,  intent(IN) :: steps
    real(r8), intent(IN) :: omega
    integer,  intent(IN) :: scope

    ! local variables
    integer :: i
    integer :: c
    integer :: f
    real(r8) :: dot
    real(r8), dimension(nfc) :: LoopTmp
    integer, parameter :: BOUNDARY = 0
    real(r8), dimension(:), pointer :: X

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    X => Ubik_values_ptr(X_vec)

    do i = 1, steps

       select case (scope)
       case (PRECOND_SCOPE_GLOBAL)      ! communicate boundary data every iteration
          call Ubik_destroy (X_vec, overlap_only = .true.)
       case (PRECOND_SCOPE_LOCAL)       ! use what we've got for boundary data
       case DEFAULT
          call TLS_panic ('SSOR: unknown preconditioner scope')
       end select

       ! gather boundary data, communicating if necessary, as determined by the status of X%Aux1
       call GATHER_BOUNDARYDATA (SOURCE=X, BOUNDARY=X_vec%aux1)
!       aux1_ptr => Ubik_aux1_ptr(X_vec)
!       call GATHER_BOUNDARYDATA (SOURCE=X, BOUNDARY=aux1_ptr)
!       call Ubik_set_aux1_ptr(X_vec, aux1_ptr)

       ! forward loop
       do c = 1, ncells
          do f = 1, nfc
             select case(Mesh(c)%Ngbr_cell(f))
             case (DEGENERATE_FACE+1:-1)
                LoopTmp(f) = X_vec%aux1(-Mesh(c)%Ngbr_cell(f))
!                LoopTmp(f) = Ubik_aux1_elem(X_vec, -Mesh(c)%Ngbr_cell(f))
             case (BOUNDARY)
                LoopTmp(f) = 0.0d0
             case (DEGENERATE_FACE)
                LoopTmp(f) = 0.0d0
             case (1:)
                LoopTmp(f) = X(Mesh(c)%Ngbr_cell(f))
             end select
          end do

          dot = DOT_PRODUCT (A(1:,c), LoopTmp)
          dot = (B(c) - dot) / A(0,c)
          X(c) = X(c) + omega*(dot - X(c))
       end do

       ! reverse loop
       do c = ncells, 1, -1
          do f = 1, nfc
             select case(Mesh(c)%Ngbr_cell(f))
             case (DEGENERATE_FACE+1:-1)
                LoopTmp(f) = X_vec%aux1(-Mesh(c)%Ngbr_cell(f))
!                LoopTmp(f) = Ubik_aux1_elem(X_vec, -Mesh(c)%Ngbr_cell(f))
             case (BOUNDARY)
                LoopTmp(f) = 0.0d0
             case (DEGENERATE_FACE)
                LoopTmp(f) = 0.0d0
             case (1:)
                LoopTmp(f) = X(Mesh(c)%Ngbr_cell(f))
             end select
          end do

          dot = DOT_PRODUCT (A(1:,c), LoopTmp)
          dot = (B(c) - dot) / A(0,c)
          X(c) = X(c) + omega*(dot - X(c))
       end do

    end do

  END SUBROUTINE xSSORx

  SUBROUTINE SSOR (A, X_vec, B, steps, omega, scope)
    !===========================================================================
    ! Purpose:
    !
    !   perform symmetric successive-over-relaxation (SSOR) iterations
    !===========================================================================
    use gs_module,            only: Gather_BoundaryData
    use linear_solution,      only: PRECOND_SCOPE_LOCAL, PRECOND_SCOPE_GLOBAL
    use mesh_module,          only: Mesh
    use parameter_module,     only: nfc
    use two_level_partition,  only: Cell_Two_Level_Partitioning, &
                                    Cell_Cell_Two_Level_Partitioned, &
                                    Cell_Cell_Two_Level_Edges_Part
    use pgslib_module,        only: Get_Num_Partitions_Available, &
                                    Get_Start_Available,          &
                                    Get_End_Available, &
                                    PARTITION_INVALID_PARTITIONS, &
                                    Get_Head_Partitions, &
                                    Get_Head_Local, &
                                    Get_Head_Available

    ! arguments
    real(r8), dimension(0:,:), intent(IN) :: A
    type(Ubik_vector_type), target, intent(INOUT) :: X_vec
    real(r8), dimension(:), intent(IN) :: B
    integer,  intent(IN) :: steps
    real(r8), intent(IN) :: omega
    integer,  intent(IN) :: scope

    ! local variables
    integer  :: i
    integer  :: c
    integer  :: f
    integer  :: p
    integer  :: nPartitions_Avail
    integer  :: Edge
    integer  :: NgbrCell
    real(r8) :: dot
    real(r8), dimension(nfc) :: LoopTmp
    integer, parameter :: BOUNDARY = 0
    integer :: neq

    real(r8), dimension(:), pointer :: x_old
    real(r8), dimension(:), pointer :: X

    ! This stuff gets loaded by accessing the partitioning
    
    integer, pointer, dimension(:) :: partition_start
    integer, pointer, dimension(:) :: partition_end
    integer, pointer, dimension(:) :: NgbrPartitions
    logical, pointer, dimension(:) :: NgbrLocal
    logical, pointer, dimension(:) :: NgbrAvailable
    integer, pointer, dimension(:) :: edge_partition_start
    integer, pointer, dimension(:) :: edge_partition_end

    !---------------------------------------------------------------------------

    X => Ubik_values_ptr(X_vec)

    neq = size(X)
    allocate(x_old(neq))

    ! This gets access to the two-level partitiong
    nPartitions_Avail = Get_Num_Partitions_Available(Cell_Two_Level_Partitioning)
    partition_start => get_Start_Available(Cell_Two_Level_Partitioning)
    partition_end   => Get_End_Available  (Cell_Two_Level_Partitioning)
    NgbrPartitions  => Get_Head_Partitions(Cell_Cell_Two_Level_Partitioned)
    NgbrLocal       => Get_Head_Local     (Cell_Cell_Two_Level_Partitioned)
    NgbrAvailable   => Get_Head_Available (Cell_Cell_Two_Level_Partitioned)
    edge_partition_start => get_Start_Available(Cell_Cell_Two_Level_Edges_Part)
    edge_partition_end   => Get_End_Available  (Cell_Cell_Two_Level_Edges_Part)

    ! We need to keep a copy of the input data around since
    ! we might access data in a partition which has already been updated.
    x_old = X

    do i = 1, steps

       select case (scope)
       case (PRECOND_SCOPE_GLOBAL)      ! communicate boundary data every iteration
          call Ubik_destroy (X_vec, overlap_only=.true.)
          x_old = X
       case (PRECOND_SCOPE_LOCAL)       ! use what we've got for boundary data
       end select

       ! gather boundary data, communicating if necessary, as determined by the status of X%Aux1
       call GATHER_BOUNDARYDATA (SOURCE=X, BOUNDARY=X_vec%aux1)
!       aux1_ptr => Ubik_aux1_ptr(X_vec)
!       call GATHER_BOUNDARYDATA (SOURCE=X, BOUNDARY=aux1_ptr)
!       call Ubik_set_aux1_ptr(X_vec, aux1_ptr)

       do p = 1, nPartitions_Avail
          
          ! We also need to count edges, since that is how we query edge partitions.
          ! The way we step through a partitioned graph is "edge" at a time.
          ! We put the edges in in the "right" order, so we just step through
          ! edges to access all the matrix elements (cell neighbors).
          edge = edge_partition_start(p)
          ! forward loop
          do c = partition_start(p), partition_end(p)
             do f = 1, nfc

                ! We also need to know the neighbor cell number.
                ! Should be able to access this from the partitioned graph, but can't yet.
                NgbrCell = Mesh(c)%Ngbr_cell(f)

                ! Now figure out whether the matrix element is between
                ! local partitions, etc...
                if (NgbrPartitions(Edge) == PARTITION_INVALID_PARTITIONS) then
                   ! If this edge isn't, then contributes nothing
                   LoopTmp(f) = 0.0d0

                else if ( NgbrLocal(Edge)) then
                   ! If the matrix element is between cells in the same partition, do
                   ! the standard SSOR thing
                   LoopTmp(f) = X(NgbrCell)

                else if ( NgbrAvailable(Edge) ) then
                   ! If the matrix element is between cells on the same processor
                   ! then we need to use an old value of the data
                   LoopTmp(f) = X_Old(NgbrCell)

                else
                   ! If we got here, then we have off-processor data
                   LoopTmp(f) = X_vec%aux1(-NgbrCell)
!                   LoopTmp(f) = Ubik_aux1_elem(X_vec, -NgbrCell)

                end if

                edge = edge + 1
             end do

             dot = DOT_PRODUCT (A(1:,c), LoopTmp)
             if(ABS(A(0,c)) > 0.0_r8) then
                 dot = (B(c) - dot) / A(0,c)
             endif
             X(c) = X(c) + omega*(dot - X(c))
          end do

          ! reverse loop
          ! Start with the last edge, work backwards.
          edge = edge_partition_end(p)
          do c = partition_end(p), partition_start(p), -1
             do f = nfc, 1, -1

                ! We also need to know the neighbor cell number.
                ! Should be able to access this from the partitioned graph, but can't yet.
                NgbrCell = Mesh(c)%Ngbr_cell(f)

                ! Now figure out whether the matrix element is between
                ! local partitions, etc...
                if (NgbrPartitions(Edge) == PARTITION_INVALID_PARTITIONS) then
                   ! If this edge isn't, then contributes nothing
                   LoopTmp(f) = 0.0d0

                else if ( NgbrLocal(Edge)) then
                   ! If the matrix element is between cells in the same partition, do
                   ! the standard SSOR thing
                   LoopTmp(f) = X(NgbrCell)

                else if ( NgbrAvailable(Edge) ) then
                   ! If the matrix element is between cells on the same processor
                   ! then we need to use an old value of the data
                   LoopTmp(f) = X_Old(NgbrCell)

                else
                   ! If we got here, then we have off-processor data
                   LoopTmp(f) = X_vec%aux1(-NgbrCell)
!                   LoopTmp(f) = Ubik_aux1_elem(X_vec, -NgbrCell)

                end if

                ! The way we step through a partitioned graph is "edge" at a time.
                ! We put the edges in in the "right" order, so we just step through
                ! edges to access all the matrix elements (cell neighbors).
                edge = edge - 1
             end do

             dot = DOT_PRODUCT (A(1:,c), LoopTmp)
             if(ABS(A(0,c)) > 0.0_r8) then
                 dot = (B(c) - dot) / A(0,c)
             endif
             X(c) = X(c) + omega*(dot - X(c))
          end do
       end do
    end do

    deallocate(x_old)

  END SUBROUTINE SSOR

  !-----------------------------------------------------------------------------

  SUBROUTINE LU_PRECONDITION (A_arg, X, B_arg, Ubik)
    !===========================================================================
    ! Purpose:
    !
    !   do a full LU factorization - FOR TESTING ONLY!!
    !   this could be really expensive
    !===========================================================================
    use linear_solution,  only: Ubik_type
    use parameter_module, only: nfc, ncells

    ! Arguments
    real(r8), dimension(0:,:), intent(IN)    :: A_arg
    real(r8), dimension(:),    intent(INOUT) :: X
    real(r8), dimension(:),    intent(IN)    :: B_arg
    type(Ubik_type), target, intent(INOUT)   :: Ubik

    ! Local Variables
    logical, save :: first_time = .true.
    real(r8), dimension(:,:), allocatable, save :: A
    real(r8), dimension(:),   allocatable, save :: B
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! if the factor flag is true, perform the LU factorization,
    ! and set the flag to false (for repeated substitutions with the
    ! same LU factorization), apply the multipliers to b, and back
    ! substitute.  if false, assume the factorization is valid, and
    ! just apply the multipliers to b and back substitute.
    !
    ! we make copies of A and b, in case the caller is doing something
    ! with them and doesn't expect them to be changed

    ! make sure we have room to work
    if (first_time) then
       allocate(A(ncells,ncells), stat=status)
       if (status /= 0) call TLS_panic ('LU_PreCondition: allocation error: A')
       allocate(B(ncells), stat=status)
       if (status /= 0) call TLS_panic ('LU_PreCondition: allocation error: B')
       first_time = .false.
    end if

    ! if FACTOR is true, this is a new matrix problem
    if (Ubik%Factor) then
       ! fill full matrix with local values from sparse matrix
       call GET_A_COEFFICIENTS_FULL (nfc, ncells, A, A_arg)
       ! factor A
       call LU_FACTOR_A_CF (A)
       ! mark it as having been factored
       Ubik%Factor = .false.
    end if

    ! fill b with rhs
    call GET_B_COEFFICIENTS (B, B_arg)

    ! factor b (apply elimination steps stored in A)
    call LU_FACTOR_B (A, B)

    ! back substitute (find new x)
    call LU_SOLVE (A, X, B)

  END SUBROUTINE LU_PRECONDITION

  SUBROUTINE ILU0_PRECONDITION (A_arg, X, B_arg, Ubik)
    !===========================================================================
    ! Purpose:
    !
    !   do ILU0 preconditioning
    !===========================================================================
    use linear_solution,  only: Ubik_type
    use parameter_module, only: nfc, ncells

    ! arguments
    real(r8), dimension(0:,:), intent(IN)    :: A_arg
    real(r8), dimension(:),    intent(INOUT) :: X
    real(r8), dimension(:),    intent(IN)    :: B_arg
    type(Ubik_type), target,   intent(INOUT) :: Ubik

    ! local variables
    logical, save :: first_time = .true.
    real(r8), dimension(:,:), allocatable, save :: A
    real(r8), dimension(:),   allocatable, save :: B
    integer :: status

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! if the factor flag is true, perform the LU factorization,
    ! and set the flag to false (for repeated substitutions with the
    ! same LU factorization), apply the multipliers to b, and back
    ! substitute.  if false, assume the factorization is valid, and
    ! just apply the multipliers to b and back substitute.
    !
    ! we make copies of A and b, in case the caller is doing something
    ! with them and doesn't expect them to be changed

    ! make sure we have room to work
    if (first_time) then
       allocate(A(0:nfc,ncells), stat=status)
       if (status /= 0) call TLS_panic ('ILU0_PreCondition: allocation error: A')
       allocate(B(ncells), stat=status)
       if (status /= 0) call TLS_panic ('ILU0_PreCondition: allocation error: B')
       first_time = .false.
    end if

    ! if FACTOR is true, this is a new matrix problem
    if (Ubik%Factor) then
       ! fill sparse matrix copy with local values from sparse matrix
       call GET_A_COEFFICIENTS_SPARSE (nfc, ncells, A, A_arg)
       ! factor A
       call ILU0_FACTOR_A (nfc, ncells, A)
       ! mark it as having been factored
       Ubik%Factor = .false.
    end if

    ! fill b with rhs
    call GET_B_COEFFICIENTS (B, B_arg)

    ! factor b (apply elimination steps stored in A)
    call ILU0_FACTOR_B (nfc, ncells, A, B)

    ! back substitute (find new x)
    call ILU0_SOLVE (nfc, ncells, A, X, B)

  END SUBROUTINE ILU0_PRECONDITION

  SUBROUTINE TM_DIAG (A, X_vec, B, steps, omega, scope)
    !===========================================================================
    ! Purpose:
    !
    !   Diagonal scaling using diagonal elements of the preconditioning matrix
    !===========================================================================
    use linear_solution,  only: PRECOND_SCOPE_LOCAL, PRECOND_SCOPE_GLOBAL

    ! arguments
    type(real_var_vector), dimension(:), intent(IN) :: A
    type(Ubik_vector_type), target, intent(INOUT) :: X_vec
    real(r8), dimension(:), intent(IN) :: B
    integer,  intent(IN) :: steps
    real(r8), intent(IN) :: omega
    integer,  intent(IN) :: scope

    ! local variables
    integer :: i,j
    real(r8), dimension(UBOUND(A,1)) :: Dot
    real(r8), dimension(:), pointer :: X

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    call start_timer("Precondition")
    X => Ubik_values_ptr(X_vec)

    do i = 1, steps

       select case (scope)
       case (PRECOND_SCOPE_GLOBAL)
          call Ubik_destroy (X_vec, overlap_only = .true.)
       case (PRECOND_SCOPE_LOCAL)
       case Default
          call TLS_panic ('TM_DIAG: unknown preconditioner scope')
       end select

       ! Diagonal scaling
       Dot = 0.0_r8
       !
       ! update iterate
       do j = 1,UBOUND(A,1)
          Dot(j) = (b(j) - Dot(j)) / A(j)%v(1)
       end do
       X = X + omega*(Dot - X)

    end do

    call stop_timer("Precondition")
    
  END SUBROUTINE TM_DIAG

  SUBROUTINE TM_SSOR (A, A_Map, X_vec, B, steps, omega, scope)
    !===========================================================================
    ! Purpose:
    !
    !  For thermo-mechanics only, perform symmetric successive-over-relaxation 
    ! (SSOR) iterations.
    !===========================================================================
    use gs_module,        only: NN_Gather_BoundaryData
    use linear_solution,  only: PRECOND_SCOPE_LOCAL, PRECOND_SCOPE_GLOBAL
    Use parameter_module, Only: ndim, nnodes

    ! arguments
    type(real_var_vector), dimension(:), intent(IN) :: A
    type(int_var_vector), dimension(:), intent(IN) :: A_Map
    type(Ubik_vector_type), target, intent(INOUT) :: X_vec
    real(r8), dimension(:), intent(IN) :: B
    integer,  intent(IN) :: steps
    real(r8), intent(IN) :: omega
    integer,  intent(IN) :: scope

    ! local variables
    integer :: i, j
    integer :: c
    integer :: f
    integer :: status
    real(r8) :: dot
    real(r8), dimension(:), pointer   :: LoopTmp
    integer, parameter                :: BOUNDARY = 0
    real(r8), dimension(:), pointer   :: X
    real(r8), dimension(:), pointer   :: A_Vec
    integer, dimension(:), pointer    :: A_Map_Vec

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    call start_timer("Precondition")

    X => Ubik_values_ptr(X_vec)

    allocate(LoopTmp(MAXVAL(SIZES(A_Map))), stat = status)
    if (status /= 0) call TLS_panic ('TM_SSOR: allocation error: LoopTmp')

    do i = 1, steps
       ! Initially X is all zeros, so we could skip the first nn_gather step and 
       ! set dot products to zero when i = 1.  Maybe later...
       select case (scope)
       case (PRECOND_SCOPE_GLOBAL)      ! communicate boundary data every iteration
          call Ubik_destroy (X_vec, overlap_only = .true.)
       case (PRECOND_SCOPE_LOCAL)       ! use what we've got for boundary data
       case DEFAULT
          call TLS_panic ('TM_SSOR: unknown preconditioner scope')
       end select

       ! gather boundary data, communicating if necessary, as determined by the status of X%Aux2
       ! Since we have ndim * nnodes degrees of freedom, we call NN_GATHER_BOUNDARYDATA ndim
       ! times and store the off processor data in X_vec%aux2
       if (.not. ASSOCIATED (X_vec%aux2)) then
          do j = 1,ndim
             call NN_GATHER_BOUNDARYDATA (SOURCE=X(j:ndim*(nnodes-1)+j:ndim), BOUNDARY=X_vec%aux1)
             if (j == 1) then
                allocate(X_vec%aux2(ndim * SIZE(X_vec%aux1)), stat = status)
                if (status /= 0) call TLS_panic ('TM_SSOR: allocation error: X_vec%aux2')
             end if
             X_vec%aux2(j:ndim*(SIZE(X_vec%aux1)-1)+j:ndim) = X_vec%aux1
             DEALLOCATE(X_vec%aux1)          
          end do
       end if

       ! forward loop
       do c = 1, UBOUND(A,1)
          LoopTmp(:) = 0.0
          A_Map_Vec => FLATTEN(A_Map(c))
          A_Vec => FLATTEN(A(c))
          do f = 2, SIZES(A_Map(c))
             if (A_Map_Vec(f) > 0) then
                LoopTmp(f) = X(A_Map_Vec(f))
             else
                LoopTmp(f) = X_vec%aux2(-A_Map_Vec(f))
             end if
          end do
          dot = DOT_PRODUCT (A_Vec(2:), LoopTmp(2:SIZES(A_Map(c))))
          dot = (B(c) - dot) / A_Vec(1)
          X(c) = X(c) + omega*(dot - X(c))
       end do

       ! reverse loop
       do c = UBOUND(A,1),1, -1
          LoopTmp(:) = 0.0
          A_Map_Vec => FLATTEN(A_Map(c))
          A_Vec => FLATTEN(A(c))
          do f = 2, SIZES(A_Map(c))
             if (A_Map_Vec(f) > 0) then
                LoopTmp(f) = X(A_Map_Vec(f))
             else
                LoopTmp(f) = X_vec%aux2(-A_Map_Vec(f))
             end if
          end do
          dot = DOT_PRODUCT (A_Vec(2:), LoopTmp(2:SIZES(A_Map(c))))
          dot = (B(c) - dot) / A_Vec(1)
          X(c) = X(c) + omega*(dot - X(c))
       end do

    end do

    deallocate(LoopTmp)
    deallocate(X_vec%aux2)
    
    call stop_timer("Precondition")

  END SUBROUTINE TM_SSOR

  SUBROUTINE PIECEWISE_LINEAR_PROLONGATION (X_f, X_c, Delta, ndomains, thisdomain)
    !===========================================================================
    ! Purpose:
    !
    !   perform piecewise linear prolongation, in a crufty manner.  This code
    !   is ugly.  It is of the use-once variety.  It works.
    !
    !   *** ONLY FOR THE 2D POISSON PROBLEM ***
    !===========================================================================
    use mesh_module,      only: Cell
    use parameter_module, only: ncells
    use pgslib_module,    only: PGSLIB_BCAST

    ! Arguments.
    real(r8), dimension(:), intent(INOUT) :: X_f
    real(r8), dimension(:), intent(INOUT) :: X_c
    real(r8), dimension(:), intent(INOUT) :: Delta
    integer, intent(IN) :: ndomains
    integer, intent(IN) :: thisdomain

    ! Local Variables.
    integer :: n, i
    real(r8)  :: h, x, y

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! set up constants
    n = IFIX(SQRT(FLOAT(nDomains)))
    h = 1.0_r8 / float(n)
    ! distribute X_c to all processors
    call PGSLib_BCAST (X_c)

    ! loop over all the cells and calculate X_f
    do i = 1, ncells
       X = (cell(i)%centroid(1) - MOD(thisDomain-1,n) * h) / h
       Y = (cell(i)%centroid(2) - (thisDomain-1)/n * h) / h

       if (thisDomain == 1) then
          ! corner 1
          if (X < 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(0.0_r8, 0.0_r8, 0.0_r8, X_c(1), X*2.0, Y*2.0)
          else if (X >= 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(0.0_r8, 0.0_r8, X_c(1), X_c(2), X-0.5, Y*2.0)
          else if (X < 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(0.0_r8, X_c(1), 0.0_r8, X_c(n+1), X*2.0, Y-0.5)
          else if (X >= 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(1), X_c(2), X_c(n+1), X_c(n+2), X-0.5, Y-0.5)
          end if
       else if (thisDomain == n) then
          ! corner 2
          if (X < 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(0.0_r8, 0.0_r8, X_c(n-1), X_c(n), X+0.5, Y*2.0)
          else if (X >= 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(0.0_r8, 0.0_r8, X_c(n), 0.0_r8, (X-0.5)*2.0, Y*2.0)
          else if (X < 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(n-1), X_c(n), X_c(2*n-1), X_c(2*n), X+0.5, Y-0.5)
          else if (X >= 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(n), 0.0_r8, X_c(2*n), 0.0_r8, (X-0.5)*2.0, Y-0.5)
          end if
       else if (thisDomain == nDomains - n + 1) then
          ! corner 3
          if (X < 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(0.0_r8, X_c(nDomains-2*n+1), 0.0_r8, X_c(nDomains-n+1), X*2.0, Y+0.5)
          else if (X >= 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(X_c(nDomains-2*n+1), X_c(nDomains-2*n+2), X_c(nDomains-n+1), X_c(nDomains-n+2), X-0.5, Y+0.5)
          else if (X < 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(0.0_r8, X_c(nDomains-n+1), 0.0_r8, 0.0_r8, X*2.0, (Y-0.5)*2.0)
          else if (X >= 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(nDomains-n+1), X_c(nDomains-n+2), 0.0_r8, 0.0_r8, X-0.5, (Y-0.5)*2.0)
          end if
       else if (thisDomain == nDomains) then
          ! corner 4
          if (X < 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(X_c(nDomains-n-1), X_c(nDomains-n), X_c(nDomains-1), X_c(nDomains), X+0.5, Y+0.5)
          else if (X >= 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(X_c(nDomains-n), 0.0_r8, X_c(nDomains), 0.0_r8, (X-0.5)*2.0, Y+0.5)
          else if (X < 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(nDomains-1), X_c(nDomains), 0.0_r8, 0.0_r8, X+0.5, (Y-0.5)*2.0)
          else if (X >= 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(nDomains), 0.0_r8, 0.0_r8, 0.0_r8, (X-0.5)*2.0, (Y-0.5)*2.0)
          end if
       else if ((thisDomain-1)/n == 0) then
          ! side 1
          if (X < 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(0.0_r8, 0.0_r8, X_c(thisDomain-1), X_c(thisDomain), X+0.5, Y*2.0)
          else if (X >= 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(0.0_r8, 0.0_r8, X_c(thisDomain), X_c(thisDomain+1), X-0.5, Y*2.0)
          else if (X < 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(thisDomain-1), X_c(thisDomain), X_c(thisDomain+n-1), X_c(thisDomain+n), X+0.5, Y-0.5)
          else if (X >= 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(thisDomain), X_c(thisDomain+1), X_c(thisDomain+n), X_c(thisDomain+n+1), X-0.5, Y-0.5)
          end if
       else if ((thisDomain-1)/n == n-1) then
          ! side 2
          if (X < 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(X_c(thisDomain-n-1), X_c(thisDomain-n), X_c(thisDomain-1), X_c(thisDomain), X+0.5, Y+0.5)
          else if (X >= 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(X_c(thisDomain-n), X_c(thisDomain-n+1), X_c(thisDomain), X_c(thisDomain+1), X-0.5, Y+0.5)
          else if (X < 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(thisDomain-1), X_c(thisDomain), 0.0_r8, 0.0_r8, X+0.5, (Y-0.5)*2.0)
          else if (X >= 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(thisDomain), X_c(thisDomain+1), 0.0_r8, 0.0_r8, X-0.5, (Y-0.5)*2.0)
          end if
       else if (MOD(thisDomain-1,n) == 0) then
          ! side 3
          if (X < 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(0.0_r8, X_c(thisDomain-n), 0.0_r8, X_c(thisDomain), X*2.0, Y+0.5)
          else if (X >= 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(X_c(thisDomain-n), X_c(thisDomain-n+1), X_c(thisDomain), X_c(thisDomain+1), X-0.5, Y+0.5)
          else if (X < 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(0.0_r8, X_c(thisDomain), 0.0_r8, X_c(thisDomain+n), X*2.0, Y-0.5)
          else if (X >= 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(thisDomain), X_c(thisDomain+1), X_c(thisDomain+n), X_c(thisDomain+n+1), X-0.5, Y-0.5)
          end if
       else if (MOD(thisDomain-1,n) == n-1) then
          ! side 4
          if (X < 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(X_c(thisDomain-n-1), X_c(thisDomain-n), X_c(thisDomain-1), X_c(thisDomain), X+0.5, Y+0.5)
          else if (X >= 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(X_c(thisDomain-n), 0.0_r8, X_c(thisDomain), 0.0_r8, (X-0.5)*2.0, Y+0.5)
          else if (X < 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(thisDomain-1), X_c(thisDomain), X_c(thisDomain+n-1), X_c(thisDomain+n), X+0.5, Y-0.5)
          else if (X >= 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(thisDomain), 0.0_r8, X_c(thisDomain+n), 0.0_r8, (X-0.5)*2.0, Y-0.5)
          end if
       else
          ! internal
          if (X < 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(X_c(thisDomain-n-1), X_c(thisDomain-n), X_c(thisDomain-1), X_c(thisDomain), X+0.5, Y+0.5)
          else if (X >= 0.5 .And. Y < 0.5) then
             Delta(i) = BLI(X_c(thisDomain-n), X_c(thisDomain-n+1), X_c(thisDomain), X_c(thisDomain+1), X-0.5, Y+0.5)
          else if (X < 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(thisDomain-1), X_c(thisDomain), X_c(thisDomain+n-1), X_c(thisDomain+n), X+0.5, Y-0.5)
          else if (X >= 0.5 .And. Y >= 0.5) then
             Delta(i) = BLI(X_c(thisDomain), X_c(thisDomain+1), X_c(thisDomain+n), X_c(thisDomain+n+1), X-0.5, Y-0.5)
          end if
       end if
    end do

    X_f = X_f + Delta

  END SUBROUTINE PIECEWISE_LINEAR_PROLONGATION

  SUBROUTINE INVERSE_DISTANCE_PROLONGATION (X_f, X_c, ndomains, thisdomain)
    !===========================================================================
    ! Purpose:
    !
    !   Perform inverse-distance prolongation.
    !
    !===========================================================================
    use cutoffs_module,       only: alittle
    use mesh_module,          only: Cell
    use parameter_module,     only: ndim, ncells
    use pgslib_module,        only: PGSLib_BCAST, PGSLib_COLLATE

    ! Arguments
    integer, intent(IN) :: ndomains
    integer, intent(IN) :: thisdomain
    real(r8), dimension(:), intent(INOUT) :: X_f
    real(r8), dimension(:), intent(INOUT) :: X_c

    ! Local Variables
    integer :: n, p
    real(r8) :: volume
    real(r8), dimension(ndim)          :: Centroid
    real(r8), dimension(ndim,ndomains) :: All_Centroids
    real(r8), dimension(ncells)        :: Distance, Weight, X_cell

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Distribute X_c to all processors.
    call PGSLib_BCAST (X_c)

    ! Compute the geometric centroid on this domain.
    volume   = SUM(Cell%Volume)
    do n = 1,ndim
       Centroid(n) = SUM(Cell%Centroid(n)*Cell%Volume) / Volume
    end do

    ! Collate and broadcast centroids to all processors.
    All_Centroids = 0.0_r8
    do n = 1,ndim
       call PGSLib_COLLATE (All_Centroids(n,:), Centroid(n))
       call PGSLib_BCAST (All_Centroids(n,:))
    end do

    ! Compute cell coarse grid corrections in cell based on the
    ! distance from each cell centroid to each processor centroid.
    X_cell = 0.0_r8
    Weight = 0.0_r8
    do p = 1,ndomains
       Distance = 0.0_r8
       do n = 1,ndim
          Distance = Distance + (Cell%Centroid(n) - All_Centroids(n,p))**2
       end do
       Distance = 1.0_r8/(SQRT(Distance) + alittle)
       X_cell = X_cell + Distance*X_c(p)
       Weight  = Weight + Distance
    end do
    X_cell = X_cell/Weight

    ! Add in the coarse grid correction.
    X_f = X_f + X_cell

  END SUBROUTINE INVERSE_DISTANCE_PROLONGATION

  REAL FUNCTION BLI (v1, v2, v3, v4, x, y)
    !===========================================================================
    ! Purpose:
    !
    !    Perform a bilinear interpolation on the unit square
    !===========================================================================
    ! Arguments
    real(r8), intent(IN) :: v1
    real(r8), intent(IN) :: v2
    real(r8), intent(IN) :: v3
    real(r8), intent(IN) :: v4
    real(r8), intent(IN) :: x
    real(r8), intent(IN) :: y

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    BLI = v1 * (1.0_r8 - x) * (1.0_r8 - y) &
                        + v2 * x       * (1.0_r8 - y) &
                        + v3 * (1.0_r8 - x) * y       &
                        + v4 * x       * y

  END FUNCTION BLI

  SUBROUTINE RESIDUAL (Res, A, X_vec, B, scope)
    !===========================================================================
    ! Purpose:
    !
    !   determine the residual from our equation system, b - Ax
    !===========================================================================
    use gs_module,        only: EE_GATHER
    use linear_solution,  only: PRECOND_SCOPE_LOCAL, PRECOND_SCOPE_GLOBAL
    use parameter_module, only: nfc, ncells

    ! Arguments
    real(r8), dimension(ncells), intent(OUT) :: Res
    real(r8), dimension(0:nfc,ncells), intent(IN) :: A
    type(Ubik_vector_type), target, intent(INOUT) :: X_vec
    real(r8), dimension(ncells), intent(IN) :: B
    integer, intent(IN) :: scope

    ! Local Variables
    integer :: i
    real(r8), dimension(nfc, ncells) :: X_Neighbors
    real(r8), dimension(:), pointer  :: X

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    X => Ubik_values_ptr(X_vec)

    ! local/global specific processing
    select case (scope)
    case (PRECOND_SCOPE_GLOBAL)
       call Ubik_destroy (X_vec, overlap_only = .true.)
    case (PRECOND_SCOPE_LOCAL)
       ! noop - use whatever we have
    case Default
       call TLS_panic ('Residual: unknown preconditioner scope')
    end select

    ! gather element face neighbors
    call EE_GATHER (X_Neighbors, X, BOUNDARY=X_vec%aux1)
!    aux1_ptr => Ubik_aux1_ptr(X_vec)
!    call EE_GATHER (X_Neighbors, X, BOUNDARY=aux1_ptr)
!    call Ubik_set_aux1_ptr(X_vec, aux1_ptr)

    Res = B - A(0,:)*X

    do i = 1, nfc
       Res = Res - A(i,:)*X_Neighbors(i,:)
    end do

  END SUBROUTINE RESIDUAL

  SUBROUTINE GET_A_COEFFICIENTS_FULL (ncoef, n, A, A_arg)
    !===========================================================================
    ! Purpose:
    !===========================================================================
    use mesh_module,      only: Mesh

    ! Arguments
    integer, intent(IN) :: ncoef
    integer, intent(IN) :: n
    real(r8), dimension(n,n), intent(INOUT) :: A
    real(r8), dimension(0:ncoef,n), intent(IN) :: A_arg

    ! Local Variables
    integer :: col, i, row

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! make a dense copy of the local parts of the sparse matrix
    A = 0.0_r8
    do row = 1, n
       A(row,row) = A_arg(0,row)
       do i = 1, ncoef
          col = Mesh(row)%Ngbr_Cell(i)
          if (col > 0) then
             A(row,col) = A_arg(i,row)
          end if
       end do
    end do

  END SUBROUTINE GET_A_COEFFICIENTS_FULL

  SUBROUTINE GET_A_COEFFICIENTS_SPARSE (ncoef, n, A, A_arg)
    !===========================================================================
    ! Purpose:
    !===========================================================================
    use mesh_module,      only: Mesh

    ! Arguments
    integer, intent(IN) :: ncoef
    integer, intent(IN) :: n
    real(r8), dimension(0:ncoef,n), intent(INOUT) :: A
    real(r8), dimension(0:ncoef,n), intent(IN)    :: A_arg

    ! Local Variables
    integer :: i, row

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! make a copy of the local parts of the sparse matrix
    A = 0.0_r8
    do row = 1, n
       A(0,row) = A_arg(0,row)
       do i = 1, ncoef
          if (Mesh(row)%Ngbr_Cell(i) > 0) then
             A(i,row) = A_arg(i,row)
          end if
       end do
    end do

  END SUBROUTINE GET_A_COEFFICIENTS_SPARSE

  SUBROUTINE GET_B_COEFFICIENTS (B, B_arg)

    ! Arguments
    real(r8), dimension(:), intent(INOUT) :: B
    real(r8), dimension(:), intent(IN)    :: B_arg

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! make a copy
    B = B_arg

!!$ RCF Note July 2, 1998
!!$ During merge, seems like Bryan removed the code below.  I'm worried
!!$ that I may have screwed something up during the merge, so I've left
!!$ in the code, but commented it away.
!!$  Bryan, if this code is indeed dead, please delete everything
!!$ which is commented.
!!$<<<<<<< preconditioner_module.F90
!!$    ! if global, fold off-processor contributions into the RHS
!!$    if (scope == PRECOND_SCOPE_GLOBAL) then
!!$       ! if SCOPE is global, ditch current boundary data, if any
!!$       if (ASSOCIATED(X%Aux1)) DEALLOCATE(X%Aux1)
!!$       ! Does this belong here?
!!$       call GATHER_BOUNDARYDATA (BOUNDARY=X%Aux1, &
!!$                                 SOURCE = X%Values)
!!$
!!$
!!$       ! fold the off-processor contributions into the RHS
!!$       do r = 1, ncells
!!$          do f = 1, nfc
!!$             i = mesh(r)%ngbr_cell(f)
!!$             if (i < 0) then
!!$                b(r) = b(r) - A_arg(f,r) * X%aux1(-i)
!!$             end if
!!$          end do
!!$       end do
!!$    end if
!!$
!!$=======
!!$>>>>>>> 1.25

  END SUBROUTINE GET_B_COEFFICIENTS

  SUBROUTINE LU_FACTOR_A_RF (A)
    !===========================================================================
    ! Purpose:
    !
    !   do a full LU factorization - go across rows first
    !===========================================================================
    use cutoffs_module,   only: alittle

    ! Arguments
    real(r8), dimension(:,:), intent(INOUT) :: A

    ! Local Variables
    integer :: row
    integer :: col
    integer :: i
    integer :: n
    real(r8)   :: m

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! factor A into L and U

    ! determine size of system
    n = SIZE(A,1)

    ! for all rows except the first
    do row = 2, n

       ! check for trouble on the diagonal of the row above
       if (Abs(A(row-1,row-1)) < alittle) then
!         write (*,*) 'row, value: ',row-1,A(row-1,row-1)
          call TLS_panic ('LU_Factor_A_RF: zero pivot')
       end if

       ! for all columns left of the diagonal
       do col = 1, row-1

          ! if not a computed zero
          if (Abs(A(row,col)) >= alittle) then

             ! calculate multiplier, and store it in "L"
             m = A(row,col) / A(col,col)
             A(row,col) = m

             ! Scoot across row k, subtracting m * upper row from lower row
             do i = col+1, n
                A(row,i) = A(row,i) - m * A(col,i)
             end do ! i

          end if

       end do ! col
    end do ! row

    ! check for trouble on the diagonal of the last row
    row = n
    if (Abs(A(row,row)) < alittle) then
!      write (*,*) 'row, value: ',row,A(row,row)
       call TLS_panic ('LU_Factor_A_RF: zero pivot')
    end if

  END SUBROUTINE LU_FACTOR_A_RF

  SUBROUTINE LU_FACTOR_A_CF (A)
    !===========================================================================
    ! Purpose:
    !
    !   do a full LU factorization - go down columns first
    !===========================================================================
    use cutoffs_module,   only: alittle

    ! Arguments
    real(r8), dimension(:,:), intent(INOUT) :: A

    ! Local Variables
    integer :: row
    integer :: col
    integer :: i
    integer :: n
    real(r8)   :: m
    real(r8)   :: colmax, amax

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! factor A into L and U

    ! determine size of system
    n = SIZE(A,1)

    ! scan for the maximum diagonal entry of A for pivot test
    ! -- assume entire matrix is on processor
    amax = 0.0
    do col = 1, n
       colmax = ABS(A(col,col))
       if (colmax > amax) amax = colmax
    end do

    ! for all columns except the last
    do col = 1, n - 1

       ! check for trouble on the diagonal
       if (Abs(A(col,col)) < alittle*amax) then
          write (*,*) 'row, value: ',col,A(col,col)
          call TLS_panic ('LU_Factor_A_CF: zero pivot')
       end if

       ! for all rows below the diagonal
       do row = col+1, n

          ! if not a computed zero
          if (Abs(A(row,col)) >= alittle*amax) then

             ! calculate multiplier, and store it in "L"
             m = A(row,col) / A(col,col)
             A(row,col) = m

             ! Scoot across row k, subtracting m * upper row from lower row
             do i = col+1, n
                A(row,i) = A(row,i) - m * A(col,i)
             end do ! i

          end if

       end do ! col
    end do ! row

    ! check for trouble on the diagonal of the last row
    row = n
    if (Abs(A(row,row)) < alittle*amax) then
       write (*,*) 'row, value: ',row,A(row,row)
       call TLS_panic ('LU_Factor_A_CF: zero pivot')
    end if

  END SUBROUTINE LU_FACTOR_A_CF

  SUBROUTINE LU_FACTOR_B (A, b)
    !===========================================================================
    ! Purpose:
    !
    !   do a full LU factorization
    !===========================================================================

    ! Arguments
    real(r8), dimension(:,:), intent(IN)    :: A
    real(r8), dimension(:),   intent(INOUT) :: b

    ! Local Variables
    integer :: row
    integer :: col

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do row = 2, SIZE(b,1)
       do col = 1, row-1
          b(row) = b(row) - A(row,col) * b(col)
       end do
    end do

  END SUBROUTINE LU_FACTOR_B

  SUBROUTINE LU_SOLVE (A, x, b)
    !===========================================================================
    ! Purpose:
    !
    !   do a full LU factorization
    !===========================================================================
    ! Arguments
    real(r8), dimension(:,:), intent(IN)    :: A
    real(r8), dimension(:),   intent(INOUT) :: x
    real(r8), dimension(:),   intent(IN)    :: b

    ! Local Variables
    integer :: row
    integer :: col
    real(r8)   :: sum

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! back substitute (find new x)
    do row = SIZE(A,1), 1, -1
       sum = 0.0_r8
       do col = row+1, SIZE(A,1)
          sum = sum + A(row,col) * x(col)
       end do
       x(row) = (b(row) - sum) / A(row,row)
    end do

  END SUBROUTINE LU_SOLVE

  SUBROUTINE ILU0_FACTOR_A (ncoef, n, A)
    !===========================================================================
    ! Purpose:
    !
    !   do an Incomplete LU factorization, fill in 0
    !===========================================================================
    use mesh_module,      only: Mesh
    use cutoffs_module,   only: alittle
    use utilities_module, only: bubble_permute

    ! Arguments
    integer, intent(IN) :: ncoef
    integer, intent(IN) :: n
    real(r8), dimension(0:ncoef,n), intent(INOUT) :: A

    ! Local Variables
    integer,  dimension(ncoef) :: p    ! permutation vector
    integer  :: i    ! index
    integer  :: j    ! index
    integer  :: k    ! index
    integer  :: face ! face number
    integer  :: row  ! row number
    integer  :: col  ! column number
    integer  :: colj ! column number
    real(r8) :: m    ! multiplier

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! factor A into L and U (sort of - in the ILU(0) way)
    do row = 2, n        ! for all rows except the first

       ! check for trouble on the diagonal of the row above
       if (Abs(A(0,row-1)) < alittle) &
          call TLS_panic ('ILU0: zero pivot')

       ! calculate an ordered permutation vector, p
       p = BUBBLE_PERMUTE(Mesh(row)%Ngbr_Cell(:))

       ! for all columns in sparse A
       do face = 1, ncoef

          i = p(face)                           ! index into map and data arrays
          col = Mesh(row)%Ngbr_Cell(i)          ! actual matrix column number

          if (col > 0                        &  ! if not a structural zero
             .And. col < row                 &  ! and left of diagonal
             .And. Abs(A(i,row)) >= alittle) &  ! and not a computed zero
          then

             ! calculate multiplier, and store it in "L"
             m = A(i,row) / A(0,col)
             A(i,row) = m

             ! Scoot across row k, subtracting m * upper row from
             ! lower row.  if there is no matching column entry in
             ! the upper row, this is a noop.  Don't forget the
             ! diagonal element, which must be special cased in our
             ! data structure
             do j = face, ncoef ! for each column left in lower row

                ! lower row column number
                colj = Mesh(row)%Ngbr_Cell(p(j))

                ! search upper row to see if there is a matching column number
                do k = 1, ncoef
                   if (Mesh(col)%Ngbr_Cell(k) == colj) then
                      ! if we find a match, multiply and subtract elements
                      A(p(j),row) = A(p(j),row) - m * A(k,col)
                      Exit           ! there can only be one match
                   end if
                end do

             end do ! j

             ! special case for diagonal element at A(0,row)
             colj = row

             ! search upper row to see if there is a matching column number
             do k = 1, ncoef
                if (Mesh(col)%Ngbr_Cell(k) == colj) then
                   ! if we find a match, multiply and subtract elements
                   A(0,row) = A(0,row) - m * A(k,col)
                   Exit              ! there can only be one match
                end if
             end do
          end if
       end do ! face
    end do ! row

  END SUBROUTINE ILU0_FACTOR_A

  SUBROUTINE ILU0_FACTOR_B (ncoef, n, A, B)
    !===========================================================================
    ! Purpose:
    !
    !   do an Incomplete LU factorization, fill in 0
    !===========================================================================
    use mesh_module,      only: Mesh
    use utilities_module, only: bubble_permute

    ! Arguments
    integer, intent(IN) :: ncoef
    integer, intent(IN) :: n
    real(r8), dimension(0:,:), intent(IN)    :: A
    real(r8), dimension(:),    intent(INOUT) :: B

    ! Local Variables
    integer :: row  ! row number
    integer :: col  ! column number
    integer :: face ! index
    integer :: i    ! index
    integer, dimension(ncoef) :: p    ! permutation vector

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do row = 2, n
       p = BUBBLE_PERMUTE(Mesh(row)%Ngbr_Cell(:))
       do face = 1, ncoef
          i = p(face)
          col = Mesh(row)%Ngbr_Cell(i)

          if (col > 0 .and. col < row) then
             b(row) = b(row) - A(i,row) * b(col)
          end if
       end do
    end do

  END SUBROUTINE ILU0_FACTOR_B

  SUBROUTINE ILU0_SOLVE (ncoef, n, A, X, B)
    !===========================================================================
    ! Purpose:
    !
    !   do an Incomplete LU factorization, fill in 0
    !===========================================================================
    use mesh_module,      only: Mesh

    ! Arguments
    integer, intent(IN) :: ncoef
    integer, intent(IN) :: n
    real(r8), dimension(0:,:), intent(IN)    :: A
    real(r8), dimension(:),    intent(INOUT) :: X
    real(r8), dimension(:),    intent(IN)    :: B

    ! Local Variables
    integer :: row  ! row number
    integer :: col  ! column number
    integer :: i    ! index
    real(r8)   :: sum  ! dot product result

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! back substitute (find new x)
    do row = n, 1, -1
       sum = 0.0_r8
       do i = 1, ncoef
          col = Mesh(row)%Ngbr_Cell(i)
          if (col > row) then
             sum = sum + A(i,row) * X(col)
          end if
       end do
       X(row) = (B(row) - sum) / A(0,row)
    end do

  END SUBROUTINE ILU0_SOLVE

END MODULE PRECONDITIONERS

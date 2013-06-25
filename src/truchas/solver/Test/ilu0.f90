! ilu0.f90

!-------------------------------------------------------------------------------
! test code for the ILU0 procedure in src/solver/preconditioner_module.F90
!
! instructions:
!    make ILU0 public in src/solver/preconditioner_module.F90, and
!    build the entire code (successfully), so that you have valid .mod files
!    and lib*.a.  Then run ./ilu0.build.  This will compile
!    and run the test.
!-------------------------------------------------------------------------------

Program TEST_ILU0

  Use kind_module,           Only: int_kind, real_kind
  Use JTpack_module,         Only: JT_vector_type
  Use parameter_module,      Only: ncells, nfc
  Use solver_data_module,    Only: JT
  Use mesh_module,           Only: Mesh
  Use cutoffs_module,        Only: alittle
  Use preconditioner_module, Only: ILU0

  Implicit None  

  Integer (KIND=int_kind),  Parameter          :: N = 9
  Integer (KIND=int_kind)                      :: i
  Real    (KIND=real_kind), Dimension(0:nfc,N) :: A
  Real    (KIND=real_kind), Dimension(N)       :: b
  Real    (KIND=real_kind), Dimension(N)       :: answer
  Type    (JT_vector_type)                     :: x

  ! dimension the test problem
  ncells = N

  ! allocate some space for x, in the proper type
  Allocate(x%Values(1:N))
  Allocate(Mesh(1:N))

  ! set the factor flag to TRUE, so the matrix gets factored
  JT%ILU_Factor = .true.

  ! build the sparse matrices, in the secret almost-ELL
  ! format - do it by hand
  Mesh(1)%Ngbr_Cell(1) = 2
  Mesh(1)%Ngbr_Cell(2) = 4
  Mesh(1)%Ngbr_Cell(3) = 0
  Mesh(1)%Ngbr_Cell(4) = 0
  Mesh(1)%Ngbr_Cell(5) = 0
  Mesh(1)%Ngbr_Cell(6) = 0

  A(0,1) = 4
  A(1,1) = -1
  A(2,1) = -1
  A(3,1) = 0
  A(4,1) = 0
  A(5,1) = 0
  A(6,1) = 0

  Mesh(2)%Ngbr_Cell(1) = 1
  Mesh(2)%Ngbr_Cell(2) = 3
  Mesh(2)%Ngbr_Cell(3) = 5
  Mesh(2)%Ngbr_Cell(4) = 0
  Mesh(2)%Ngbr_Cell(5) = 0
  Mesh(2)%Ngbr_Cell(6) = 0

  A(0,2) = 4
  A(1,2) = -1
  A(2,2) = -1
  A(3,2) = -1
  A(4,2) = 0
  A(5,2) = 0
  A(6,2) = 0

  Mesh(3)%Ngbr_Cell(1) = 2
  Mesh(3)%Ngbr_Cell(2) = 4
  Mesh(3)%Ngbr_Cell(3) = 6
  Mesh(3)%Ngbr_Cell(4) = 0
  Mesh(3)%Ngbr_Cell(5) = 0
  Mesh(3)%Ngbr_Cell(6) = 0

  A(0,3) = 4
  A(1,3) = -1
  A(2,3) = -1
  A(3,3) = -1
  A(4,3) = 0
  A(5,3) = 0
  A(6,3) = 0

  Mesh(4)%Ngbr_Cell(1) = 3
  Mesh(4)%Ngbr_Cell(2) = 5
  Mesh(4)%Ngbr_Cell(3) = 7
  Mesh(4)%Ngbr_Cell(4) = 1
  Mesh(4)%Ngbr_Cell(5) = 0
  Mesh(4)%Ngbr_Cell(6) = 0

  A(0,4) = 4
  A(1,4) = -1
  A(2,4) = -1
  A(3,4) = -1
  A(4,4) = -2
  A(5,4) = 0
  A(6,4) = 0

  Mesh(5)%Ngbr_Cell(1) = 4
  Mesh(5)%Ngbr_Cell(2) = 6
  Mesh(5)%Ngbr_Cell(3) = 8
  Mesh(5)%Ngbr_Cell(4) = 2
  Mesh(5)%Ngbr_Cell(5) = 0
  Mesh(5)%Ngbr_Cell(6) = 0

  A(0,5) = 4
  A(1,5) = -1
  A(2,5) = -1
  A(3,5) = -1
  A(4,5) = -2
  A(5,5) = 0
  A(6,5) = 0

  Mesh(6)%Ngbr_Cell(1) = 5
  Mesh(6)%Ngbr_Cell(2) = 7
  Mesh(6)%Ngbr_Cell(3) = 9
  Mesh(6)%Ngbr_Cell(4) = 3
  Mesh(6)%Ngbr_Cell(5) = 0
  Mesh(6)%Ngbr_Cell(6) = 0

  A(0,6) = 4
  A(1,6) = -1
  A(2,6) = -1
  A(3,6) = -1
  A(4,6) = -2
  A(5,6) = 0
  A(6,6) = 0

  Mesh(7)%Ngbr_Cell(1) = 6
  Mesh(7)%Ngbr_Cell(2) = 8
  Mesh(7)%Ngbr_Cell(3) = 4
  Mesh(7)%Ngbr_Cell(4) = 2
  Mesh(7)%Ngbr_Cell(5) = 0
  Mesh(7)%Ngbr_Cell(6) = 0

  A(0,7) = 4
  A(1,7) = -1
  A(2,7) = -1
  A(3,7) = -2
  A(4,7) = 3
  A(5,7) = 0
  A(6,7) = 0

  Mesh(8)%Ngbr_Cell(1) = 7
  Mesh(8)%Ngbr_Cell(2) = 9
  Mesh(8)%Ngbr_Cell(3) = 5
  Mesh(8)%Ngbr_Cell(4) = 0
  Mesh(8)%Ngbr_Cell(5) = 0
  Mesh(8)%Ngbr_Cell(6) = 0

  A(0,8) = 4
  A(1,8) = -1
  A(2,8) = -1
  A(3,8) = -2
  A(4,8) = 0
  A(5,8) = 0
  A(6,8) = 0

  Mesh(9)%Ngbr_Cell(1) = 8
  Mesh(9)%Ngbr_Cell(2) = 6
  Mesh(9)%Ngbr_Cell(3) = 0
  Mesh(9)%Ngbr_Cell(4) = 0
  Mesh(9)%Ngbr_Cell(5) = 0
  Mesh(9)%Ngbr_Cell(6) = 0

  A(0,9) = 4
  A(1,9) = -1
  A(2,9) = -2
  A(3,9) = 0
  A(4,9) = 0
  A(5,9) = 0
  A(6,9) = 0

  ! and now a right hand side
  b(1) = 1.0
  b(2) = 2.0
  b(3) = 3.0
  b(4) = 4.0
  b(5) = 5.0
  b(6) = 6.0
  b(7) = 7.0
  b(8) = 8.0
  b(9) = 9.0

  ! here's the right answer
  answer(1) = 2.6907559061276136d0
  answer(2) = 3.7939216219129723d0
  answer(3) = 4.5926649353013644d0
  answer(4) = 5.9691020025974817d0
  answer(5) = 7.3845411468722819d0
  answer(6) = 7.5768470891942759d0
  answer(7) = 6.4441635400945776d0
  answer(8) = 7.8475679056452083d0
  answer(9) = 6.9766384455691473d0

  ! call the procedure
  Call ILU0 (A, x, b)

  ! check the results
  If (Any(Abs(answer-x%Values) > alittle)) Then
     Write (*,*) 'ilu0 test failed'
  Else
     Write (*,*) 'ilu0 test passed'
  End If

End Program TEST_ILU0

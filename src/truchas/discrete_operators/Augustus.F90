!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Module: Augustus
! Author: Michael L. Hall, <Hall@LANL.gov>
!
! This module consists of procedures adapted from the Augustus Diffusion 
! Package (http://www.LANL.gov/Augustus). This particular incarnation was 
! developed specifically for the Truchas code and currently contains 
! only the routines necessary to generate the cell-based Support Operator 
! discretization matrix. See the file Augustus_COPYING for copyright and 
! license restrictions.
! 
! The Augustus Package is based on the methods derived in the following 
! references:
!
! Michael L. Hall and Jim E. Morel, "Diffusion Discretization Schemes in 
!   Augustus: A New Hexahedral Symmetric Support Operator Method", in 
!   Proceedings of the 1998 Nuclear Explosives Code Developers Conference 
!   (NECDC), Las Vegas, NV, October 25-30, 1998. LA-UR-98-3146. Available 
!   online at http://www.LANL.gov/Augustus.
!
! Michael L. Hall, Jim E. Morel, and Mikhail J. Shashkov, "A Local Support 
!   Operator Diffusion Discretization Scheme for Hexahedral Meshes",
!   JOWOG 42 Presentation, October 21, 1999. LA-UR-99-5834. Available online 
!   at http://www.LANL.gov/Augustus.
!
! Jim E. Morel, Michael L. Hall & Mikhail J. Shashkov (2001), "A Local 
!   Support-Operators Diffusion Discretization Scheme for Hexahedral Meshes", 
!   Journal of Computational Physics 170(1):338-372, June 2001. LA-UR-99-4358.
!   Available online at http://www.LANL.gov/Augustus.
!
! This module contains the following public procedure:
!
!   Augustus_Set_Support_Op_Matrix - Sets the cell-based matrix for the 
!                                    Support Operator diffusion discretization.
!
! This module contains the following private procedures:
!
!   Average_Nodes_to_Faces                 - Averages node values to get face 
!                                            values for a cell.
!   Cross_Product (.Cross.)                - Returns the cross-product of two 
!                                            vectors.
!   Determinant_Small_Matrix (Determinant) - Returns the determinant of a 
!                                            small matrix.
!   Hex_Volume_Face_Integral               - Returns an integral used in the 
!                                            calculation of the volume of a 
!                                            hexahedron.
!   Hexahedron_Volume                      - Returns the volume of a 
!                                            hexahedron.
!   Invert_Small_Matrix (Invert)           - Inverts a small matrix.
!   Parallepiped_Volume                    - Returns the volume of a 
!                                            parallepiped.
!   Quadrilateral_Volume                   - Returns the volume associated 
!                                            with a quadrilateral in 2-D.
!   Segment_Volume                         - Returns the volume associated 
!                                            with a line segment in 1-D.
!   SO_Node_Volume_Weights                 - Calculates the node volume 
!                                            weightings for the Support 
!                                            Operator method.
!   Triangle_Area                          - Returns the area of a triangle.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module Augustus

  ! Global use associations.

  use kind_module     ! This is the only dependence on the main Truchas code.
                      ! It defines the intrinsic types (real_kind, int_kind, 
                      ! log_kind) which are used below.

  ! Set counter-clockwise numbering option for 2-D.

# define NUMBER2D_COUNTER_CLOCKWISE

  ! Start up with everything untyped and private.

  implicit none
  private
 
  ! Public procedure.

  public :: Augustus_Set_Support_Op_Matrix

  interface Augustus_Set_Support_Op_Matrix
    module procedure Augustus_Set_Support_Op_Matrix
  end interface

  ! Private procedures.

  interface OPERATOR (.Cross.)
    module procedure Cross_Product
  end interface

  interface Determinant
    module procedure Determinant_Small_Matrix
  end interface

  interface Invert
    module procedure Invert_Small_Matrix
  end interface

  ! Global class variables.

  ! Define number parameters.
  real(real_kind), parameter :: zero=0.d0, one=1.d0,  two=2.d0, three=3.d0, &
                                four=4.d0, twelve=12.d0
  real(real_kind), parameter :: half=one/two, fourth=one/four
  real(real_kind), parameter :: pi=3.141592653589793238462643383279d0, &
                                twothirdspi=two*pi/three, &
                                fourthirdspi=four*pi/three

contains

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! Procedure: Augustus_Set_Support_Op_Matrix
  ! Author:    Michael L. Hall, <Hall@LANL.gov>
  !
  ! This procedure sets the individual cell matrix for the Support Operator
  ! Method, as implemented in the Augustus code package. It does not invert 
  ! this matrix.
  !
  ! The matrix returned by this procedure, S_Cell, is the matrix in the 
  ! following equation,
  !
  !         [      ]   [      ]
  !   S     [ F .A ] = [ dPhi ]
  !    Cell [  f  f]   [     f]
  !
  ! where f signifies a face value, Cell signifies a cell value, and
  !
  !   F_f.A_f  =  a vector of fluxes (F = -D grad Phi) dotted with face 
  !               areas for each face, 
  !   dPhi_f   =  a vector of intensity differences (dPhi_f = Phi_Cell - Phi_f)
  !               for each face,
  !   D        =  a diffusion coefficient defined on each face.
  !
  ! S_Cell is a square matrix with dimensions equal to the number of faces 
  ! per cell (NLocal_Faces).
  !
  ! Several geometric quantities are input rather than calculated here. It 
  ! would be safer and cleaner to calculate them here, but faster to use the 
  ! calculated values from the host code. This may be changed in the future.
  !
  ! Input variables:
  !
  !   NCells_PE                   Number of cells on this PE.
  !   NDimensions                 Number of dimensions (1-3).
  !   Cell_Volume                 Volumes of the cells.
  !   Adjacent_Faces              Local face numbers for the NDimensions faces
  !                               around each node (see Numbering section).
  !   Diff_Coeff_Faces_of_Cells   Diffusion coefficient for each local face
  !                               for all cells (see Numbering section).
  !   Area_Faces_of_Cells         Scalar area of each local face for all cells
  !                               (see Numbering section).
  !   Unit_Normal_Faces_of_Cells  Unit normal vector of each local face
  !                               for all cells (see Numbering section).
  !   Coordinates_Nodes_of_Cells  Coordinates of the local nodes of each cell
  !                               (see Numbering section).
  ! 
  ! Output variable:
  !
  !   S_Cell                      Support Operator matrix in "S F.A = dPhi"
  !                               equation for each cell (see Numbering
  !                               section).
  !
  ! Numbering:
  !
  !   Some of the internal procedures in this module have specific numbering
  !   requirements, which are described in the procedures themselves. However,
  !   the Augustus_Set_Support_Op_Matrix procedure, which is the only public
  !   interface of this module, handles all internal numbering requirements
  !   and has fairly simple numbering requirements itself.
  !
  !   First, the nodes must be numbered in this manner in the
  !   Coordinates_Nodes_of_Cells variable:
  !
  !     1-D geometry.
  !
  !        +-----------+  
  !       1|           |2 
  !        +-----------+
  !
  !     2-D geometry.
  !
  !       4             3
  !        +-----------+       
  !        |           |
  !        |           |
  !        |           |
  !        |           |
  !        |           |
  !        +-----------+
  !       1             2
  !
  !     3-D geometry.
  !
  !           8----------7
  !          /|         /|
  !         / |        / |
  !        /  |       /  |
  !       5----------6   |
  !       |   |      |   |
  !       |   4------|---3
  !       |  /       |  /     
  !       | /        | /
  !       |/         |/
  !       1----------2
  !
  !   Next, the Adjacent_Faces variable can be set using whatever face
  !   numbering the host code uses. It must only specify which face numbers
  !   are adjacent to each node (it is essentially a local connectivity
  !   matrix).
  !
  !   Then, the Diff_Coeff_Faces_of_Cells, Area_Faces_of_Cells, and
  !   Unit_Normal_Faces_of_Cells variables must be specified using a numbering
  !   which is consistent with the host code numbering used by Adjacent_Faces.
  !   The resultant S_Cell variable is also consistent with the specified
  !   host code face numbering.

  subroutine Augustus_Set_Support_Op_Matrix ( &
    S_Cell, Diff_Coeff_Faces_of_Cells, Area_Faces_of_Cells, &
    Unit_Normal_Faces_of_Cells, Adjacent_Faces, Coordinates_Nodes_of_Cells, &
    Cell_Volume, NCells_PE, NDimensions)

    ! Input variables.

    integer(int_kind), intent(in) :: NCells_PE   ! Number of cells on this PE.
    integer(int_kind), intent(in) :: NDimensions ! Number of dimensions.
    ! Volumes of the cells.
    real(real_kind), dimension(NCells_PE), intent(in) :: Cell_Volume
    ! Local face numbers for the NDimensions faces around each node.
    integer(int_kind), dimension(2**NDimensions,NDimensions), intent(in) :: &
      Adjacent_Faces
    ! Diffusion coefficient for each local face for all cells on this PE.
    real(real_kind), dimension(2*NDimensions,NCells_PE), intent(in) :: &
      Diff_Coeff_Faces_of_Cells
    ! Scalar area of each local face for all cells on this PE.
    real(real_kind), dimension(2*NDimensions,NCells_PE), intent(in) :: &
      Area_Faces_of_Cells
    ! Unit normal vector of each local face for all cells on this PE.
    real(real_kind), dimension(NDimensions,2*NDimensions,NCells_PE), &
      intent(in) :: Unit_Normal_Faces_of_Cells
    ! Coordinates of the local nodes of each cell on this PE.
    real(real_kind), dimension(NDimensions,2**NDimensions,NCells_PE), &
      intent(in) :: Coordinates_Nodes_of_Cells

    ! Internal variables.

    integer(int_kind) :: c             ! Cell loop counter.
    integer(int_kind) :: n             ! Local node loop counter.
    integer(int_kind) :: d             ! Dimension loop counter.
    integer(int_kind) :: f             ! Face loop counter.
    integer(int_kind) :: i, j          ! Matrix row and column loop counters.
    integer(int_kind) :: NLocal_Faces  ! Number of local faces in a cell.
    integer(int_kind) :: NLocal_Nodes  ! Number of local nodes in a cell.
    logical(log_kind) :: degenerate    ! Degenerate node flag.
    ! The Node Volume/Weights.
    real(real_kind), dimension(2**NDimensions) :: Node_Volume
    ! The Jacobian, its inverse, & the transpose of the inverse for a cell.
    real(real_kind), dimension(NDimensions,NDimensions) :: Jacobian 
    real(real_kind), dimension(NDimensions,NDimensions) :: J_minus
    real(real_kind), dimension(NDimensions,NDimensions) :: J_minus_Transpose
    ! Unexpanded Support Operator matrix in "S F.A = dPhi" for a single node.
    real(real_kind), dimension(NDimensions,NDimensions) :: S_Node
    ! Diffusion Coefficient evaluated at a node.
    real(real_kind) :: Diff_Coeff_Node

    ! Output variable.

    real(real_kind), dimension(2*NDimensions,2*NDimensions,NCells_PE), &
      intent(out) :: S_Cell  ! Support Operator matrix in "S F.A = dPhi" eqn.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Diagnostic output.
    ! write (6,*) 'Calling Augustus SO'

    ! First, some definitions.

    NLocal_Faces = 2*NDimensions
    NLocal_Nodes = 2**NDimensions

    !--------------------------------------------------
    ! Set the S matrix for each cell (don't invert it).
    !--------------------------------------------------

    S_Cell = zero

    do c = 1, NCells_PE

      ! Calculate node volumes (or weights).

      call SO_Node_Volume_Weights (Node_Volume, Cell_Volume(c), NDimensions, &
                                   Coordinates_Nodes_of_Cells(:,:,c))

      ! Start loop over local nodes in the cell.

      do n = 1, NLocal_Nodes

        ! Check for degeneracy.

        degenerate = Node_Volume(n).eq.zero
        do d = 1, NDimensions
          degenerate = degenerate .or. &
            Area_Faces_of_Cells(Adjacent_Faces(n,d),c) == zero
        end do

        ! If this node is degenerate, skip it.

        if (.not. degenerate) then

          ! Calculate Jacobian matrix for this node. Columns of
          ! the Jacobian are equal to the area vectors for the
          ! adjacent faces.

          do d = 1, NDimensions
            Jacobian(:,d) = Area_Faces_of_Cells(  Adjacent_Faces(n,d),c) * &
                     Unit_Normal_Faces_of_Cells(:,Adjacent_Faces(n,d),c)
          end do

          ! Calculate S_Node = J^-1 J^-T.

          J_minus = Jacobian
          call Invert (J_minus)
          J_minus_Transpose = TRANSPOSE(J_minus)
          S_Node = MATMUL(J_minus, J_minus_Transpose)

          ! Scale by Node_Volume/D_Node. D_Node is defined to be the 
          ! average of the adjacent face values.
  
          Diff_Coeff_Node = zero
          do d = 1, NDimensions
            Diff_Coeff_Node = Diff_Coeff_Node + &
              Diff_Coeff_Faces_of_Cells(Adjacent_Faces(n,d),c)
          end do
          Diff_Coeff_Node = Diff_Coeff_Node / REAL(NDimensions, real_kind)
          ! Correct division by zero problems. 
          Diff_Coeff_Node = MAX(Diff_Coeff_Node, 1.d-100)
          S_Node = S_Node * Node_Volume(n) / Diff_Coeff_Node

          ! Scatter each node matrix, S_Node, to the overall cell 
          ! matrix, S_Cell. This is equivalent to premultiplying
          ! by P^T and postmultiplying by P.
  
          do i = 1, NDimensions
            do j = 1, NDimensions
              S_Cell(Adjacent_Faces(n,i),Adjacent_Faces(n,j),c) = &
              S_Cell(Adjacent_Faces(n,i),Adjacent_Faces(n,j),c) + S_Node(i,j)
            end do
          end do
        end if
      end do

      ! Set diagonals for degenerate faces.

      do f = 1, NLocal_Faces
        if (Area_Faces_of_Cells(f,c) == zero .or. &
            S_Cell(f,f,c) == zero) then
          S_Cell(f,f,c) = one
        end if
      end do
    end do
    return
  end subroutine Augustus_Set_Support_Op_Matrix
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! Procedure: Average_Nodes_to_Faces
  ! Author:    Michael L. Hall, <Hall@LANL.gov>
  !
  ! The Average_Nodes_to_Faces procedure averages a variable that is defined 
  ! on the cell nodes to derive cell face values. This is one means of 
  ! interpolation, but there are others (conservation of quantity*volume, for 
  ! instance).
  !
  ! This procedure depends upon a particular node and cell numbering, which 
  ! are explained in the comments below.
  
  subroutine Average_Nodes_to_Faces (Value_Faces, Value_Nodes, NDimensions)

    ! Input variables.

    integer(int_kind), intent(in) :: NDimensions  ! Number of dimensions.
    ! Values on the nodes, numbered as shown in the comments below.
    real(real_kind), dimension(:,:), intent(in) :: Value_Nodes

    ! Output variable.

    real(real_kind), dimension(:,:), intent(out) :: Value_Faces

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Verify requirements.

    if (SIZE(Value_Nodes,2) /= 2**NDimensions) then
      write (6,*) 'Error (Average_Nodes_to_Faces): Values_Nodes ', &
                  'is incorrectly dimensioned. '
    end if
    if (SIZE(Value_Faces,2) /= 2*NDimensions) then
      write (6,*) 'Error (Average_Nodes_to_Faces): Values_Faces ', &
                  'is incorrectly dimensioned. '
    end if

    ! Average from Nodes to Faces. 

    select case (NDimensions)     ! Toggle on dimensions.
    case (1)

      ! 1-D geometry - Faces and Nodes refer to the same location, 
      !                and NLocal_Faces = NLocal_Nodes.
      !
      !    +-----------+       1 - left face (-k)
      !   1|           |2      2 - right face (+k)
      !    +-----------+       

      Value_Faces = Value_Nodes

    case (2)

#     ifdef NUMBER2D_COUNTER_CLOCKWISE

        ! 2-D geometry (Counter-Clockwise Numbering).
        !
        !   4      3      3
        !    +-----------+       
        !    |           |       4 - left face (-k)
        !    |           |       2 - right face (+k)
        !   4|           |2      1 - bottom face (-l)
        !    |           |       3 - top face (+l)
        !    |           |
        !    +-----------+
        !   1      1      2

        Value_Faces(:,1) = half*(Value_Nodes(:,1) + Value_Nodes(:,2))
        Value_Faces(:,2) = half*(Value_Nodes(:,2) + Value_Nodes(:,3))
        Value_Faces(:,3) = half*(Value_Nodes(:,3) + Value_Nodes(:,4))
        Value_Faces(:,4) = half*(Value_Nodes(:,4) + Value_Nodes(:,1))

#     else

        ! 2-D geometry (Generalized Multi-Dimensional Numbering).
        !
        ! Note that the GMD numbering is used for the faces here but 
        ! not for the nodes, which are still counter-clockwise.
        !
        !   4      4      3
        !    +-----------+       
        !    |           |       1 - left face (-k)
        !    |           |       2 - right face (+k)
        !   1|           |2      3 - bottom face (-l)
        !    |           |       4 - top face (+l)
        !    |           |
        !    +-----------+
        !   1      3      2

        Value_Faces(:,1) = half*(Value_Nodes(:,4) + Value_Nodes(:,1))
        Value_Faces(:,2) = half*(Value_Nodes(:,2) + Value_Nodes(:,3))
        Value_Faces(:,3) = half*(Value_Nodes(:,1) + Value_Nodes(:,2))
        Value_Faces(:,4) = half*(Value_Nodes(:,3) + Value_Nodes(:,4))

#     endif

    case (3)

      ! 3-D geometry.
      !
      !       8----------7
      !      /|         /|    1 - left face (-k)
      !     / |  6     / |    2 - right face (+k)
      !    /  |    4  /  |    3 - front face (-l)
      !   5----------6   |    4 - back face (+l)
      !   | 1 |      | 2 |    5 - bottom face (-m)
      !   |   4------|---3    6 - top face (+m)
      !   |  / 3     |  /     
      !   | /    5   | /
      !   |/         |/
      !   1----------2

      Value_Faces(:,1) = fourth*(Value_Nodes(:,1) + Value_Nodes(:,4) + &
                                 Value_Nodes(:,5) + Value_Nodes(:,8))
      Value_Faces(:,2) = fourth*(Value_Nodes(:,2) + Value_Nodes(:,3) + &
                                 Value_Nodes(:,6) + Value_Nodes(:,7))
      Value_Faces(:,3) = fourth*(Value_Nodes(:,1) + Value_Nodes(:,2) + &
                                 Value_Nodes(:,5) + Value_Nodes(:,6))
      Value_Faces(:,4) = fourth*(Value_Nodes(:,3) + Value_Nodes(:,4) + &
                                 Value_Nodes(:,7) + Value_Nodes(:,8))
      Value_Faces(:,5) = fourth*(Value_Nodes(:,1) + Value_Nodes(:,2) + &
                                 Value_Nodes(:,3) + Value_Nodes(:,4))
      Value_Faces(:,6) = fourth*(Value_Nodes(:,5) + Value_Nodes(:,6) + &
                                 Value_Nodes(:,7) + Value_Nodes(:,8))

    end select

    return
  end subroutine Average_Nodes_to_Faces
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! Procedure: Cross_Product
  ! Author:    Michael L. Hall, <Hall@LANL.gov>
  !
  ! The Cross_Product procedure returns the cross product of two vectors
  ! (c = a X b). Either 2-D or 3-D vectors may be used. For the 2-D
  ! case, the scalar length of the vector in the orthogonal direction,
  ! also known as the "right-hand normal", is returned in c(1).
  !
  ! The procedure has an operator interface, so it may be called in the 
  ! following ways:
  !
  !   c = a .Cross. b
  !   c = Cross_Product(a, b)

  function Cross_Product (a, b) result(c)

    ! Input variables.

    real(real_kind), dimension(:), intent(in) :: a, b  ! Vectors to be crossed.

    ! Output variable.

    real(real_kind), dimension(SIZE(a)) :: c           ! Cross product result.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Verify requirements.

    if (SIZE(a) /= SIZE(b)) then
      write (6,*) 'Error (Cross_Product): vectors of unequal size.'
    end if

    ! Calculate the cross-product.

    select case (SIZE(a))
    case (1)
      write (6,*) 'Error (Cross_Product): vector size is unity.'
    case (2)
      c(1) = a(1)*b(2) - b(1)*a(2)
      c(2) = zero  ! Not used, so set it to zero.
    case (3)
      c(1) = a(2)*b(3) - b(2)*a(3)
      c(2) = a(3)*b(1) - b(3)*a(1)
      c(3) = a(1)*b(2) - b(1)*a(2)
    end select

    return
  end function Cross_Product
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! Procedure: Determinant_Small_Matrix
  ! Author:    Michael L. Hall, <Hall@LANL.gov>
  !
  ! The Determinant_Small_Matrix procedure returns the determinant of a 
  ! small (size<=3) matrix.

  function Determinant_Small_Matrix (M) result(Determinant)

    ! Input variable.

    real(real_kind), dimension(:,:) :: M

    ! Output variable.

    real(real_kind) :: Determinant

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Verify requirements.

    if (SIZE(M,1) > 3 .or. SIZE(M,2) > 3) then
      write (6,*) 'Error (Determinant_Small_Matrix): matrix too large.'
    end if

    ! Toggle on dimensionality.

    select case (SIZE(M,1))

    case (1)

      Determinant = M(1,1)

    case (2)

      Determinant = M(1,1)*M(2,2) - M(1,2)*M(2,1)

    case (3)

      Determinant =   M(1,1)*(M(2,2)*M(3,3) - M(3,2)*M(2,3)) &
                    - M(1,2)*(M(2,1)*M(3,3) - M(3,1)*M(2,3)) &
                    + M(1,3)*(M(2,1)*M(3,2) - M(3,1)*M(2,2))
    end select

    return
  end function Determinant_Small_Matrix
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! Procedure: Hexahedron_Volume
  ! Author:    Michael L. Hall, <Hall@LANL.gov>
  !
  ! The Hexahedron_Volume procedure returns the volume of the hexahedron 
  ! formed by eight points in 3-D. The points should be given in the order 
  ! in this figure to yield the correct volume:
  !
  !           8            7
  !            +----------+
  !           /|         /|
  !          / |        / |
  !       5 /  |      6/  |
  !        +----------+   |
  !        |   |      |   |
  !        | 4 +------|---+ 3
  !        |  /       |  /
  !        | /        | /
  !        |/         |/
  !        +----------+
  !       1            2
  !
  !  This routine will give the correct volume of a degenerate hexahedron, 
  !  such as a wedge (prism), a pyramid, or a tetrahedron, if the points are 
  !  specified in the hexahedral format.
  !
  !  This routine makes use of the method outlined in "Efficient Volume
  !  Computation for Three-Dimensional Hexahedral Cells", John K. Dukowicz, 
  !  J. Comp. Phys., February 1988, Vol. 74, No. 2, pp. 493-496. Briefly, that
  !  method starts from the identity
  !
  !    div r = 3 
  !
  !  (where r is the coordinate vector in three-space), integrates both sides 
  !  over a volume, and then applies Gauss' Theorem to yield:
  !
  !    3 V = Int  r.n dA
  !             S
  !
  !  This translates the volume integral into a series of surface integrals 
  !  which are easily computed. The individual surface integrals are 
  !  computed (and derived) in the routine Hex_Volume_Face_Integral.

  function Hexahedron_Volume (R1, R2, R3, R4, R5, R6, R7, R8) &
                              result(Volume)
    ! Input variables.

    ! Node points of the hexahedron, given in the order shown in the figure 
    ! above.
    real(real_kind), dimension(3), intent(in) :: R1, R2, R3, R4, R5, R6, R7, R8

    ! Output variable.

    real(real_kind) :: Volume  ! Volume of the hexahedron.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Volume = &
      ( &                    
        Hex_Volume_Face_Integral(R1, R2, R3, R4) + &  ! Bottom (-z) face
        Hex_Volume_Face_Integral(R5, R8, R7, R6) + &  ! Top (+z) face
        Hex_Volume_Face_Integral(R1, R4, R8, R5) + &  ! Left (-x) face
        Hex_Volume_Face_Integral(R2, R6, R7, R3) + &  ! Right (+x) face
        Hex_Volume_Face_Integral(R1, R5, R6, R2) + &  ! Front (-y) face
        Hex_Volume_Face_Integral(R3, R7, R8, R4)   &  ! Back (+y) face
      ) / twelve     ! The factor of 12 comes from a factor of 3 shown 
                     ! in the comments above and a factor of 4 from the 
                     ! surface integral calculation.

    ! Verify guarantees.

    if (Volume < zero) then
      write (6,*) 'Error (Hexahedron_Volume): negative volume -- points ', &
                  'may be in wrong order or cell may be malformed.'
    end if

    return
  end function Hexahedron_Volume
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! Procedure: Hex_Volume_Face_Integral
  ! Author:    Michael L. Hall, <Hall@LANL.gov>
  !
  ! The Hex_Volume_Face_Integral procedure returns a surface integral used in 
  ! the volume calculation (Hexahedron_Volume). For a correct calculation of 
  ! the integral, the four points of a hexahedral face must be input in a 
  ! counter-clockwise order when viewed from the interior of the hexahedron.
  !
  ! The integral that is calculated by this routine is:
  !
  !   S = 4 Int  r.n dA
  !            S
  !
  ! When all the surface integrals are summed, the cell volume is then
  !
  !          ---
  !        1 \
  !   V = -- /    S
  !       12 ---
  !          sides
  !
  ! This routine makes use of the method outlined in "Efficient Volume
  ! Computation for Three-Dimensional Hexahedral Cells", John K. Dukowicz, 
  ! J. Comp. Phys., February 1988, Vol. 74, No. 2, pp. 493-496. 
  !
  ! The following comments derive several equations that are stated in the 
  ! aforementioned paper.
  !
  ! The volume that we wish to calculate is
  !
  !   V = 0.333 Int  r.n dA
  !                S
  !
  ! which can be written as (following equation 5 from the paper)
  !
  !                1    1   
  !   V = 0.333 Int  Int  r.(dr/dk X dr/dl) dk dl
  !                0    0
  !
  ! where 
  !
  !   r = R1 k (1-l) + R2 k l + R3 (1-k) l + R4 (1-k) (1-l)
  !
  ! Integrating and applying some vector algebra (details given in the F77 
  ! version of this Augustus procedure) show that the volume integral (using 
  ! triple product notation) is
  !
  !   V = (1/12) { [R3,R2,R1] + [R2,R1,R4] + [R1,R4,R3] + [R4,R3,R2] }
  !
  ! which is the form given in equation 6 of the paper. Further manipulation 
  ! yields
  !
  !   V = (1/12) [R2+R3, R1+R2, R3+R4]
  !
  ! which is the first part of equation 7 (and the equation used in
  ! this procedure). The factor of 1/12 is added by the calling routine.

  function Hex_Volume_Face_Integral (R1, R2, R3, R4) result(Integral)

    ! Input variables.

    ! Four node points of the 3-D face, in counter-clockwise 
    ! order viewed from the interior of the hexahedron.
    real(real_kind), dimension(3), intent(in) :: R1, R2, R3, R4

    ! Output variable.

    real(real_kind) :: Integral  ! Face integral.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Calculate the face integral.

    Integral = DOT_PRODUCT( (R2 + R3),  (R1 + R2).Cross.(R3 + R4) )

    ! Verify guarantees - none.

    return
  end function Hex_Volume_Face_Integral
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! Procedure: Invert_Small_Matrix
  ! Author:    Michael L. Hall, <Hall@LANL.gov>
  !
  ! This routine returns the inverse of a small (size<=3) matrix. This 
  ! routine overwrites the matrix with the inverse (M = M^-1). 

  subroutine Invert_Small_Matrix (M)

    ! Input/Output variable.

    ! Matrix to be transformed.
    real(real_kind), dimension(:,:), intent(inout) :: M

    ! Internal variable.

    ! Copy of original matrix.
    real(real_kind), dimension(SIZE(M,1),SIZE(M,2)) :: Mo

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Verify requirements.

    if (SIZE(M,1) > 3 .or. SIZE(M,2) > 3) then
      write (6,*) 'Error (Invert_Small_Matrix): matrix too large.'
    end if

    ! Save original matrix.

    Mo = M

    ! Toggle on dimensionality.

    select case (SIZE(M,1))

    case (1)

      M(1,1) = one / M(1,1)

    case (2)

      ! Find transposed cofactors.

      M(1,1) = Mo(2,2)
      M(2,2) = Mo(1,1)
      M(1,2) = -Mo(1,2)
      M(2,1) = -Mo(2,1)

      M = M / Determinant (Mo)

    case (3)

      ! Find cofactors.

      M(1,1) =  (Mo(2,2)*Mo(3,3)-Mo(3,2)*Mo(2,3))
      M(1,2) = -(Mo(2,1)*Mo(3,3)-Mo(3,1)*Mo(2,3))
      M(1,3) =  (Mo(2,1)*Mo(3,2)-Mo(3,1)*Mo(2,2))
      M(2,1) = -(Mo(1,2)*Mo(3,3)-Mo(3,2)*Mo(1,3))
      M(2,2) =  (Mo(1,1)*Mo(3,3)-Mo(3,1)*Mo(1,3))
      M(2,3) = -(Mo(1,1)*Mo(3,2)-Mo(3,1)*Mo(1,2))
      M(3,1) =  (Mo(1,2)*Mo(2,3)-Mo(2,2)*Mo(1,3))
      M(3,2) = -(Mo(1,1)*Mo(2,3)-Mo(2,1)*Mo(1,3))
      M(3,3) =  (Mo(1,1)*Mo(2,2)-Mo(2,1)*Mo(1,2))

      ! Transpose and scale.

      M = TRANSPOSE(M) / Determinant (Mo)

    end select

    return
  end subroutine Invert_Small_Matrix
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! Procedure: Parallepiped_Volume 
  ! Author:    Michael L. Hall, <Hall@LANL.gov>
  ! 
  ! The Parallepiped_Volume procedure returns the volume of a parallelepiped, 
  ! when given the four points (a node and three edge points) that define it. 
  ! The center node should be in node point one, and the other three node 
  ! points should be in counter-clockwise order, viewed from the outside of 
  ! the parallelepiped, in order to assure positive volume.
  !
  ! The volume of the parallelepiped is given by the scalar triple product of 
  ! vectors between the node point coordinates (R):
  !
  !  Volume = (R  - R ) . (R  - R ) X (R  - R )
  !           ( 4    1)   ( 3    1)   ( 2    1)

  function Parallepiped_Volume (R1, R2, R3, R4) result(Volume)

    ! Input variables.

    ! Center point of the parallelepiped corner.
    real(real_kind), dimension(3), intent(in) :: R1
    ! Three edge points of the parallelepiped, in counter-clockwise 
    ! order looking from the outside.
    real(real_kind), dimension(3), intent(in) :: R2, R3, R4

    ! Output variable.

    real(real_kind) :: Volume  ! Volume of the parallelepiped.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Calculate the volume.

    Volume = DOT_PRODUCT( (R4 - R1),  (R3 - R1).Cross.(R2 - R1) )

    ! Verify guarantees.

    if (Volume < zero) then
      write (6,*) 'Error (Parallepiped_Volume): negative volume -- points ', &
                  'may be in wrong order or cell may be malformed.'
    end if

    return
  end function Parallepiped_Volume
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! Procedure: Quadrilateral_Volume 
  ! Author:    Michael L. Hall, <Hall@LANL.gov>
  ! 
  ! The Quadrilateral_Volume procedure returns the volume associated with a 
  ! quadrilateral in 2-D. 
  !
  ! For cartesian geometry, the quadrilateral "volume" is equal to the
  ! area (suppressed z), and the area is the sum of two triangular 
  ! areas. 
  ! 
  ! For cylindrical geometry, the *exact* (not approximate) volume of
  ! revolution around the z-axis for a triangle is:
  ! 
  !        2 pi
  !   V = ------ (R  + R  + R ) (Area of the triangle) 
  !         3      1    2    3
  !
  ! where R_i are the r-coordinates of the triangle points. The cell
  ! volume is then given by the sum of the triangular volumes. 
  !
  ! The points should be given in a counter-clockwise order.

  function Quadrilateral_Volume (Geometry, R1, R2, R3, R4) result(Volume)

    ! Input variables.

    ! Four points of the quadrilateral, in counter-clockwise order.
    real(real_kind), dimension(2), intent(in) :: R1, R2, R3, R4
    integer(int_kind), intent(in) :: Geometry           ! Geometry flag: 
                                                        !   0 - Cartesian, 
                                                        !   1 - Cylindrical.

    ! Internal variables.

    ! Output variable.

    real(real_kind) :: Volume  ! Volume of the parallelepiped.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    select case (Geometry)
       case (0)
          Volume = Triangle_Area (R1, R2, R4) + Triangle_Area (R2, R3, R4)
       case (1)
          Volume = twothirdspi * ( (R1(1) + R2(1) + R4(1)) * Triangle_Area (R1, R2, R4) &
                                 + (R2(1) + R3(1) + R4(1)) * Triangle_Area (R2, R3, R4) )
       case default
          write (6,*) 'Error (Quadrilateral_Volume): Geometry ', geometry, ' is not valid.'
    end select

    ! Verify guarantees.

    if (Volume < zero) then
      write (6,*) 'Error (Quadrilateral_Volume): negative volume -- points ', &
                  'may be in wrong order or cell may be malformed.'
    end if

    return
  end function Quadrilateral_Volume
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! Procedure: Segment_Volume
  ! Author:    Michael L. Hall, <Hall@LANL.gov>
  !
  ! The Segment_Volume procedure returns the volume associated with a line 
  ! segment in 1-D. The segment is denoted by two node points (in any order). 
  ! 
  ! In cartesian geometry, the "volume" (representing a slab) is equal to 
  ! the length of the segment:
  ! 
  !   V = | R  - R  |   (suppressed y, z)
  !          2    1
  ! 
  ! In cylindrical geometry, the "volume" (representing a cylindrical shell) 
  ! is equal to the difference in the area of two circles:
  ! 
  !             2    2
  !   V = pi | R  - R  |   (suppressed z)
  !             2    1
  ! 
  ! In spherical geometry, the volume is equal to the spherical shell
  ! volume:
  ! 
  !        4        3    3
  !   V = --- pi | R  - R  |
  !        3        2    1
  !
  ! where R_i are the r-coordinates of the two points.

  function Segment_Volume (Geometry, R1, R2) result(Volume)

    ! Input variables.
    
    real(real_kind), dimension(1), intent(in) :: R1, R2 ! Segment end points.
    integer(int_kind), intent(in) :: Geometry           ! Geometry flag: 
                                                        !   0 - Cartesian, 
                                                        !   1 - Cylindrical, 
                                                        !   2 - Spherical.

    ! Output variable.

    real(real_kind) :: Volume    ! Volume of the segment.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Calculate the volume.

    select case (Geometry)
    case (0)

      ! Cartesian segment (represents a slab).
      Volume = ABS(R1(1) - R2(1))

    case (1)

      ! Cylindrical (represents a cylindrical shell).
      Volume = pi * ABS(R1(1)**2 - R2(1)**2)

    case (2)

      ! Spherical (represents a spherical shell).
      Volume = fourthirdspi * ABS(R1(1)**3 - R2(1)**3)

    end select

    return
  end function Segment_Volume
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! Procedure: SO_Node_Volume_Weights
  ! Author:    Michael L. Hall, <Hall@LANL.gov>
  !
  ! This procedure calculates the nodal volumes or weights for
  ! the Support Operator Method. While strictly these values should
  ! be the volumes associated with each node (segment, quadrant or
  ! octant of a cell), they may be defined in several ways, and are
  ! sometimes therefore referred to as nodal weights. The definitions
  ! available in Augustus are:
  !
  !   Weighting_Scheme = "True Volume"
  !     - Weights are set to the actual segment/quadrant/octant volume for 
  !       the node.
  !
  !   Weighting_Scheme = "Parallel Unscaled"
  !     - Weights are set to the volume of a parallelepiped defined at 
  !       each node which is roughly the same size as the true volume.
  !
  !   Weighting_Scheme = "Parallel Scaled"
  !     - Weights are set to the volume of a parallelepiped defined at 
  !       each node which is roughly the same size as the true volume.
  !       Then, the weights for a cell are scaled so that the sum is equal 
  !       to the true volume for the cell.

  subroutine SO_Node_Volume_Weights (Node_Volumes, Cell_Volume, NDimensions, &
                                     R_node)

    ! Input variables.

    integer(int_kind), intent(in) :: NDimensions  ! Number of dimensions.
    real(real_kind), intent(in) :: Cell_Volume    ! Volume of the cell.
    ! Coordinates for each node of a cell. Should really be called
    ! Coordinates_Nodes_of_Cell, but is called R_node here due to heavy use 
    ! in this procedure.
    real(real_kind), dimension(NDimensions,2**NDimensions), intent(in) :: &
      R_node

    ! Internal variables.

    integer(int_kind) :: NLocal_Nodes      ! Number of local nodes in a cell.
    character(len=17) :: Weighting_Scheme  ! Chooses one of three possible 
                                           ! weighting schemes.
    integer(int_kind) :: Geometry          ! Geometry flag: 
                                           !   0 - Cartesian, 
                                           !   1 - Cylindrical, 
                                           !   2 - Spherical.
    ! Coordinates of the cell center.
    real(real_kind), dimension(NDimensions) :: R_cell
    ! Coordinates of the face centers.
    real(real_kind), dimension(NDimensions,2*NDimensions) :: R_face
    ! Coordinates of the three edge (yes, edge) centers around a node in 3-D.
    real(real_kind), dimension(NDimensions,NDimensions) :: R_edge

    ! Output variable.

    real(real_kind), dimension(2**NDimensions), intent(out) :: Node_Volumes

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    ! First, some definitions.

    NLocal_Nodes = 2**NDimensions
    Geometry = 0   ! Cartesian geometry.

    ! Choose a weighting scheme (all of these schemes are operational).

      Weighting_Scheme = "Parallel Scaled"
    ! Weighting_Scheme = "Parallel Unscaled"
    ! Weighting_Scheme = "True Volume"

    !-----------------------
    ! True Volume Weighting.
    !-----------------------
  
    if (TRIM(Weighting_Scheme) == "True Volume") then

      ! Determine face midpoints.

      call Average_Nodes_to_Faces (R_face, R_node, NDimensions)

      ! Determine cell center value.

      R_cell(:) = SUM(R_node(:,:), DIM=2) / REAL(NLocal_Nodes, real_kind)

      ! Determine volumes associated with each node.

      select case (NDimensions)

      ! 1-D geometry. 
      !
      !    +-----------+       1 - left face (-k)
      !   1|           |2      2 - right face (+k)
      !    +-----------+       

      case(1)

        Node_Volumes(1) = Segment_Volume (Geometry, R_node(:,1), R_cell)
        Node_Volumes(2) = Segment_Volume (Geometry, R_cell, R_node(:,2))

      ! 2-D cases.

      case(2)

        ! 2-D geometry (Counter-Clockwise Numbering).
        !
        !   4      3      3
        !    +-----------+       
        !    |           |       4 - left face (-k)
        !    |           |       2 - right face (+k)
        !   4|           |2      1 - bottom face (-l)
        !    |           |       3 - top face (+l)
        !    |           |
        !    +-----------+
        !   1      1      2

#       ifdef NUMBER2D_COUNTER_CLOCKWISE

          Node_Volumes(1) = Quadrilateral_Volume ( &
            Geometry, R_face(:,4), R_node(:,1), R_face(:,1), R_cell)
          Node_Volumes(2) = Quadrilateral_Volume ( &
            Geometry, R_face(:,1), R_node(:,2), R_face(:,2), R_cell)
          Node_Volumes(3) = Quadrilateral_Volume ( &
            Geometry, R_face(:,2), R_node(:,3), R_face(:,3), R_cell)
          Node_Volumes(4) = Quadrilateral_Volume ( &
            Geometry, R_face(:,3), R_node(:,4), R_face(:,4), R_cell)

        ! 2-D geometry (Generalized Multi-Dimensional Numbering).
        !
        ! Note that the GMD numbering is used for the faces here but 
        ! not for the nodes, which are still counter-clockwise.
        !
        !   4      4      3
        !    +-----------+       
        !    |           |       1 - left face (-k)
        !    |           |       2 - right face (+k)
        !   1|           |2      3 - bottom face (-l)
        !    |           |       4 - top face (+l)
        !    |           |
        !    +-----------+
        !   1      3      2

#       else

          Node_Volumes(1) = Quadrilateral_Volume ( &
            Geometry, R_face(:,1), R_node(:,1), R_face(:,3), R_cell)
          Node_Volumes(2) = Quadrilateral_Volume ( &
            Geometry, R_face(:,3), R_node(:,2), R_face(:,2), R_cell)
          Node_Volumes(3) = Quadrilateral_Volume ( &
            Geometry, R_face(:,2), R_node(:,3), R_face(:,4), R_cell)
          Node_Volumes(4) = Quadrilateral_Volume ( &
            Geometry, R_face(:,4), R_node(:,4), R_face(:,1), R_cell)

#       endif

      ! 3-D geometry.
      !
      !       8----------7
      !      /|         /|    1 - left face (-k)
      !     / |  6     / |    2 - right face (+k)
      !    /  |    4  /  |    3 - front face (-l)
      !   5----------6   |    4 - back face (+l)
      !   | 1 |      | 2 |    5 - bottom face (-m)
      !   |   4------|---3    6 - top face (+m)
      !   |  / 3     |  /     
      !   | /    5   | /
      !   |/         |/
      !   1----------2

      case(3)

        R_edge(:,1) = half * (R_node(:,1) + R_node(:,2))
        R_edge(:,2) = half * (R_node(:,1) + R_node(:,4))
        R_edge(:,3) = half * (R_node(:,1) + R_node(:,5))
        Node_Volumes(1) = Hexahedron_Volume ( &
          R_node(:,1), R_edge(:,1), R_face(:,5), R_edge(:,2), &
          R_edge(:,3), R_face(:,3), R_cell,      R_face(:,1))

        R_edge(:,1) = half * (R_node(:,2) + R_node(:,3))
        R_edge(:,2) = half * (R_node(:,2) + R_node(:,1))
        R_edge(:,3) = half * (R_node(:,2) + R_node(:,6))
        Node_Volumes(2) = Hexahedron_Volume ( &
          R_node(:,2), R_edge(:,1), R_face(:,5), R_edge(:,2), &
          R_edge(:,3), R_face(:,2), R_cell,      R_face(:,3))

        R_edge(:,1) = half * (R_node(:,3) + R_node(:,4))
        R_edge(:,2) = half * (R_node(:,3) + R_node(:,2))
        R_edge(:,3) = half * (R_node(:,3) + R_node(:,7))
        Node_Volumes(3) = Hexahedron_Volume ( &
          R_node(:,3), R_edge(:,1), R_face(:,5), R_edge(:,2), &
          R_edge(:,3), R_face(:,4), R_cell,      R_face(:,2))

        R_edge(:,1) = half * (R_node(:,4) + R_node(:,1))
        R_edge(:,2) = half * (R_node(:,4) + R_node(:,3))
        R_edge(:,3) = half * (R_node(:,4) + R_node(:,8))
        Node_Volumes(4) = Hexahedron_Volume ( &
          R_node(:,4), R_edge(:,1), R_face(:,5), R_edge(:,2), &
          R_edge(:,3), R_face(:,1), R_cell,      R_face(:,4))

        R_edge(:,1) = half * (R_node(:,5) + R_node(:,8))
        R_edge(:,2) = half * (R_node(:,5) + R_node(:,6))
        R_edge(:,3) = half * (R_node(:,5) + R_node(:,1))
        Node_Volumes(5) = Hexahedron_Volume ( &
          R_node(:,5), R_edge(:,1), R_face(:,6), R_edge(:,2), &
          R_edge(:,3), R_face(:,1), R_cell,      R_face(:,3))

        R_edge(:,1) = half * (R_node(:,6) + R_node(:,5))
        R_edge(:,2) = half * (R_node(:,6) + R_node(:,7))
        R_edge(:,3) = half * (R_node(:,6) + R_node(:,2))
        Node_Volumes(6) = Hexahedron_Volume ( &
          R_node(:,6), R_edge(:,1), R_face(:,6), R_edge(:,2), &
          R_edge(:,3), R_face(:,3), R_cell,      R_face(:,2))

        R_edge(:,1) = half * (R_node(:,7) + R_node(:,6))
        R_edge(:,2) = half * (R_node(:,7) + R_node(:,8))
        R_edge(:,3) = half * (R_node(:,7) + R_node(:,3))
        Node_Volumes(7) = Hexahedron_Volume ( &
          R_node(:,7), R_edge(:,1), R_face(:,6), R_edge(:,2), &
          R_edge(:,3), R_face(:,2), R_cell,      R_face(:,4))

        R_edge(:,1) = half * (R_node(:,8) + R_node(:,7))
        R_edge(:,2) = half * (R_node(:,8) + R_node(:,5))
        R_edge(:,3) = half * (R_node(:,8) + R_node(:,4))
        Node_Volumes(8) = Hexahedron_Volume ( &
          R_node(:,8), R_edge(:,1), R_face(:,6), R_edge(:,2), &
          R_edge(:,3), R_face(:,4), R_cell,      R_face(:,1))

      end select

      ! Check to see if sum of node volumes equals cell volume.

      if (ABS((SUM(Node_Volumes) - Cell_Volume) / Cell_Volume) > 1.d-12) then
        write (6,*) 'SUM(Node_Volumes) / Cell_Volume mismatch: '
        write (6,*) '  SUM(Node_Volumes) = ', SUM(Node_Volumes)
        write (6,*) '  Cell_Volume       = ', Cell_Volume
      end if

    end if

    !--------------------------------------
    ! Parallel Scaled and Unscaled Weights.
    !--------------------------------------

    if (TRIM(Weighting_Scheme) == "Parallel Scaled" .or. &
        TRIM(Weighting_Scheme) == "Parallel Unscaled" ) then

      ! Toggle on number of dimensions.

      select case (NDimensions)

      ! Cell Volume. (probably not valid for SO method unless cartesian)

      case (1)

        Node_Volumes(1) = Segment_Volume (Geometry, R_node(:,1), R_node(:,2))
        Node_Volumes(2) = Segment_Volume (Geometry, R_node(:,1), R_node(:,2))

      ! Parallelogram Volume.

      case (2)

        Node_Volumes(1) = &
          two * Triangle_Area (R_node(:,1), R_node(:,2), R_node(:,4))
        Node_Volumes(2) = &
          two * Triangle_Area (R_node(:,2), R_node(:,3), R_node(:,1))
        Node_Volumes(3) = &
          two * Triangle_Area (R_node(:,3), R_node(:,4), R_node(:,2))
        Node_Volumes(4) = &
          two * Triangle_Area (R_node(:,4), R_node(:,1), R_node(:,3))

      ! Parallelepiped Volume.

      case (3)

        Node_Volumes(1) = Parallepiped_Volume ( &
                            R_node(:,1), R_node(:,2), R_node(:,5), R_node(:,4))
        Node_Volumes(2) = Parallepiped_Volume ( &
                            R_node(:,2), R_node(:,1), R_node(:,3), R_node(:,6))
        Node_Volumes(3) = Parallepiped_Volume ( &
                            R_node(:,3), R_node(:,2), R_node(:,4), R_node(:,7))
        Node_Volumes(4) = Parallepiped_Volume ( &
                            R_node(:,4), R_node(:,1), R_node(:,8), R_node(:,3))
        Node_Volumes(5) = Parallepiped_Volume ( &
                            R_node(:,5), R_node(:,1), R_node(:,6), R_node(:,8))
        Node_Volumes(6) = Parallepiped_Volume ( &
                            R_node(:,6), R_node(:,2), R_node(:,7), R_node(:,5))
        Node_Volumes(7) = Parallepiped_Volume ( &
                            R_node(:,7), R_node(:,3), R_node(:,8), R_node(:,6))
        Node_Volumes(8) = Parallepiped_Volume ( &
                            R_node(:,8), R_node(:,4), R_node(:,5), R_node(:,7))
      end select

      ! Apply the scalings.

      if (TRIM(Weighting_Scheme) == "Parallel Unscaled" ) then

        ! This is really a null scaling. The segment, parallelogram or 
        ! parallelepiped volumes were calculated on a cell-size basis 
        ! instead of a node-size basis, so they must be scaled down by
        ! the number of nodes per cell.

        Node_Volumes = Node_Volumes / REAL(NLocal_Nodes, real_kind)

      else if (TRIM(Weighting_Scheme) == "Parallel Scaled") then

        ! This scaling assures that the sum of the node volumes is equal 
        ! to the volume of the cell as a whole. This makes sense intuitively.
        ! Note that the true volumes (weight=True_Volume) sum to the cell 
        ! volume, but the unscaled parallel volumes (weight=Parallel_Unscaled)
        ! do not do so (but they are close).

        Node_Volumes = Node_Volumes * Cell_Volume/SUM(Node_Volumes)

      end if
    end if

    ! Verify guarantees.

    if (ANY(Node_Volumes(:) < zero)) then
      write (6,*) 'Error (SO_Node_Volume_Weights): negative node volume ', &
                  'weighting -- cell may be malformed.'
    end if

    return
  end subroutine SO_Node_Volume_Weights
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !
  ! Procedure: Triangle_Area 
  ! Author:    Michael L. Hall, <Hall@LANL.gov>
  ! 
  ! The Triangle_Area procedure returns the area of the triangle formed by 
  ! three node points, represented as 2-D vectors. The coordinates (R) for the 
  ! node points should be given in a counter-clockwise order to assure 
  ! positive area. The area is determined by taking the determinant of this 
  ! matrix:
  ! 
  !    | R1(1) R1(2)  1  |
  !    | R2(1) R2(2)  1  |
  !    | R3(1) R3(2)  1  |
  ! 
  ! See _Introduction to Linear Algebra_, by Strang, p. 228 for more info.

  function Triangle_Area (R1, R2, R3) result(Area)

    ! Input variables.

    ! Three points of the triangle, in counter-clockwise order.
    real(real_kind), dimension(2), intent(in) :: R1, R2, R3

    ! Output variable.

    real(real_kind) :: Area  ! Area of the triangle.

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Calculate the area.

    Area = half * (R1(1)*R2(2) - R2(1)*R1(2) + &
                   R2(1)*R3(2) - R3(1)*R2(2) + &
                   R3(1)*R1(2) - R1(1)*R3(2) )

    ! Verify guarantees.

    if (Area <= zero) then
      write (6,*) 'Error (Triangle_Area): negative or zero area -- ', &
                  'points may be in wrong order.'
    end if

    return
  end function Triangle_Area
  !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end module Augustus

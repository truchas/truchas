!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE NODE_OPERATOR_MODULE
  !======================================================================
  ! Purpose:
  !
  !   Define data structures and procedures associated with the node based
  !   operator for thermo-mechanics
  !=======================================================================
  ! Nxtot is for the boundary kludge
  use kinds, only: r8
  !use parameter_module, only: ndim, nvc, nvf
  use solid_mechanics_mesh, only: ndim, nvc, nvf
  use bc_data_types,    only: BC_MAX_OPERATORS 
  implicit none
  private

  ! Public procedures
  public :: LINEAR_PROP, LINEAR_GRAD

  ! Public types and variables
  public :: CV_Face_Internal, &
            CV_Face_Boundary, &
            CV_Internal,      &
            CV_Boundary,      &
            nipc, nipbf,      &
            cv_init, &
            Nodal_Volume, &
            mech_precond_init
  ! Boundary variables
  public :: nbface, nbnode, nmechbc
  ! Number of valid BC operators - currently x,y,z tractions and x,y,z displacements,
  ! normal displacement, normal traction, free (unconstrained) interface, normally
  ! constrained interface and contact.  Not all of these are fully implemented.
  integer, parameter :: nmechbc = 11
  ! Boundary face and node counters, one for each BC operator
  integer, save,  dimension(BC_MAX_OPERATORS) :: nbface, nbnode
  ! Integration points per cell
  integer, parameter :: nipc = nvc*ndim/2
  ! Integration points per boundary face
  integer, parameter :: nipbf = nvf
  !
  ! Flag to call setup routines only once
  logical, save :: cv_init=.false.
  ! Flag to call precondition setup routines only when needed
  logical, save :: mech_precond_init=.false.
  ! Control volume faces internal to the mesh cells
  type CV_Face_Internal
     ! Normal vector for control volume face
     real(r8), pointer, Dimension(:,:,:) :: Face_Normal
     ! Area of control volume face
     real(r8), pointer, Dimension(:,:) :: Face_Area
     ! IP coordinates for control volume face
     real(r8), pointer, Dimension(:,:,:) :: Face_Coord
     ! Inverse jacobian for gradient
     real(r8), pointer, Dimension(:,:,:,:) :: Face_Ijac
  end type CV_Face_Internal
  ! Control volume faces on boundaries
  type CV_Face_Boundary
     ! Cell that contains these integration points
     integer, pointer, Dimension(:) :: Cell
     ! Nodes associated with these integration points
     integer, pointer, Dimension(:,:) :: Node
     ! BC values for these nodes, only one component
     real(r8), pointer, Dimension(:) :: Node_Value
     ! Surface normal components for these nodes
     real(r8), pointer, Dimension(:,:,:) :: Node_Normal
     ! Number of interfaces (and normal vectors) for these nodes
     integer, pointer, Dimension(:) :: NN_Count
     ! Cell faces that contains these CV faces
     integer, pointer, Dimension(:) :: Face
     ! Node numbers associated with these CV faces
     integer, pointer, Dimension(:,:) :: Face_Node
     ! Normal vectors for control volume faces
     real(r8), pointer, Dimension(:,:,:) :: Face_Normal
     ! Areas of control volume faces
     real(r8), pointer, Dimension(:,:) :: Face_Area
     ! IP coordinates for control volume faces
     real(r8), pointer, Dimension(:,:,:) :: Face_Coord
     ! Inverse jacobian for gradient
     real(r8), pointer, Dimension(:,:,:,:) :: Face_Ijac
  end type CV_Face_Boundary
  ! Create instances of control volume derived types
  type(CV_Face_Internal), save :: CV_Internal
  ! Control volume faces on boundaries
  type(CV_Face_Boundary), save, allocatable, dimension(:) :: CV_Boundary
  ! Nodal volume array for preconditioner
  real(r8), pointer, Dimension(:) :: Nodal_Volume

  ! Components of Linear Interpolation Coefficients used by LINEAR_PROP and LINEAR_GRAD
  real(r8), parameter :: LC1(3,8) = reshape( &
      [ 0.0_r8,  1.0_r8,  1.0_r8, &
        0.0_r8,  0.0_r8,  1.0_r8, &
        1.0_r8,  0.0_r8,  1.0_r8, &
        1.0_r8,  1.0_r8,  1.0_r8, &
        0.0_r8,  1.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8,  0.0_r8, &
        1.0_r8,  0.0_r8,  0.0_r8, &
        1.0_r8,  1.0_r8,  0.0_r8 ], shape=[3,8])

  real(r8), parameter :: LC2(3,8) = reshape( &
      [ 1.0_r8, -1.0_r8, -1.0_r8, &
        1.0_r8,  1.0_r8, -1.0_r8, &
       -1.0_r8,  1.0_r8, -1.0_r8, &
       -1.0_r8, -1.0_r8, -1.0_r8, &
        1.0_r8, -1.0_r8,  1.0_r8, &
        1.0_r8,  1.0_r8,  1.0_r8, &
       -1.0_r8,  1.0_r8,  1.0_r8, &
       -1.0_r8, -1.0_r8,  1.0_r8 ], shape=[3,8])

CONTAINS

  SUBROUTINE LINEAR_PROP (num, Xi, Vrtx, Prop, id)
    !=======================================================================
    ! Purpose(s):
    !   Property evaluation using 2-D or 3-D linear interpolation at
    !   the logical coordinate vector location.
    !
    !               nvc
    !               ---
    !               \
    !   Prop_p(:) =  |  Coef_v(:) * Vrtx(v,:)
    !               /
    !               ---
    !               v=1
    !
    !               ndim
    !               ___
    !               | |
    !   Coef_v(:) = | | [LC1(n,v) + LC2(n,v)*Xi(n,:)]
    !               | |
    !               n=1
    !
    !=======================================================================

    ! Argument List
    integer, optional, intent(IN)  :: id
    integer,           intent(IN)  :: num

    real(r8), dimension(num),      intent(OUT) :: Prop
    real(r8), dimension(ndim,num), intent(IN)  :: Xi
    real(r8), dimension(nvc,num),  intent(IN)  :: Vrtx

    ! Local Variables
    integer :: e, s, n, v

    real(r8), dimension(num) :: Coef

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Extent Bounds
    if (present(id)) then
       s = id ! Vector Element
       e = id
    else
       s = 1  ! Entire Vector
       e = num
    end if

    ! Initialize Property
    Prop(s:e) = 0.0_r8

    ! Evaluate Property
    Vrtx_Loop : do v = 1, nvc
       ! Initialize Coefficient
       Coef(s:e) = 1.0_r8

       ! Accumulate Coefficient
       Ndim_Loop : do n = 1, ndim
          Coef(s:e) = Coef(s:e) * (LC1(n,v) + LC2(n,v)*Xi(n,s:e))
       end do Ndim_Loop

       ! Accumulate Property
       Prop(s:e) = Prop(s:e) + Coef(s:e)*Vrtx(v,s:e)
    end do Vrtx_Loop

  END SUBROUTINE LINEAR_PROP

  SUBROUTINE LINEAR_GRAD (num, Xi, Vrtx, Grad, id)
    !=======================================================================
    ! Purpose(s):
    !   Gradient evaluation using 2-D or 3-D linear interpolation at
    !   the logical coordinate vector location.
    !
    !               nvc
    !               ---
    !               \
    !   Grad(n,:) =  |  Coef_v(:) * Vrtx(v,:)
    !               /
    !               ---
    !               v=1
    !
    !                         ndim-1
    !                          ___
    !                          | |
    !   Coef_v(:) = LC2(n,v) * | | [LC1(nlc,v) + LC2(nlc,v)*Xi(nlc,:)]
    !                          | |
    !                          nn=1
    !
    !
    !   The index 'nlc' is never equivalent to n, and is derived from
    !   the Tensor matrix, nlc = Tensor(nn+1,n).
    !=======================================================================
    use tensor_module,    only: Tensor

    ! Argument List
    integer, intent(IN), optional :: id
    integer, intent(IN)           :: num

    real(r8), dimension(ndim,num), intent(IN)  :: Xi
    real(r8), dimension(ndim,num), intent(OUT) :: Grad
    real(r8), dimension(nvc,num),  intent(IN)  :: Vrtx

    ! Local Variables
    integer :: e, s, v, n
    integer :: nn, nlc

    real(r8), dimension(num) :: Coef

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Extent Bounds
    if (present(id)) then
       s = id ! Vector Element
       e = id
    else
       s = 1  ! Entire Vector
       e = num
    end if

    ! Initialize Gradient
    Grad(:,s:e) = 0.0_r8

    ! Evaluate Gradient
    Ndim_Loop : do n = 1, ndim
       Vrtx_Loop : do v = 1, nvc
          ! Initialize Coefficient
          Coef(s:e) = 1.0_r8

          ! Coefficient Factor Terms
          do nn = 1, ndim-1
             nlc = Tensor(nn+1,n) ! Logical Coordinate
             Coef(s:e) = Coef(s:e) * (LC1(nlc,v) + LC2(nlc,v)*Xi(nlc,s:e))
          end do

          ! Coefficient Sign
          Coef(s:e) = LC2(n,v)*Coef(s:e)

          ! Accumulate Gradient
          Grad(n,s:e) = Grad(n,s:e) + Coef(s:e)*Vrtx(v,s:e)
       end do Vrtx_Loop
    end do Ndim_Loop

  END SUBROUTINE LINEAR_GRAD

END MODULE NODE_OPERATOR_MODULE


MODULE LINEAR_MODULE
  !=======================================================================
  ! Purpose(s):
  !   Define various procedures for evaluating a functions property
  !   and gradient using 2-D or 3-D linear interpolation.
  !
  ! Contains: LINEAR_GRAD        <- Interface(s)
  !           LINEAR_PROP
  !
  !           LINEAR_PROP_FACE   <- Procedure(s)
  !           LINEAR_PROP_VECTOR
  !           LINEAR_GRAD_VECTOR
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: ndim, nvc
  implicit none
  private

  ! Public Subroutines
  public :: LINEAR_GRAD
  public :: LINEAR_PROP

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Interface Procedures
  INTERFACE LINEAR_PROP ! Linear Property Evaluation
     MODULE PROCEDURE LINEAR_PROP_FACE
     MODULE PROCEDURE LINEAR_PROP_FACE_INT
     MODULE PROCEDURE LINEAR_PROP_VECTOR
  END INTERFACE

  INTERFACE LINEAR_GRAD ! Linear Gradient Evaluation
     MODULE PROCEDURE LINEAR_GRAD_VECTOR
  END INTERFACE

  ! Components of Linear Interpolation Coefficients
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

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE LINEAR_PROP_FACE (f, Vrtx, Prop)
    !=======================================================================
    ! Purpose(s):
    !   Property evaluation using 2-D or 3-D linear interpolation at
    !   the cell-face centroid vector location.
    !
    !               nvc
    !               ---
    !               \
    !   Prop_p(:) =  |  Coef_v(:) * Vrtx_v(:)
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
    use mesh_module,      only: Cell
    use parameter_module, only: ncells, ndim, nvc

    ! Argument List
    integer, intent(IN) :: f ! Cell-Face Number

    real(r8), dimension(ncells),     intent(OUT) :: Prop
    real(r8), dimension(nvc,ncells), intent(IN)  :: Vrtx

    ! Local Variables
    integer :: n, v

    real(r8), dimension(ncells) :: Coef

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize Property
    Prop(:) = 0.0_r8

    ! Evaluate Property
    Vrtx_Loop : do v = 1, nvc
       ! Initialize Coefficient
       Coef(:) = 1.0_r8

       ! Accumulate Coefficient
       Ndim_Loop : do n = 1, ndim
          Coef(:) = Coef(:) * (LC1(n,v) + LC2(n,v)*Cell(:)%Face_Centroid_L(n,f))
       end do Ndim_Loop

       ! Accumulate Property
       Prop(:) = Prop(:) + Coef(:)*Vrtx(v,:)
    end do Vrtx_Loop

  END SUBROUTINE LINEAR_PROP_FACE

  SUBROUTINE LINEAR_PROP_FACE_INT (f, Vrtx, Prop)
    !=======================================================================
    ! Purpose(s):
    !   Property evaluation using 2-D or 3-D linear interpolation at
    !   the cell-face centroid vector location.
    !
    !               nvc
    !               ---
    !               \
    !   Prop_p(:) =  |  Coef_v(:) * Vrtx_v(:)
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
    use mesh_module,      only: Cell
    use parameter_module, only: ncells, ndim, nvc

    ! Argument List
    integer, intent(IN) :: f ! Cell-Face Number

    real(r8), dimension(ncells),     intent(OUT) :: Prop
    integer,  dimension(nvc,ncells), intent(IN)  :: Vrtx

    ! Local Variables
    integer :: n, v

    real(r8), dimension(ncells) :: Coef

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize Property
    Prop(:) = 0.0_r8

    ! Evaluate Property
    Vrtx_Loop : do v = 1, nvc
       ! Initialize Coefficient
       Coef(:) = 1.0_r8

       ! Accumulate Coefficient
       Ndim_Loop : do n = 1, ndim
          Coef(:) = Coef(:) * (LC1(n,v) + LC2(n,v)*Cell(:)%Face_Centroid_L(n,f))
       end do Ndim_Loop

       ! Accumulate Property
       Prop(:) = Prop(:) + Coef(:)*Vrtx(v,:)
    end do Vrtx_Loop

  END SUBROUTINE LINEAR_PROP_FACE_INT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE LINEAR_PROP_VECTOR (num, Xi, Vrtx, Prop, id)
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
    use parameter_module, only: ndim, nvc

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

  END SUBROUTINE LINEAR_PROP_VECTOR

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE LINEAR_GRAD_VECTOR (num, Xi, Vrtx, Grad, id)
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
    use parameter_module, only: ndim, nvc
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

  END SUBROUTINE LINEAR_GRAD_VECTOR

END MODULE LINEAR_MODULE

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
  !           LINEAR_COEF
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kind_module,      only: real_kind
  use parameter_module, only: ndim, nvc

  implicit none

  ! Pivate Module
  private

  ! Public Subroutines
  public :: LINEAR_COEF
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
  real(KIND = real_kind), dimension(ndim,nvc) :: LC1, LC2

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
    use constants_module, only: one, zero
    use kind_module,      only: int_kind, real_kind
    use mesh_module,      only: Cell
    use parameter_module, only: ncells, ndim, nvc

    implicit none

    ! Argument List
    integer(KIND = int_kind), intent(IN) :: f ! Cell-Face Number

    real(KIND = real_kind), dimension(ncells),     intent(OUT) :: Prop
    real(KIND = real_kind), dimension(nvc,ncells), intent(IN)  :: Vrtx

    ! Local Variables
    integer(KIND = int_kind) :: n, v

    real(KIND = real_kind), dimension(ncells) :: Coef

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize Property
    Prop(:) = zero

    ! Evaluate Property
    Vrtx_Loop : do v = 1, nvc
       ! Initialize Coefficient
       Coef(:) = one

       ! Accumulate Coefficient
       Ndim_Loop : do n = 1, ndim
          Coef(:) = Coef(:) * (LC1(n,v) + LC2(n,v)*Cell(:)%Face_Centroid_L(n,f))
       end do Ndim_Loop

       ! Accumulate Property
       Prop(:) = Prop(:) + Coef(:)*Vrtx(v,:)
    end do Vrtx_Loop

    return

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
    use constants_module, only: one, zero
    use kind_module,      only: int_kind, real_kind
    use mesh_module,      only: Cell
    use parameter_module, only: ncells, ndim, nvc

    implicit none

    ! Argument List
    integer(KIND = int_kind), intent(IN) :: f ! Cell-Face Number

    real(KIND = real_kind),   dimension(ncells),     intent(OUT) :: Prop
    integer(KIND = int_kind), dimension(nvc,ncells), intent(IN)  :: Vrtx

    ! Local Variables
    integer(KIND = int_kind) :: n, v

    real(KIND = real_kind), dimension(ncells) :: Coef

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize Property
    Prop(:) = zero

    ! Evaluate Property
    Vrtx_Loop : do v = 1, nvc
       ! Initialize Coefficient
       Coef(:) = one

       ! Accumulate Coefficient
       Ndim_Loop : do n = 1, ndim
          Coef(:) = Coef(:) * (LC1(n,v) + LC2(n,v)*Cell(:)%Face_Centroid_L(n,f))
       end do Ndim_Loop

       ! Accumulate Property
       Prop(:) = Prop(:) + Coef(:)*Vrtx(v,:)
    end do Vrtx_Loop

    return

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
    use constants_module, only: zero, one
    use kind_module,      only: int_kind, real_kind
    use parameter_module, only: ndim, nvc

    implicit none

    ! Argument List
    integer(KIND = int_kind), optional, intent(IN)  :: id
    integer(KIND = int_kind),           intent(IN)  :: num

    real(KIND = real_kind), dimension(num),      intent(OUT) :: Prop
    real(KIND = real_kind), dimension(ndim,num), intent(IN)  :: Xi
    real(KIND = real_kind), dimension(nvc,num),  intent(IN)  :: Vrtx

    ! Local Variables
    integer(KIND = int_kind) :: e, s, n, v

    real(KIND = real_kind), dimension(num) :: Coef

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
    Prop(s:e) = zero

    ! Evaluate Property
    Vrtx_Loop : do v = 1, nvc
       ! Initialize Coefficient
       Coef(s:e) = one

       ! Accumulate Coefficient
       Ndim_Loop : do n = 1, ndim
          Coef(s:e) = Coef(s:e) * (LC1(n,v) + LC2(n,v)*Xi(n,s:e))
       end do Ndim_Loop

       ! Accumulate Property
       Prop(s:e) = Prop(s:e) + Coef(s:e)*Vrtx(v,s:e)
    end do Vrtx_Loop

    return

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
    use constants_module, only: one, zero
    use kind_module,      only: int_kind, real_kind
    use parameter_module, only: ndim, nvc
    use tensor_module,    only: Tensor

    implicit none

    ! Argument List
    integer(KIND = int_kind), intent(IN), optional :: id
    integer(KIND = int_kind), intent(IN)           :: num

    real(KIND = real_kind), dimension(ndim,num), intent(IN)  :: Xi
    real(KIND = real_kind), dimension(ndim,num), intent(OUT) :: Grad
    real(KIND = real_kind), dimension(nvc,num),  intent(IN)  :: Vrtx

    ! Local Variables
    integer(KIND = int_kind) :: e, s, v, n
    integer(KIND = int_kind) :: nn, nlc

    real(KIND = real_kind), dimension(num) :: Coef

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
    Grad(:,s:e) = zero

    ! Evaluate Gradient
    Ndim_Loop : do n = 1, ndim
       Vrtx_Loop : do v = 1, nvc
          ! Initialize Coefficient
          Coef(s:e) = one

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

    return

  END SUBROUTINE LINEAR_GRAD_VECTOR

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  SUBROUTINE LINEAR_COEF ()
    !=======================================================================
    ! Purpose(s):
    !   Evaluate components of linear interpolation coefficients.
    !
    !               ndim
    !               ___
    !               | |
    !   Coef_v(:) = | | [LC1(n,v) + LC2(n,v)*Xi(n,:)]
    !               | |
    !               n=1
    !
    !   LC1: v\n|  1  2  3       LC2: v\n|  1  2  3
    !        ___|___________          ___|___________
    !         1 |  0  1  1             1 |  1 -1 -1
    !         2 |  0  0  1             2 |  1  1 -1
    !         3 |  1  0  1             3 | -1  1 -1
    !         4 |  1  1  1             4 | -1 -1 -1
    !         5 |  0  1  0             5 |  1 -1  1
    !         6 |  0  0  0             6 |  1  1  1
    !         7 |  1  0  0             7 | -1  1  1
    !         8 |  1  1  0             8 | -1 -1  1
    !
    !=======================================================================
    use constants_module, only: zero, one
    use kind_module,      only: int_kind
    use parameter_module, only: ndim, nvc

    implicit none

    ! Argument List

    ! Local Variables
    integer(KIND = int_kind) :: n, v

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Evaluate Components
    Vrtx_Loop : do v = 1, nvc
       select case (v)
       case (1) ! Vertex One
          do n = 1, ndim
             select case (n)
             case (1) ! 1-D Coefficients
                LC1(n,v) = + zero
                LC2(n,v) = + one
             case (2) ! 2-D Coefficients
                LC1(n,v) = + one
                LC2(n,v) = - one
             case (3) ! 3-D Coefficients
                LC1(n,v) = + one
                LC2(n,v) = - one
             end select
          end do
       case (2) ! Vertex Two
          do n = 1, ndim
             select case (n)
             case (1) ! 1-D Coefficients
                LC1(n,v) = + zero
                LC2(n,v) = + one
             case (2) ! 2-D Coefficients
                LC1(n,v) = + zero
                LC2(n,v) = + one
             case (3) ! 3-D Coefficients
                LC1(n,v) = + one
                LC2(n,v) = - one
             end select
          end do
       case (3) ! Vertex Three
          do n = 1, ndim
             select case (n)
             case (1) ! 1-D Coefficients
                LC1(n,v) = + one
                LC2(n,v) = - one
             case (2) ! 2-D Coefficients
                LC1(n,v) = + zero
                LC2(n,v) = + one
             case (3) ! 3-D Coefficients
                LC1(n,v) = + one
                LC2(n,v) = - one
             end select
          end do
       case (4) ! Vertex Four
          do n = 1, ndim
             select case (n)
             case (1) ! 1-D Coefficients
                LC1(n,v) = + one
                LC2(n,v) = - one
             case (2) ! 2-D Coefficients
                LC1(n,v) = + one
                LC2(n,v) = - one
             case (3) ! 3-D Coefficients
                LC1(n,v) = + one
                LC2(n,v) = - one
             end select
          end do
       case (5) ! Vertex Five
          do n = 1, ndim
             select case (n)
             case (1) ! 1-D Coefficients
                LC1(n,v) = + zero
                LC2(n,v) = + one
             case (2) ! 2-D Coefficients
                LC1(n,v) = + one
                LC2(n,v) = - one
             case (3) ! 3-D Coefficients
                LC1(n,v) = + zero
                LC2(n,v) = + one
             end select
          end do
       case (6) ! Vertex Six
          do n = 1, ndim
             select case (n)
             case (1) ! 1-D Coefficients
                LC1(n,v) = + zero
                LC2(n,v) = + one
             case (2) ! 2-D Coefficients
                LC1(n,v) = + zero
                LC2(n,v) = + one
             case (3) ! 3-D Coefficients
                LC1(n,v) = + zero
                LC2(n,v) = + one
             end select
          end do
       case (7) ! Vertex Seven
          do n = 1, ndim
             select case (n)
             case (1) ! 1-D Coefficients
                LC1(n,v) = + one
                LC2(n,v) = - one
             case (2) ! 2-D Coefficients
                LC1(n,v) = + zero
                LC2(n,v) = + one
             case (3) ! 3-D Coefficients
                LC1(n,v) = + zero
                LC2(n,v) = + one
             end select
          end do
       case (8) ! Vertex Eight
          do n = 1, ndim
             select case (n)
             case (1) ! 1-D Coefficients
                LC1(n,v) = + one
                LC2(n,v) = - one
             case (2) ! 2-D Coefficients
                LC1(n,v) = + one
                LC2(n,v) = - one
             case (3) ! 3-D Coefficients
                LC1(n,v) = + zero
                LC2(n,v) = + one
             end select
          end do
       end select
    end do Vrtx_Loop

    return

  END SUBROUTINE LINEAR_COEF

END MODULE LINEAR_MODULE

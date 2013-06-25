! -*- Mode: f90 -*-

module ubik_vector_type_module
  use ubik_types_module

  implicit none
  save

  type vector_type
!     private
     real(real_type), dimension(:), pointer :: values => NULL()
     real(real_type), dimension(:), pointer :: aux1 => NULL()
     real(real_type), dimension(:), pointer :: aux2 => NULL()
     real(real_type), dimension(:), pointer :: aux3 => NULL()
     integer(int_type) :: status
  end type vector_type

  interface create
     module procedure create_only
     module procedure create_and_size
  end interface

  interface destroy
     module procedure destroy_vec
  end interface

  interface values
     module procedure values_all
     module procedure value_single
  end interface

  interface spew
     module procedure Write1Vector
     module procedure Write2Vectors
     module procedure Write3Vectors
  end interface

  interface set_status_alloc_failure
     module procedure set_status_alloc_failure_vec
  end interface

  interface alloc_failure
     module procedure alloc_failure_vec
  end interface

  interface SIZE
     module procedure SIZE_vec
  end interface

  interface Norm
     module procedure Norm_vec
  end interface

  interface assignment(=)
     module procedure set_values_int
     module procedure set_values_real
     module procedure set_values_real_vector
  end interface

#ifndef dec
  interface operator(-)
     module procedure vector_minus_vector_type
     module procedure vector_type_minus_vector
  end interface

  interface operator(+)
     module procedure vector_type_plus_vector_type
  end interface
#endif

  private
  public :: &
       vector_type, &
       create, &
       destroy, &
       assignment(=), &
       set_values_ptr, &
       set_aux1_ptr, &
       values, &
       values_ptr, &
       aux1_ptr, &
       aux1_elem, &
       spew, &
       set_status_alloc_failure, &
       alloc_failure, &
       SIZE, &
       Norm

#ifdef dec
  public :: &
       vector_type_plus_vector_type, &
       vector_minus_vector_type, &
       vector_type_minus_vector
#else
  public :: &
       operator(+), &
       operator(-)
#endif

contains

  subroutine set_values_int (vec, i)
    integer(int_type), intent(in) :: i
    type(vector_type), intent(inout) :: vec

    vec%values = i

    return
  end subroutine set_values_int

  subroutine set_values_real (vec, r)
    real(real_type), intent(in) :: r
    type(vector_type), intent(inout) :: vec

    vec%values = r

    return
  end subroutine set_values_real

  subroutine set_values_real_vector (vec, v)
    real(real_type), dimension(:), intent(in) :: v
    type(vector_type), intent(inout) :: vec

    vec%values = v

    return
  end subroutine set_values_real_vector

  subroutine set_values_ptr (vec, v)
    real(real_type), dimension(:), target, intent(in) :: v
    type(vector_type), intent(inout) :: vec

    vec%values => v

    return
  end subroutine set_values_ptr

  function values_ptr(vec) result(ptr)
    type(vector_type), intent(in) :: vec
    real(real_type), dimension(:), pointer :: ptr

    ptr => vec%values

    return
  end function values_ptr

  function values_all (vec)
    type(vector_type), intent(in) :: vec
    real(real_type), dimension(SIZE(vec%values)) :: values_all

    values_all = vec%values

    return
  end function values_all

  function value_single (vec, i)
    integer(int_type), intent(in) :: i
    type(vector_type), intent(in) :: vec
    real(real_type) :: value_single

    value_single = vec%values(i)

    return
  end function value_single

  subroutine set_aux1_ptr (vec, ptr)
    real(real_type), dimension(:), pointer :: ptr
    type(vector_type), intent(inout) :: vec

    vec%aux1 => ptr

    return
  end subroutine set_aux1_ptr

  function aux1_ptr(vec) result(ptr)
    type(vector_type), intent(in) :: vec
    real(real_type), dimension(:), pointer :: ptr

    ptr => vec%aux1

    return
  end function aux1_ptr

  function aux1_elem (vec, i)
    integer(int_type), intent(in) :: i
    type(vector_type), intent(in) :: vec
    real(real_type) :: aux1_elem

    aux1_elem = vec%aux1(i)

    return
  end function aux1_elem

  !*******************************************************************
  !
  ! constructors and destructors
  !
  !**********************************************************************
  subroutine create_only (vec, overlap_only)
    logical(logical_type), intent(in), optional :: overlap_only
    type(vector_type), intent(inout) :: vec

    if (.not.PRESENT(overlap_only)) then
       nullify (vec%values)
    end if
    nullify (vec%aux1)
    nullify (vec%aux2)
    nullify (vec%aux3)

    return
  end subroutine create_only

  subroutine create_and_size (vec, n, overlap_only)
    logical(logical_type), intent(in), optional :: overlap_only
    integer(int_type), intent(in) :: n
    type(vector_type), intent(inout) :: vec
    integer(int_type) :: status

    if (.not.PRESENT(overlap_only)) then
       allocate (vec%values(1:n), stat=status)
       if (status /= 0) then
          call set_status_alloc_failure(vec)
       end if
       return
    end if

    allocate (vec%aux1(1:n), vec%aux2(1:n), vec%aux3(1:n), stat=status)
    if (status /= 0) then
       call set_status_alloc_failure(vec)
    end if
 
    return
  end subroutine create_and_size

  subroutine destroy_vec (vec, overlap_only)
    logical(logical_type), intent(in), optional :: overlap_only
    type(vector_type), intent(inout) :: vec
    integer(int_type) :: status

    if (.not.PRESENT(overlap_only)) then
       if (associated(vec%values)) then
          deallocate (vec%values, stat=status)
          if (status /= 0) then
             call set_status_alloc_failure(vec)
          end if
       end if
       return
    end if
    if (associated(vec%aux1)) then
       deallocate (vec%aux1, stat=status)
       if (status /= 0) then
          call set_status_alloc_failure(vec)
       end if
    end if
    if (associated(vec%aux2)) then
       deallocate (vec%aux2, stat=status)
       if (status /= 0) then
          call set_status_alloc_failure(vec)
       end if
    end if
    if (associated(vec%aux3)) then
       deallocate (vec%aux3, stat=status)
       if (status /= 0) then
          call set_status_alloc_failure(vec)
       end if
    end if

    return
  end subroutine destroy_vec

  subroutine Write1Vector (lu, v, title)
    use Ubik_Vector_module

    character(len=*), intent(in) :: title
    integer(int_type), intent(in) :: lu
    type(vector_type), intent(in) :: v

    call spew (lu, v%values, title)

    return
  end subroutine Write1Vector

  subroutine Write2Vectors (lu, u, v, title)
    use Ubik_Vector_module

    character(len=*), intent(in) :: title
    integer(int_type), intent(in) :: lu
    type(vector_type), intent(in) :: u, v

    call spew (lu, u%values, v%values, title)

    return
  end subroutine Write2Vectors

  subroutine Write3Vectors (lu, u, v, w, title)
    use Ubik_Vector_module

    character(len=*), intent(in) :: title
    integer(int_type), intent(in) :: lu
    type(vector_type), intent(in) :: u, v, w

    call spew (lu, u%values, v%values, w%values, title)

    return
  end subroutine Write3Vectors

  !*****************************************************************
  !
  ! status setting and accessing
  !
  !*****************************************************************

  subroutine set_status_alloc_failure_vec(vec)
    use ubik_constants_module

    type(vector_type), intent(inout) :: vec

    vec%status = alloc_failure__
    return
  end subroutine set_status_alloc_failure_vec

  function alloc_failure_vec(vec)
    use ubik_constants_module

    logical(logical_type) :: alloc_failure_vec
    type(vector_type), intent(in) :: vec

    alloc_failure_vec = (vec%status == alloc_failure__)
    return
  end function alloc_failure_vec

  function Norm_vec(inorm, vec)
    use Ubik_Vector_module

    real(real_type) :: Norm_vec
    integer(int_type), intent(in) :: inorm
    type(vector_type), intent(in) :: vec

    Norm_vec = Norm(inorm, vec%values)

    return
  end function Norm_vec

  function SIZE_vec (vec)
    use Ubik_Vector_module

    integer(int_type) :: SIZE_vec
    type(vector_type), intent(in) :: vec

    SIZE_vec = SIZE(vec%values)

    return
  end function SIZE_vec

  function vector_type_plus_vector_type (vec1, vec2)
    type(vector_type), intent(in) :: vec1, vec2
    real(real_type), dimension(SIZE(vec1%values)) :: vector_type_plus_vector_type

    vector_type_plus_vector_type = vec1%values + vec2%values

    return
  end function vector_type_plus_vector_type

  function vector_minus_vector_type (v, vec)
    real(real_type), dimension(:), intent(in) :: v
    type(vector_type), intent(in) :: vec
    real(real_type), dimension(SIZE(vec%values)) :: vector_minus_vector_type

    vector_minus_vector_type = v - vec%values

    return
  end function vector_minus_vector_type

  function vector_type_minus_vector (vec, v)
    real(real_type), dimension(:), intent(in) :: v
    type(vector_type), intent(in) :: vec
    real(real_type), dimension(SIZE(vec%values)) :: vector_type_minus_vector

    vector_type_minus_vector = vec%values - v

    return
  end function vector_type_minus_vector

end module ubik_vector_type_module

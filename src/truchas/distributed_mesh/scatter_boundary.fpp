# /* include file for index_partitioning.F90 */

#ifdef _OP_
#undef _OP_
#endif

#ifdef _SUM_REDUCTION_
#  undef _SUM_REDUCTION_
#  define _OP_(a,b)   a + b
#  ifdef _INTEGER_DATA_
#    undef _INTEGER_DATA_
#    define _TYPE_  integer
#    define _ID_(a) 0
#    define _PROC1_ scatter_boundary_sum_I1
#    define _PROC2_ scatter_boundary_sum_XI1
#    define _AUX_   scatter_boundary_sum_I_aux
#  elif defined(_SINGLE_DATA_)
#    undef _SINGLE_DATA_
#    define _TYPE_  real
#    define _ID_(a) 0.0
#    define _PROC1_ scatter_boundary_sum_S1
#    define _PROC2_ scatter_boundary_sum_XS1
#    define _AUX_   scatter_boundary_sum_S_aux
#  elif defined(_DOUBLE_DATA_)
#    undef _DOUBLE_DATA_
#    define _TYPE_  real(kind(1.0d0))
#    define _ID_(a) 0.0d0
#    define _PROC1_ scatter_boundary_sum_D1
#    define _PROC2_ scatter_boundary_sum_XD1
#    define _AUX_   scatter_boundary_sum_D_aux
#  endif
#endif

#ifdef _MIN_REDUCTION_
#  undef _MIN_REDUCTION_
#  define _OP_(a,b)   min(a, b)
#  define _ID_(a)     huge(a)
#  ifdef _INTEGER_DATA_
#    undef _INTEGER_DATA_
#    define _TYPE_  integer
#    define _PROC1_ scatter_boundary_min_I1
#    define _PROC2_ scatter_boundary_min_XI1
#    define _AUX_   scatter_boundary_min_I_aux
#  elif defined(_SINGLE_DATA_)
#    undef _SINGLE_DATA_
#    define _TYPE_  real
#    define _PROC1_ scatter_boundary_min_S1
#    define _PROC2_ scatter_boundary_min_XS1
#    define _AUX_   scatter_boundary_min_S_aux
#  elif defined(_DOUBLE_DATA_)
#    undef _DOUBLE_DATA_
#    define _TYPE_  real(kind(1.0d0))
#    define _PROC1_ scatter_boundary_min_D1
#    define _PROC2_ scatter_boundary_min_XD1
#    define _AUX_   scatter_boundary_min_D_aux
#  endif
#endif

#ifdef _MAX_REDUCTION_
#  undef _MAX_REDUCTION_
#  define _OP_(a,b)   max(a, b)
#  define _ID_(a)     -huge(a)
#  ifdef _INTEGER_DATA_
#    undef _INTEGER_DATA_
#    define _TYPE_  integer
#    define _PROC1_ scatter_boundary_max_I1
#    define _PROC2_ scatter_boundary_max_XI1
#    define _AUX_   scatter_boundary_max_I_aux
#  elif defined(_SINGLE_DATA_)
#    undef _SINGLE_DATA_
#    define _TYPE_  real
#    define _PROC1_ scatter_boundary_max_S1
#    define _PROC2_ scatter_boundary_max_XS1
#    define _AUX_   scatter_boundary_max_S_aux
#  elif defined(_DOUBLE_DATA_)
#    undef _DOUBLE_DATA_
#    define _TYPE_  real(kind(1.0d0))
#    define _PROC1_ scatter_boundary_max_D1
#    define _PROC2_ scatter_boundary_max_XD1
#    define _AUX_   scatter_boundary_max_D_aux
#  endif
#endif

#ifdef _OR_REDUCTION_
#undef _OR_REDUCTION_
#define _OP_(a,b)   a .or. b
#define _ID_(a)     .false.
#define _TYPE_      logical
#define _PROC1_ scatter_boundary_or_L1
#define _PROC2_ scatter_boundary_or_XL1
#define _AUX_   scatter_boundary_or_aux
#endif

#ifdef _AND_REDUCTION_
#undef _AND_REDUCTION_
#define _OP_(a,b)   a .and. b
#define _ID_(a)     .true.
#define _TYPE_      logical
#define _PROC1_ scatter_boundary_and_L1
#define _PROC2_ scatter_boundary_and_XL1
#define _AUX_   scatter_boundary_and_aux
#endif

#ifndef _OP_
#error "One of {_SUM_, _MIN_, _MAX_, _AND_, _OR_}_REDUCTION_ must be defined"
#endif

  subroutine _PROC1_ (this, onP_data, offP_data)

    type(ip_desc), intent(in) :: this    ! partition descriptor
    _TYPE_, intent(inout) :: onP_data(:)   ! on-process data
    _TYPE_, intent(in)    :: offP_data(:)  ! off-process data

    ASSERT( defined(this) )
    ASSERT( associated(this%offP_index) )
    ASSERT( size(onP_data) == this%onP_size_ )
    ASSERT( size(offP_data) == this%offP_size_ )

    call _AUX_ (this, onP_data, offP_data)

  end subroutine _PROC1_

  subroutine _PROC2_ (this, local_data)

    type(ip_desc), intent(in) :: this   ! partition descriptor
    _TYPE_, intent(inout) :: local_data(:)

    ASSERT( defined(this) )
    ASSERT( associated(this%offP_index) )
    ASSERT( size(local_data) == this%local_size_ )

    call _AUX_ (this, local_data(:this%onP_size_), local_data(this%onP_size_+1:))

  end subroutine _PROC2_

  subroutine _AUX_ (this, onP_data, offP_data)

    type(ip_desc), intent(in) :: this    ! partition descriptor
    _TYPE_, intent(inout) :: onP_data(:)   ! on-process data
    _TYPE_, intent(in)    :: offP_data(:)  ! off-process data

    integer :: j
    _TYPE_ :: dup_buffer(size(this%dup_index))
    _TYPE_, allocatable :: sup_buffer(:)

    if (associated(this%sup_index)) then ! SUP_BUFFER /= OFFP_DATA

      !! Scatter the off-process data into the supplement buffer.
      allocate(sup_buffer(pgslib_size_of_sup(this%trace)))
      sup_buffer = _ID_(sup_buffer)
      do j = 1, size(this%sup_index)
        sup_buffer(this%sup_index(j)) = _OP_(sup_buffer(this%sup_index(j)), offP_data(j))
      end do

      !! Global communication: SUP_BUFFER -> DUP_BUFFER.
      dup_buffer = pgslib_scatter_buffer(sup_buffer, this%trace)
      deallocate(sup_buffer)

    else  ! SUP_BUFFER == OFFP_DATA

      !! Global communication: SUP_BUFFER (= OFFP_DATA) -> DUP_BUFFER.
      dup_buffer = pgslib_scatter_buffer(offP_data, this%trace)

    end if

    !! Accumulate the duplicate buffer values to the on-process data.
    do j = 1, size(this%dup_index)
      onP_data(this%dup_index(j)) = _OP_(onP_data(this%dup_index(j)), dup_buffer(j))
    end do

  end subroutine _AUX_

#undef _OP_
#undef _ID_
#undef _TYPE_
#undef _PROC1_
#undef _PROC2_
#undef _AUX_

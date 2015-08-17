# /* include file for index_partitioning.F90 */

#ifdef _TYPE_
#undef _TYPE_
#endif

#ifdef _LOGICAL_DATA_
#undef _LOGICAL_DATA_
#define _TYPE_ logical
#define _PROC_GB_A1_ gather_boundary_L1
#define _PROC_GB_A2_ gather_boundary_L2
#define _PROC_GB_A3_ gather_boundary_L3
#define _PROC_GB_B1_ gather_boundary_XL1
#define _PROC_GB_B2_ gather_boundary_XL2
#define _PROC_GB_B3_ gather_boundary_XL3
#define _AUX_   gather_boundary_L_aux
#define _PROC_BIC1_ boundary_is_current_L1
#define _PROC_BIC2_ boundary_is_current_L2
#define _PROC_BIC3_ boundary_is_current_L3
#define _EQUAL_ .eqv.
#endif

#ifdef _INTEGER_DATA_
#undef _INTEGER_DATA_
#define _TYPE_ integer
#define _PROC_GB_A1_ gather_boundary_I1
#define _PROC_GB_A2_ gather_boundary_I2
#define _PROC_GB_A3_ gather_boundary_I3
#define _PROC_GB_B1_ gather_boundary_XI1
#define _PROC_GB_B2_ gather_boundary_XI2
#define _PROC_GB_B3_ gather_boundary_XI3
#define _AUX_   gather_boundary_I_aux
#define _PROC_BIC1_ boundary_is_current_I1
#define _PROC_BIC2_ boundary_is_current_I2
#define _PROC_BIC3_ boundary_is_current_I3
#define _EQUAL_ ==
#endif

#ifdef _SINGLE_DATA_
#undef _SINGLE_DATA_
#define _TYPE_ real
#define _PROC_GB_A1_ gather_boundary_S1
#define _PROC_GB_A2_ gather_boundary_S2
#define _PROC_GB_A3_ gather_boundary_S3
#define _PROC_GB_B1_ gather_boundary_XS1
#define _PROC_GB_B2_ gather_boundary_XS2
#define _PROC_GB_B3_ gather_boundary_XS3
#define _AUX_   gather_boundary_S_aux
#define _PROC_BIC1_ boundary_is_current_S1
#define _PROC_BIC2_ boundary_is_current_S2
#define _PROC_BIC3_ boundary_is_current_S3
#define _EQUAL_ ==
#endif

#ifdef _DOUBLE_DATA_
#undef _DOUBLE_DATA_
#define _TYPE_ real(kind(1.0d0))
#define _PROC_GB_A1_ gather_boundary_D1
#define _PROC_GB_A2_ gather_boundary_D2
#define _PROC_GB_A3_ gather_boundary_D3
#define _PROC_GB_B1_ gather_boundary_XD1
#define _PROC_GB_B2_ gather_boundary_XD2
#define _PROC_GB_B3_ gather_boundary_XD3
#define _AUX_   gather_boundary_D_aux
#define _PROC_BIC1_ boundary_is_current_D1
#define _PROC_BIC2_ boundary_is_current_D2
#define _PROC_BIC3_ boundary_is_current_D3
#define _EQUAL_ ==
#endif

#ifndef _TYPE_
#error "one of LOGICAL_DATA, INTEGER_DATA, SINGLE_DATA, DOUBLE_DATA must be defined"
#endif

  subroutine _PROC_GB_A1_ (this, onP_data, offP_data)

    type(ip_desc), intent(in) :: this ! partition descriptor
    _TYPE_, intent(in)  :: onP_data(:)   ! on-process data array
    _TYPE_, intent(out) :: offP_data(:)  ! off-process data array

    ASSERT( this%defined() )
    ASSERT( associated(this%offP_index) )
    ASSERT( size(onP_data) == this%onP_size_ )
    ASSERT( size(offP_data) == this%offP_size_ )

    call _AUX_ (this, onP_data, offP_data)

  end subroutine _PROC_GB_A1_

  subroutine _PROC_GB_A2_ (this, onP_data, offP_data)

    type(ip_desc), intent(in)  :: this  ! partition descriptor
    _TYPE_, intent(in)  :: onP_data(:,:)  ! on-process data array
    _TYPE_, intent(out) :: offP_data(:,:) ! off-process data array

    integer :: k

    ASSERT( this%defined() )
    ASSERT( associated(this%offP_index) )
    ASSERT( size(offP_data,1) == size(onP_data,1) )
    ASSERT( size(onP_data,2) == this%onP_size_ )
    ASSERT( size(offP_data,2) == this%offP_size_ )

    do k = 1, size(onP_data,dim=1)   ! Ugh...
      call _AUX_ (this, onP_data(k,:), offP_data(k,:))
    end do

  end subroutine _PROC_GB_A2_

  subroutine _PROC_GB_A3_ (this, onP_data, offP_data)

    type(ip_desc), intent(in)  :: this  ! partition descriptor
    _TYPE_, intent(in)  :: onP_data(:,:,:)  ! on-process data array
    _TYPE_, intent(out) :: offP_data(:,:,:) ! off-process data array

    integer :: k, j

    ASSERT( this%defined() )
    ASSERT( associated(this%offP_index) )
    ASSERT( size(offP_data,1) == size(onP_data,1) )
    ASSERT( size(offP_data,2) == size(onP_data,2) )
    ASSERT( size(onP_data,3) == this%onP_size_ )
    ASSERT( size(offP_data,3) == this%offP_size_ )

    do j = 1, size(onP_data,dim=2)   ! Double ugh...
      do k = 1, size(onP_data,dim=1)
        call _AUX_ (this, onP_data(k,j,:), offP_data(k,j,:))
      end do
    end do

  end subroutine _PROC_GB_A3_

  subroutine _PROC_GB_B1_ (this, local_data)

    type(ip_desc), intent(in) :: this   ! partition descriptor
    _TYPE_, intent(inout) :: local_data(:)  ! local data array

    ASSERT( this%defined() )
    ASSERT( associated(this%offP_index) )
    ASSERT( size(local_data) == this%local_size_ )

    call _AUX_ (this, local_data(:this%onP_size_), local_data(this%onP_size_+1:))

  end subroutine _PROC_GB_B1_

  subroutine _PROC_GB_B2_ (this, local_data)

    type(ip_desc), intent(in) :: this   ! partition descriptor
    _TYPE_, intent(inout) :: local_data(:,:)  ! local data array

    integer :: k

    ASSERT( this%defined() )
    ASSERT( associated(this%offP_index) )
    ASSERT( size(local_data,2) == this%local_size_ )

    do k = 1, size(local_data,dim=1)   ! Ugh...
      call _AUX_ (this, local_data(k,:this%onP_size_), local_data(k,this%onP_size_+1:))
    end do

  end subroutine _PROC_GB_B2_

  subroutine _PROC_GB_B3_ (this, local_data)

    type(ip_desc), intent(in) :: this   ! partition descriptor
    _TYPE_, intent(inout) :: local_data(:,:,:)  ! local data array

    integer :: k, j

    ASSERT( this%defined() )
    ASSERT( associated(this%offP_index) )
    ASSERT( size(local_data,3) == this%local_size_ )

    do j = 1, size(local_data,dim=2)   ! Double ugh...
      do k = 1, size(local_data,dim=1)
        call _AUX_ (this, local_data(k,j,:this%onP_size_), local_data(k,j,this%onP_size_+1:))
      end do
    end do

  end subroutine _PROC_GB_B3_

  subroutine _AUX_ (this, onP_data, offP_data)

    type(ip_desc), intent(in) :: this ! partition descriptor
    _TYPE_, intent(in)  :: onP_data(:)   ! on-process data array
    _TYPE_, intent(out) :: offP_data(:)  ! off-process data array

    integer :: j
    _TYPE_ :: dup_buffer(size(this%dup_index))
    _TYPE_, allocatable :: sup_buffer(:)

    !! Gather the on-process boundary data into the duplicate buffer
    do j = 1, size(this%dup_index)
      dup_buffer(j) = onP_data(this%dup_index(j))
    end do

    if (associated(this%sup_index)) then ! SUP_BUFFER /= OFFP_DATA

      !! Global communication: DUP_BUFFER -> SUP_BUFFER.
      allocate(sup_buffer(pgslib_size_of_sup(this%trace)))
      sup_buffer = pgslib_gather_buffer(dup_buffer, this%trace)

      !! Gather the off-process data from the supplement buffer.
      do j = 1, size(this%sup_index)
        offP_data(j) = sup_buffer(this%sup_index(j))
      end do
      deallocate(sup_buffer)

    else  ! SUP_BUFFER == OFFP_DATA

      !! Global communication: DUP_BUFFER -> SUP_BUFFER (= OFFP_DATA).
      offP_data = pgslib_gather_buffer(dup_buffer, this%trace)

    end if

  end subroutine _AUX_

  logical function _PROC_BIC1_ (this, local_data) result (l)
    type(ip_desc), intent(in) :: this
    _TYPE_, intent(in) :: local_data(:)
    _TYPE_:: offP_data(this%offP_size_)
    call gather_boundary (this, local_data(:this%onP_size_), offP_data)
    l = global_all(local_data(this%onP_size_+1:) _EQUAL_ offP_data)
  end function _PROC_BIC1_
  
  logical function _PROC_BIC2_ (this, local_data) result (l)
    type(ip_desc), intent(in) :: this
    _TYPE_, intent(in) :: local_data(:,:)
    _TYPE_ :: offP_data(size(local_data,1),this%offP_size_)
    call gather_boundary (this, local_data(:,:this%onP_size_), offP_data)
    l = global_all(local_data(:,this%onP_size_+1:) _EQUAL_ offP_data)
  end function _PROC_BIC2_

! There is no specific for GLOBAL_ALL to handle rank-3 matrices.  Argh...
!  logical function _PROC_BIC3_ (this, local_data) result (l)
!    type(ip_desc), intent(in) :: this
!    _TYPE_, intent(in) :: local_data(:,:,:)
!    _TYPE_ :: offP_data(size(local_data,1),size(local_data,2),this%offP_size_)
!    call gather_boundary (this, local_data(:,:,:this%onP_size_), offP_data)
!    l = global_all(local_data(:,:,this%onP_size_+1:) _EQUAL_ offP_data)
!  end function _PROC_BIC3_
  
#undef _TYPE_
#undef _PROC_GB_A1_
#undef _PROC_GB_A2_
#undef _PROC_GB_B1_
#undef _PROC_GB_B2_
#undef _PROC_GB_A3_
#undef _PROC_GB_B3_
#undef _AUX_
#undef _PROC_BIC1_
#undef _PROC_BIC2_
#undef _PROC_BIC3_
#undef _EQUAL_


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!  DANU
!     Fortran Unit Test Utility module
!

module funit_utils

    use, intrinsic :: iso_c_binding

    implicit none

    public

! ==============================================================================
! Interfaces
! ==============================================================================

    ! C utility functions
    interface
        subroutine exit_now(code) bind(c)
        use, intrinsic :: iso_c_binding
        integer(C_INT), value  ::code
        end subroutine exit_now
    end interface

    interface
        subroutine fail_exit_now() bind(c)
        end subroutine fail_exit_now
    end interface

    interface
        subroutine pass_exit_now() bind(c)
        end subroutine pass_exit_now
    end interface

    interface
        subroutine create_test_h5_file() bind(c)
        end subroutine create_test_h5_file
    end interface    

    interface
        subroutine delete_test_file() bind(c)
        end subroutine delete_test_file
    end interface    

    interface
        integer(C_INT64_T) function open_test_file() bind(c)
        use, intrinsic :: iso_c_binding
        end function open_test_file
    end interface

    interface
        subroutine close_test_file(id) bind(c)
        use, intrinsic :: iso_c_binding
        integer(C_INT64_T), value :: id
        end subroutine close_test_file
    end interface

    interface
        integer(C_INT) function generate_random_int(seed) &
                                bind(c)
        use, intrinsic :: iso_c_binding
        integer(C_INT), intent(inout) :: seed
        end function generate_random_int
    end interface

    interface
        integer(C_INT) function generate_random_bound_int(imin,imax,seed) &
                                bind(c)
        use, intrinsic :: iso_c_binding
        integer(C_INT), value         :: imin
        integer(C_INT), value         :: imax
        integer(C_INT), intent(inout) :: seed
        end function generate_random_bound_int
    end interface

    interface
        real(C_FLOAT) function generate_random_float(seed) &
                                bind(c)
        use, intrinsic :: iso_c_binding
        integer(C_INT), intent(inout) :: seed
        end function generate_random_float
    end interface

    interface
        real(C_FLOAT) function generate_random_bound_float(fmin,fmax,seed) &
                                bind(c)
        use, intrinsic :: iso_c_binding
        real(C_FLOAT), value         :: fmin
        real(C_FLOAT), value         :: fmax
        integer(C_INT), intent(inout) :: seed
        end function generate_random_bound_float
    end interface

    interface
        real(C_DOUBLE) function generate_random_double(seed) &
                                bind(c)
        use, intrinsic :: iso_c_binding
        integer(C_INT), intent(inout) :: seed
        end function generate_random_double
    end interface

    interface
        real(C_DOUBLE) function generate_random_bound_double(dmin,dmax,seed) &
                                bind(c)
        use, intrinsic :: iso_c_binding
        real(C_DOUBLE), value         :: dmin
        real(C_DOUBLE), value         :: dmax
        integer(C_INT), intent(inout) :: seed
        end function generate_random_bound_double
    end interface

    interface
         type(C_PTR) function create_hid_struct(id) bind(c)
         use, intrinsic :: iso_c_binding
         integer(C_INT64_T), value     :: id
         end function create_hid_struct
    end interface     

    interface
         subroutine destroy_hid_struct(ptr) bind(c)
         use, intrinsic :: iso_c_binding
         type(C_PTR), value     :: ptr
         end subroutine destroy_hid_struct
    end interface   

    interface
         integer(C_INT) function int_array_cmp(a,b,num) bind(c)
         use, intrinsic :: iso_c_binding
         integer(C_INT), intent(in) :: a(*)      
         integer(C_INT), intent(in) :: b(*)
         integer(C_INT), value, intent(in) :: num
         end function int_array_cmp
    end interface     

    interface
         integer(C_INT) function double_array_cmp(a,b,num) bind(c)
         use, intrinsic :: iso_c_binding
         real(C_DOUBLE), intent(in) :: a(*)      
         real(C_DOUBLE), intent(in) :: b(*)      
         integer(C_INT), value, intent(in) :: num
         end function double_array_cmp
    end interface     

    interface
         integer(C_INT) function float_array_cmp(a,b,num) bind(c)
         use, intrinsic :: iso_c_binding
         real(C_FLOAT), intent(in) :: a(*)      
         real(C_FLOAT), intent(in) :: b(*)      
         integer(C_INT), value, intent(in) :: num
         end function float_array_cmp
    end interface   

    interface generate_random_array
      module procedure gen_array_int_rank1
      module procedure gen_array_int_rank2
      module procedure gen_array_int_rank3
      module procedure gen_array_real4_rank1
      module procedure gen_array_real4_rank2
      module procedure gen_array_real4_rank3
      module procedure gen_array_real8_rank1
      module procedure gen_array_real8_rank2
      module procedure gen_array_real8_rank3
    end interface generate_random_array  

    interface generate_random_bound_array
      module procedure gen_array_bound_int_rank1
      module procedure gen_array_bound_int_rank2
      module procedure gen_array_bound_int_rank3
    end interface generate_random_bound_array  

! ==============================================================================    
    contains
! ==============================================================================    

! ------------------------------------------------------------------------------

subroutine gen_array_int_rank1(array,iseed)

! --- Calling arguments

      integer,dimension(:),intent(inout) :: array
      integer,             intent(inout) :: iseed

! --- Local variables     
     
      integer,dimension(1) :: dims
      integer :: i

      dims = shape(array)
      i = 1
      do while ( i .le. dims(1) )
        array(i) = generate_random_int(iseed)
        i = i + 1
      end do  

end subroutine gen_array_int_rank1

! ------------------------------------------------------------------------------

subroutine gen_array_int_rank2(array,iseed)

! --- Calling arguments

      integer,dimension(:,:),intent(inout) :: array
      integer,               intent(inout) :: iseed

! --- Local variables     
     
      integer,dimension(2) :: dims
      integer :: i, j

      dims = shape(array)
      i = 1
      do while ( i .le. dims(1) )
        j = 1
        do while ( j .le. dims(2) )
          array(i,j) = generate_random_int(iseed)
          j = j + 1
        end do  
        i = i + 1
      end do  

end subroutine gen_array_int_rank2

! ------------------------------------------------------------------------------

subroutine gen_array_int_rank3(array,iseed)

! --- Calling arguments

      integer,dimension(:,:,:),intent(inout) :: array
      integer,                 intent(inout) :: iseed

! --- Local variables     
     
      integer,dimension(3) :: dims
      integer :: i, j, k

      dims = shape(array)
      i = 1
      do while ( i .le. dims(1) )
        j = 1
        do while ( j .le. dims(2) )
          k = 1
          do while ( k .le. dims(3) )
            array(i,j,k) = generate_random_int(iseed)
            k = k + 1
          end do  
          j = j + 1
        end do  
        i = i + 1
      end do  

end subroutine gen_array_int_rank3

! ------------------------------------------------------------------------------

subroutine gen_array_real4_rank1(array,iseed)

! --- Calling arguments

      real(C_FLOAT),dimension(:),intent(inout) :: array
      integer,             intent(inout) :: iseed

! --- Local variables     
     
      integer,dimension(1) :: dims
      integer :: i

      dims = shape(array)
      i = 1
      do while ( i .le. dims(1) )
        array(i) = generate_random_float(iseed)
        i = i + 1
      end do  

end subroutine gen_array_real4_rank1

! ------------------------------------------------------------------------------

subroutine gen_array_real4_rank2(array,iseed)

! --- Calling arguments

      real(C_FLOAT),dimension(:,:),intent(inout) :: array
      integer,               intent(inout) :: iseed

! --- Local variables     
     
      integer,dimension(2) :: dims
      integer :: i, j

      dims = shape(array)
      i = 1
      do while ( i .le. dims(1) )
        j = 1
        do while ( j .le. dims(2) )
          array(i,j) = generate_random_float(iseed)
          j = j + 1
        end do  
        i = i + 1
      end do  

end subroutine gen_array_real4_rank2

! ------------------------------------------------------------------------------

subroutine gen_array_real4_rank3(array,iseed)

! --- Calling arguments

      real(C_FLOAT),dimension(:,:,:),intent(inout) :: array
      integer,                 intent(inout) :: iseed

! --- Local variables     
     
      integer,dimension(3) :: dims
      integer :: i, j, k

      dims = shape(array)
      i = 1
      do while ( i .le. dims(1) )
        j = 1
        do while ( j .le. dims(2) )
          k = 1
          do while ( k .le. dims(3) )
            array(i,j,k) = generate_random_float(iseed)
            k = k + 1
          end do  
          j = j + 1
        end do  
        i = i + 1
      end do  

end subroutine gen_array_real4_rank3

! ------------------------------------------------------------------------------

subroutine gen_array_real8_rank1(array,iseed)

! --- Calling arguments

      real(C_DOUBLE),dimension(:),intent(inout) :: array
      integer,             intent(inout) :: iseed

! --- Local variables     
     
      integer,dimension(1) :: dims
      integer :: i

      dims = shape(array)
      i = 1
      do while ( i .le. dims(1) )
        array(i) = generate_random_double(iseed)
        i = i + 1
      end do  

end subroutine gen_array_real8_rank1

! ------------------------------------------------------------------------------

subroutine gen_array_real8_rank2(array,iseed)

! --- Calling arguments

      real(C_DOUBLE),dimension(:,:),intent(inout) :: array
      integer,               intent(inout) :: iseed

! --- Local variables     
     
      integer,dimension(2) :: dims
      integer :: i, j

      dims = shape(array)
      i = 1
      do while ( i .le. dims(1) )
        j = 1
        do while ( j .le. dims(2) )
          array(i,j) = generate_random_double(iseed)
          j = j + 1
        end do  
        i = i + 1
      end do  

end subroutine gen_array_real8_rank2

! ------------------------------------------------------------------------------

subroutine gen_array_real8_rank3(array,iseed)

! --- Calling arguments

      real(C_DOUBLE),dimension(:,:,:),intent(inout) :: array
      integer,                 intent(inout) :: iseed

! --- Local variables     
     
      integer,dimension(3) :: dims
      integer :: i, j, k

      dims = shape(array)
      i = 1
      do while ( i .le. dims(1) )
        j = 1
        do while ( j .le. dims(2) )
          k = 1
          do while ( k .le. dims(3) )
            array(i,j,k) = generate_random_double(iseed)
            k = k + 1
          end do  
          j = j + 1
        end do  
        i = i + 1
      end do  

end subroutine gen_array_real8_rank3

! ------------------------------------------------------------------------------

subroutine gen_array_bound_int_rank1(imin,imax,array,iseed)

! --- Calling arguments

      integer,             intent(in)    :: imin
      integer,             intent(in)    :: imax
      integer,dimension(:),intent(inout) :: array
      integer,             intent(inout) :: iseed

! --- Local variables     
     
      integer,dimension(1) :: dims
      integer :: i

      dims = shape(array)
      i = 1
      do while ( i .le. dims(1) )
        array(i) = generate_random_bound_int(imin,imax,iseed)
        i = i + 1
      end do  

end subroutine gen_array_bound_int_rank1

! ------------------------------------------------------------------------------

subroutine gen_array_bound_int_rank2(imin,imax,array,iseed)

! --- Calling arguments

      integer,               intent(in)    :: imin
      integer,               intent(in)    :: imax
      integer,dimension(:,:),intent(inout) :: array
      integer,               intent(inout) :: iseed

! --- Local variables     
     
      integer,dimension(2) :: dims
      integer :: i, j

      dims = shape(array)
      i = 1
      do while ( i .le. dims(1) )
        j = 1
        do while ( j .le. dims(2) )
          array(i,j) = generate_random_bound_int(imin,imax,iseed)
          j = j + 1
        end do  
        i = i + 1
      end do  

end subroutine gen_array_bound_int_rank2

! ------------------------------------------------------------------------------

subroutine gen_array_bound_int_rank3(imin,imax,array,iseed)

! --- Calling arguments

      integer,                 intent(in)    :: imin
      integer,                 intent(in)    :: imax
      integer,dimension(:,:,:),intent(inout) :: array
      integer,                 intent(inout) :: iseed

! --- Local variables     
     
      integer,dimension(3) :: dims
      integer :: i, j, k

      dims = shape(array)
      i = 1
      do while ( i .le. dims(1) )
        j = 1
        do while ( j .le. dims(2) )
          k = 1
          do while ( k .le. dims(3) )
            array(i,j,k) = generate_random_bound_int(imin,imax,iseed)
            k = k + 1
          end do  
          j = j + 1
        end do  
        i = i + 1
      end do  

end subroutine gen_array_bound_int_rank3

! =============================================================================
end module funit_utils    
! =============================================================================

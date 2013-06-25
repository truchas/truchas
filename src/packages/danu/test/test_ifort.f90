program test_ifort

    implicit none

    integer :: flag

    flag = 2

    call dummy_fort_interface(flag)

end program    

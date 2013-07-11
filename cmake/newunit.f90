subroutine conftest
    open(newunit=n,file='conftest.out')
end subroutine
program main
    call conftest
end program main

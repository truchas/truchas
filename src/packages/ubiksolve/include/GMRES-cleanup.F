
    ! deallocate local arrays
    call destroy (x_vec, overlap_only=.true.)
    call destroy (xt_vec, overlap_only=.true.)
    if (allocated(v)) then
       deallocate(v)
    end if
    if (allocated(h)) then
       deallocate(h)
    end if
    if (allocated(rt)) then
       deallocate(rt)
    end if
    if (allocated(r)) then
       deallocate(r)
    end if
    if (allocated(s)) then
       deallocate(s)
    end if
    if (allocated(c)) then
       deallocate(c)
    end if
    if (allocated(xt)) then
       deallocate(xt)
    end if
    if (allocated(w)) then
       deallocate(w)
    end if

    ! exit w/o cleanup
999 continue

    ! Subtract time at start of solution from each element of cpu array
    ! to obtain true CPU time for each iteration.
    call normalize_cpu (control)

#include "formats.F"

The gcc-* configuration files, which use the GNU Fortran compiler gfortran,
are provided for testing purposes only.  The gfortran compiler is not supported
by Truchas, and does not currently work (versions 5.4 and 6.1 have been tested).

At this time gfortran has incomplete and/or flawed support of derived type
finalization (PR37336), unlimited polymorphic variables (PR67564), and
deferred-length character variables (PR68241).  The latter two particularly
impact the Petaca library, while the first impacts all of Truchas.  These
shortcomings may result in code not executing properly and memory leaks.
Visit <https://gcc.gnu.org/bugzilla/query.cgi> and search on the above bug
numbers (omit the "PR") for their latest status.

There is also a bug in both 5.4 and 6.1 that trigger an internal compiler
error in exodus_truchas_hack.F90.  This has been reported; see
<https://gcc.gnu.org/bugzilla/show_bug.cgi?id=77310>

NNC, August 2016

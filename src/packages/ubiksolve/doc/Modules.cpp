/* -*-Mode: Fundamental;-*- */

Notes about JTpack90 modules
____________________________

Last modified: 11/04/95, 23:17:14

[Page numbers refer to Fortran 90 Programming, Ellis, et al.]

- Modules must have unique names.  That is, a routine named MatrixNorm
  cannot be contained in a module named MatrixNorm.  I have chosen the
  convention that names of modules containing routines are typically
  the routine names with _module appended.

- A single "implicit none" in each module is sufficient.

- Some modules contain routines, others contain only groups of other
  modules.

- In modules containing routines, the following rules have been
  followed:

 - Use association is explicit for any entities explicitly needed in the module.
   For example, modules containing routines always need JT_types_module.  A use 
   statement for it appears even if a use statement for some other module that
   uses it also appears.
 - Use statements for specific routines appear, rather than for module groups,
   unless all members of a module group are used.  That is, if a module needs
   several BLAS-1 routines, use statements for the specific routines appear
   rather than a single use statement for JT_BLAS1_module.

4) Need to think about PUBLIC and PRIVATE attributes (p. 405).

5) Encapsulation of routines in modules and making them available via
   use association makes interfaces explicit, and hence allows use of
   assumed-shape arrays (p. 214).  This eliminates the need to pass in
   leading dimensions of 2-D arrays, etc.  But it leads to questions
   about how to implement routines that might operate on only a part
   of an array.  Options:

  - dimension as assumed-shape arrays, use whole array operations in the
    subroutine
   > requires passing in appropriate array sections

  real, dimension(5,5) :: x, y, z
    :
  call madd (x(1:3,:)), y(1:3,:)), z(1:3,:)))

  subroutine madd (x, y, z)
    real, dimension(:,:) :: x, y, z
    z = x + y
  return
  end

   > makes caller messy, callee clean

  - use explicit-shape dimensioning
   > requires passing in dimensions as well as info about elements to
     be operated upon

  real, dimension(5,5) :: x, y, z
    :
  call madd (5, 3, 5, x, y, z)

  subroutine madd (ld, nrows, ncols, x, y, z)
    real, dimension(ld,ncols) :: x, y, z
    z(1:nrows,:) = x(1:nrows,:) + y(1:nrows,:)
  return
  end

   > this example makes the caller clean, and the callee messy
   > it's easy to think of examples that would require array sections
     to be passed in, making the caller just as messy as before, e.g.:

  call madd (5, 3, 2, x(1:3,1:2)), y(2:4,2:3)), z(3:5,3:4)))

   > the limited flexibility of this approach is obvious (really should
     pass in a leading dimension for each matrix, etc.)

  - use assumed-shape dimensioning, but pass in, e.g., number of elements
    to be operated on

  real, dimension(5,5) :: x, y, z
    :
  call madd (3, 5, x, y, z)

  subroutine madd (nrows, ncols, x, y, z)
    real, dimension(:,:) :: x, y, z
    z(1:nrows,1:ncols) = x(1:nrows,1:ncols) + y(1:nrows,1:ncols)
  return
  end

   > this example makes the caller clean, and the callee messy
   > as with the previous option, it's easy to think of examples that
     would require array sections to be passed in, making the caller just
     as messy as before

  - looks like option 1 is best
   > option 2 is far too inflexible
   > option 3 adds to argument list, and doesn't really buy anything
   > consider doing the following with either option 2 or option 3

  real, dimension(3,5) :: x
  real, dimension(4,5) :: y
  real, dimension(5,5) :: z
    :
  call madd (x(:,1:2)), y(2:4,2:3)), z(3:5,3:4)))

  subroutine madd (x, y, z)
    real, dimension(:,:) :: x, y, z
    z = x + y
  return
  end

6) Need to think about opportunities for array-valued functions (p.220).

7) Need to add tests for conformability in a few routines.

Details about modules in JTpack90:

#include "Modules.list"

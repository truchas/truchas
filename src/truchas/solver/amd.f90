! takes a quick look at the T matrix dump (in Harwell-Boeing format)
! compile with almost any options: f90 -O -o amd amd.f90
! invoke as "amd infile"

! we know what HB format options we dump with and take advantage of that!

program AnalyzeMatrixDump

   implicit none

   character (len=1024)              :: filename
   real, allocatable, dimension(:,:) :: A
   real, allocatable, dimension(:)   :: b

!   integer                           :: neqn
!   integer                           :: column
!   integer                           :: count
!   real                              :: value
   integer                           :: IARGC
   integer                           :: i
   integer                           :: j
   integer                           :: status

   ! harwell boeing variables
   character (len=72) :: title
   character (len=8)  :: key
   character (len=3)  :: mxtype
   character (len=3)  :: rhstyp
   character (len=16) :: ptrfmt
   character (len=16) :: indfmt
   character (len=20) :: valfmt
   character (len=20) :: rhsfmt
   integer :: totcrd
   integer :: ptrcrd
   integer :: indcrd
   integer :: valcrd
   integer :: rhscrd
   integer :: nrow
   integer :: ncol
   integer :: nnzero
   integer :: neltvl
   integer :: nrhs
   integer :: nrhsix
   integer, dimension(:), allocatable :: colptr
   integer, dimension(:), allocatable :: rowind
   real,    dimension(:), allocatable :: values
   real,    dimension(:), allocatable :: rhsval

   !----------------------------------------------------------------------------

   ! see if we supplied an input filename argument
   if (IARGC() /= 1) then
      write (*,*) 'usage: amd filename'
      stop
   end if

   ! get the input filename argument
   call GETARG (1, filename)

   ! open the input file
   open (10, FILE=trim(filename), status='OLD', iostat=status)
   if (status /= 0) then
      write (*,*) 'can''t open '//TRIM(filename)
      stop
   end if

   ! read header block
   read (10, 1001) title, key
   read (10, 1002) totcrd, ptrcrd, indcrd, valcrd, rhscrd
   read (10, 1003) mxtype, nrow, ncol, nnzero, neltvl
   read (10, 1004) ptrfmt, indfmt, valfmt, rhsfmt
   read (10, 1005) rhstyp, nrhs, nrhsix

   ! allocate space to store the system
   allocate (colptr(ncol+1), rowind(nnzero), values(nnzero), rhsval(nrow), STAT=status)
   if (status /= 0) then
      write (*,*) 'insufficient memory to allocate working arrays for Harwell-Boeing format'
      stop
   end if
   allocate (A(nrow,ncol), b(nrow), STAT=status)
   if (status /= 0) then
      write (*,*) 'insufficient memory to allocate working arrays for ', nrow, ' equations'
      stop
   end if

   ! read matrix
   read (10, ptrfmt) (colptr(i), i = 1, ncol+1)
   read (10, indfmt) (rowind(i), i = 1, nnzero)
   read (10, valfmt) (values(i), i = 1, nnzero)
   read (10, rhsfmt) (b(i), i = 1, nrow)

   ! put the HB matrix in the full matrix
   A = 0.0
   do i = 1, ncol
      do j = colptr(i), colptr(i+1)-1
         A(rowind(j),i) = values(j)
      end do
   end do

   ! examine matrix we've just read in
   status = Analyze(A)

   ! deallocate the matrix and vectors
   deallocate (A,b)

   ! close the input file
   close (10)

1001 format (a72, a8)
1002 format (5i14)
1003 format (a3, 11x, 4i14)
1004 format (2a16, 2a20)
1005 format (a3, 11x, 2i14)

contains

   !----------------------------------------------------------------------------
   ! examination routines
   !----------------------------------------------------------------------------

   function Analyze (A)

      implicit none

      real, dimension(:,:) :: A
      integer :: Analyze

      integer :: r
      integer :: c
      integer :: n
      integer :: nmax
      real, parameter :: TOLERANCE = 1.0e-8 ! this is a single precision computation!

      Analyze = 0
      write (*,*) 'size: ', SIZE(A,1)

      ! check for maximum number of coefficients in columns
      nmax = 0
      do c = 1, SIZE(A,1)
         n = 0
         do r = 1, SIZE(A,1)
            if (A(r,c) /= 0.0) n = n+1
         end do
         if (n > nmax) nmax = n
      end do
      write (*,*) 'maximum coefficients in a column: ', nmax

      ! check for maximum number of coefficients in rows
      nmax = 0
      do r = 1, SIZE(A,1)
         n = 0
         do c = 1, SIZE(A,1)
            if (A(r,c) /= 0.0) n = n+1
         end do
         if (n > nmax) nmax = n
      end do
      write (*,*) 'maximum coefficients in a row: ', nmax

      ! check for structural symmetry
      ! this check could fail if we calculate a clean zero somewhere
      do c = 1, SIZE(A,1)
         do r = c, SIZE(A,1)
            if ((A(r,c) == 0.0 .and. A(r,c) /= 0.0) .or. (A(r,c) == 0.0 .and. A(r,c) /= 0.0)) then
               write (*,*) 'asymmetric structure at (r,c): ', r, c
               Analyze = 1
               return
            end if
         end do
      end do
      write (*,*) 'structurally symmetric'

      ! check for numerical symmetry
      do c = 1, SIZE(A,1)
         do r = c, SIZE(A,1)
            if (A(r,c) /= 0.0) then
               if (ABS(A(r,c)-A(c,r))/ABS(A(r,c)) >= TOLERANCE) then
                  Analyze = 2
                  write (*,*) 'asymmetric values (r,c): ', r, c
                  write (*,*) 'values: ', A(r,c), A(c,r), ABS(A(r,c)-A(c,r))/ABS(A(r,c))
                  return
               end if
            end if
         end do
      end do
      write (*,*) 'numerically symmetric'

      return

   end function Analyze

   !----------------------------------------------------------------------------

end program AnalyzeMatrixDump

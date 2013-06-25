! meshtranslate.f90

module translate

   implicit none
   private
   public :: readpatran,  &
             readideas,   &
             readavs,     &
             analyzemesh, &
             writegmv,    &
             writeensight

   integer :: cells
   integer :: nodes
   integer :: materials
   real*8,  pointer, dimension(:,:) :: coordinates
   integer, pointer, dimension(:,:) :: connectivity
   integer, pointer, dimension(:)   :: material

contains

   !----------------------------------------------------------------------------

   function readpatran (fname)

      implicit none

      ! arguments
      character (LEN=*) :: fname

      ! return
      integer :: readpatran

      ! local variables
      integer, parameter :: lun = 10
      integer            :: status
      integer            :: n           ! temporary
      integer            :: i           ! loop index
      integer            :: j           ! loop index
      integer            :: node        ! node number
      real*8             :: r           ! temporary
      integer :: cell
      integer :: nodesincell

      !-------------------------------------------------------------------------

      write (*,*) 'reading patran ...'

      ! open input file
      open (lun, FILE=fname, STATUS='old', IOSTAT=status)
      if (status /= 0) then
         write (*,*) 'can''t open input file '//fname
         stop
      end if
      write (*,*) '  from file '//fname

      ! ignore two lines
      read (lun, *)
      read (lun, *)

      ! get the mesh sizes
      read (lun, *) n, n, n, n, nodes, cells, n, n, n

      ! ignore one line
      read (lun, *)

      ! allocate data structures
      allocate (coordinates(nodes,3))
      allocate (connectivity(cells,8))
      allocate (material(cells))

      connectivity = 0

      ! read node data
      do i = 1, nodes
         read (lun, *) n, node
         read (lun, 100) (coordinates(node,j), j=1,3)
         read (lun, *)
      end do

      ! read connectivity and material data
      do i = 1, cells
         read (lun, *) n, cell, n, n, n
         read (lun, *) nodesincell, n, material(cell), n, r, r, r
         read (lun, *) (connectivity(cell,j), j=1,nodesincell)
      end do

      close (lun)

      readpatran = 0
      return

100   format (3(1pe16.9))

   end function readpatran

   !----------------------------------------------------------------------------

   function readideas (fname)

      implicit none

      ! arguments
      character (LEN=*) :: fname

      ! return
      integer :: readideas

      ! local variables
      integer, parameter :: lun = 10
      integer            :: status
      integer            :: n           ! temporary
      integer            :: i           ! loop index
      integer            :: j           ! loop index
      integer            :: k           ! loop index
      integer :: cell
      integer :: nodesincell
      logical :: done
      character (LEN=256) :: string

      !-------------------------------------------------------------------------

      write (*,*) 'reading ideas ...'

      ! open input file
      open (lun, FILE=fname, STATUS='old', IOSTAT=status)
      if (status /= 0) then
         write (*,*) 'can''t open input file '//fname
         stop
      end if
      write (*,*) '  from file '//fname

      ! count nodes and cells
      cells = 0
      nodes = 0
      do i = 1, 2
         done = .false.

         ! ignore one line
         read (lun, *)

         ! read a string
         read (lun, *) string

         ! count nodes and cells
         select case (TRIM(ADJUSTL(string)))
         case ('2411')
            ! nodes
            do 
               if (done) exit
               read (lun, *) string
               if (TRIM(ADJUSTL(string)) == '-1') then
                  done = .true.
               else
                  nodes = nodes + 1
               end if
            end do
            nodes = nodes / 2
         case ('2412')
            ! cells
            do 
               if (done) exit
               read (lun, *) string
               if (TRIM(ADJUSTL(string)) == '-1') then
                  done = .true.
               else
                  cells = cells + 1
               end if
            end do
            cells = cells / 2
         end select
      end do

      rewind (lun)

      ! allocate data structures
      allocate (coordinates(nodes,3))
      allocate (connectivity(cells,8))
      allocate (material(cells))
      connectivity = 0
      material = 0

      ! read nodes and cells
      do i = 1, 2
         done = .false.

         ! ignore one line
         read (lun, *)

         ! read a string
         read (lun, *) string

         ! count nodes and cells
         select case (TRIM(ADJUSTL(string)))
         case ('2411')
            ! nodes
            do j = 1, nodes
               read (lun, *) n
               read (lun, *) (coordinates(n,k), k = 1, 3)
            end do
            read (lun,*)
         case ('2412')
            ! cells
            do j = 1, cells
               read (lun, *) cell, n, n, n, n, nodesincell
               read (lun, *) (connectivity(cell,k), k = 1, nodesincell)
            end do
            read (lun,*)
         end select
      end do

      ! look for optional materials
      do
         read (lun,*,end=1000) string
         if (TRIM(ADJUSTL(string)) == 'matl') then
            do i = 1, cells
               read (lun,*) material(i)
            end do
         end if
      end do

1000  close (lun)

      readideas = 0
      return

   end function readideas

   !----------------------------------------------------------------------------

   function readavs (fname)

      implicit none

      ! arguments
      character (LEN=*) :: fname

      ! return
      integer :: readavs

      ! local variables
      integer, parameter :: lun = 10
      integer            :: status
      integer            :: n           ! temporary
      integer            :: i           ! loop index
      integer            :: j           ! loop index
      character (LEN=8)  :: c           ! character temporary

      !-------------------------------------------------------------------------

      write (*,*) 'reading avs ...'

      ! open input file
      open (lun, FILE=fname, STATUS='old', IOSTAT=status)
      if (status /= 0) then
         write (*,*) 'can''t open input file '//fname
         stop
      end if
      write (*,*) '  from file '//fname

      ! read number of nodes and cells
      read (lun, '(i10,4(x,i10))') nodes, cells

      ! allocate data structures
      allocate (coordinates(nodes,3))
      allocate (connectivity(cells,8))
      allocate (material(cells))
      connectivity = 0
      material = 0

      ! read the node coordinates
      do i = 1, nodes
         read (lun, '(i10,x,3(x,es19.12))') n, (coordinates(i,j), j = 1,3)
      end do

      ! read the connectivity
      do i = 1, cells
         read (lun, '(i10,x,i5,x,a8,8(x,i10))') n, n, c, (connectivity(i,j), j = 1, 8)
      end do

      ! don't know how to do avs materials

      close (lun)

      readavs = 0
      return

   end function readavs

   !----------------------------------------------------------------------------

   function writegmv (fname)

      implicit none

      ! arguments
      character (LEN=*) :: fname

      ! return value
      integer :: writegmv

      ! local variables
      integer, parameter :: lun = 10
      integer            :: status
      integer :: i
      integer :: j

      !-------------------------------------------------------------------------

      write (*,*) 'writing gmv ...'

      ! open output file
      open (lun, FILE=fname//'.gmv', STATUS='unknown', IOSTAT=status)
      if (status /= 0) then
         write (*,*) 'can''t open output file '//fname//'.gmv'
         stop
      end if
      write (*,*) '  to file '//fname//'.gmv'

      ! write prologue
      write (lun,100) 'gmvinput ascii'

      ! write node data
      write (lun,101) 'nodes ', nodes
      do i = 1, 3
         write (lun,102) (coordinates(j,i), j = 1, nodes)
      end do

      ! write connectivity
      write (lun,101) 'cells ', cells
      do i = 1, cells
         select case (COUNT(connectivity(i,:) /= 0))
         case (4)
            write (lun,104) 'hex 8 ', &
               connectivity(i,1), &
               connectivity(i,1), &
               connectivity(i,2), &
               connectivity(i,3), &
               connectivity(i,4), &
               connectivity(i,4), &
               connectivity(i,4), &
               connectivity(i,4)
         case (5)
            write (lun,104) 'hex 8 ', &
               connectivity(i,1), &
               connectivity(i,2), &
               connectivity(i,3), &
               connectivity(i,4), &
               connectivity(i,5), &
               connectivity(i,5), &
               connectivity(i,5), &
               connectivity(i,5)
         case (6)
            write (lun,104) 'hex 8 ', &
               connectivity(i,1), &
               connectivity(i,2), &
               connectivity(i,3), &
               connectivity(i,4), &
               connectivity(i,5), &
               connectivity(i,5), &
               connectivity(i,6), &
               connectivity(i,6)
         case (8)
            write (lun,104) 'hex 8 ', &
               connectivity(i,1), &
               connectivity(i,2), &
               connectivity(i,3), &
               connectivity(i,4), &
               connectivity(i,5), &
               connectivity(i,6), &
               connectivity(i,7), &
               connectivity(i,8)
         case default
            write (*,*) 'invalid number of vertices: ', COUNT(connectivity(i,:) /= 0)
            stop
         end select
      end do

      ! write materials
      write (lun,105) 'material ', materials, ' 0'

      do i = 1, materials
         write (lun,106) 'p', i
      end do

      write (lun,107) material

      ! write trailer
      write (lun,100) 'endgmv'

      close (lun)

      writegmv = 0
      return

100   format (a)
101   format (a,i8)
102   format (8(1es12.5,1x))
103   format (8(1es16.9,1x))
104   format (a,8(i8,1x))
105   format (a,i8,a)
106   format (a,i3.3)
107   format (8(i3.3,1x))

   end function writegmv

   !----------------------------------------------------------------------------

   function writeensight (fname)

      implicit none

      ! arguments
      character (LEN=*) :: fname

      ! return value
      integer :: writeensight

      ! local variables
      integer, parameter :: lun1 = 10
      integer, parameter :: lun2 = 11
      integer            :: status
      integer            :: i
      integer            :: j

      !-------------------------------------------------------------------------

      write (*,*) 'writing ensight ...'

      ! open output file
      open (lun1, FILE=fname//'.case', STATUS='unknown', IOSTAT=status)
      if (status /= 0) then
         write (*,*) 'can''t open output file '//fname//'.case'
         stop
      end if

      ! open output file
      open (lun2, FILE=fname//'.ensight', STATUS='unknown', IOSTAT=status)
      if (status /= 0) then
         write (*,*) 'can''t open output file '//fname//'.ensight'
         close (lun1)
         stop
      end if

      ! write case file
      write (*,*) '  to file '//fname//'.case'

      write (lun1,100) 'FORMAT'
      write (lun1,100) 'type:	ensight'
      write (lun1,*)
      write (lun1,100) 'GEOMETRY'
      write (lun1,100) 'model:	'//fname//'.ensight'

      close (lun1)

      ! write geometry file
      write (*,*) '  to file '//fname//'.ensight'
      write (lun2,100) 'EnSight geometry file by meshtranslate'
      write (lun2,*)
      write (lun2,100) 'node id assign'
      write (lun2,100) 'element id assign'

      ! write node data
      write (lun2,100) 'coordinates'
      write (lun2,101) nodes

      do i = 1, nodes
         write (lun2,102) (coordinates(i,j), j = 1, 3)
      end do

      do j = 1, materials

         write (lun2,104) 'part ', j
         write (lun2,101) j

         ! write connectivity
         write (lun2,100) 'hexa8'
         write (lun2,101) COUNT(material == j)

         do i = 1, cells
            if (material(i) == j) then
               select case (COUNT(connectivity(i,:) /= 0))
               case (4)
                  write (lun2,103)      &
                     connectivity(i,1), &
                     connectivity(i,1), &
                     connectivity(i,2), &
                     connectivity(i,3), &
                     connectivity(i,4), &
                     connectivity(i,4), &
                     connectivity(i,4), &
                     connectivity(i,4)
               case (5)
                  write (lun2,103)      &
                     connectivity(i,1), &
                     connectivity(i,2), &
                     connectivity(i,3), &
                     connectivity(i,4), &
                     connectivity(i,5), &
                     connectivity(i,5), &
                     connectivity(i,5), &
                     connectivity(i,5)
               case (6)
                  write (lun2,103)      &
                     connectivity(i,1), &
                     connectivity(i,2), &
                     connectivity(i,3), &
                     connectivity(i,4), &
                     connectivity(i,5), &
                     connectivity(i,5), &
                     connectivity(i,6), &
                     connectivity(i,6)
               case (8)
                  write (lun2,103)      &
                     connectivity(i,1), &
                     connectivity(i,2), &
                     connectivity(i,3), &
                     connectivity(i,4), &
                     connectivity(i,5), &
                     connectivity(i,6), &
                     connectivity(i,7), &
                     connectivity(i,8)
               case default
                  write (*,*) 'invalid number of vertices: ', COUNT(connectivity(i,:) /= 0)
                  stop
               end select
            end if
         end do
      end do

      close (lun2)

      writeensight = 0
      return

100   format (a)
101   format (i8)
102   format (3es12.5)
103   format (8i8)
104   format (a,i8)

   end function writeensight

   !----------------------------------------------------------------------------

   function analyzemesh ()
      implicit none

      ! return value
      integer :: analyzemesh

      ! local variables
      real*8 :: maxx
      real*8 :: minx
      real*8 :: maxy
      real*8 :: miny
      real*8 :: maxz
      real*8 :: minz
      integer :: i
      integer :: j
      integer, allocatable, dimension(:) :: mattmp

      !-------------------------------------------------------------------------

      write (*,*) 'analyzing mesh ...'

      write (*,*) '  # of cells:  ', cells
      write (*,*) '  # of nodes:  ', nodes

      ! find coordinate extent
      maxx = coordinates(1,1)
      minx = coordinates(1,1)
      maxy = coordinates(1,2)
      miny = coordinates(1,2)
      maxz = coordinates(1,3)
      minz = coordinates(1,3)

      do i = 2, nodes
         if (coordinates(i,1) > maxx) maxx = coordinates(i,1)
         if (coordinates(i,1) < minx) minx = coordinates(i,1)
         if (coordinates(i,2) > maxy) maxy = coordinates(i,2)
         if (coordinates(i,2) < miny) miny = coordinates(i,2)
         if (coordinates(i,3) > maxz) maxz = coordinates(i,3)
         if (coordinates(i,3) < minz) minz = coordinates(i,3)
      end do

      ! report coordinate extent
      write (*,'(a,2es12.4)') '   x range:     ', minx, maxx
      write (*,'(a,2es12.4)') '   y range:     ', miny, maxy
      write (*,'(a,2es12.4)') '   z range:     ', minz, maxz

      ! count materials
      materials = 0
      allocate (mattmp(cells))
      mattmp = material

      do i = 1, cells
         if (mattmp(i) /= -1) then
            materials = materials + 1
            do j = i+1, cells
               if (mattmp(j) == mattmp(i)) mattmp(j) = -1
            end do
         end if
      end do

      ! report materials data
      write (*,*) '  # materials: ', materials
      write (*,*) '  material #s:'

      do i = 1, cells
         if (mattmp(i) /= -1) write (*,*) '    ', mattmp(i)
      end do

      ! normalize materials to start at 1
      i = MINVAL(material)
      ! the "if" is a slight optimization
      if (i /= 1) material = material - i + 1

      deallocate (mattmp)

      analyzemesh = 0
      return
   end function analyzemesh

   !----------------------------------------------------------------------------

end module translate

!-------------------------------------------------------------------------------

program meshtranslate

   use translate

   implicit none

   integer :: status                    ! status return
   integer :: IARGC                     ! count of arguments on the command line
   integer :: imodecnt                  ! count of input modes specified
   logical :: pmode                     ! patran
   logical :: imode                     ! ideas
   logical :: amode                     ! avs
   integer :: omodecnt                  ! count of output modes specified
   logical :: gmode                     ! gmv
   logical :: emode                     ! ensight
   integer :: i                         ! loop index
   character (LEN=1024) :: mode         ! string to hold options
   character (LEN=1024) :: ifname       ! input file name
   character (LEN=1024) :: ofname       ! output file name

   !----------------------------------------------------------------------------

   ! initialize flags and counters
   imodecnt = 0
   pmode    = .false.
   imode    = .false.
   amode    = .false.
   omodecnt = 0
   gmode    = .false.
   emode    = .false.

   ! parse command line
   if (IARGC() < 4) then
      call usage ()
   end if

   do i = 1, IARGC()-2
      call GETARG (i, mode)

      select case (TRIM(mode))
      case ('-patran')
         pmode    = .true.
         imodecnt = imodecnt + 1
      case ('-ideas')
         imode    = .true.
         imodecnt = imodecnt + 1
      case ('-avs')
         amode    = .true.
         imodecnt = imodecnt + 1
      case ('-gmv')
         gmode    = .true.
         omodecnt = omodecnt + 1
      case ('-ensight')
         emode    = .true.
         omodecnt = omodecnt + 1
      case default
         write (*,*) 'invalid option '//TRIM(mode)
         call usage ()
      end select
   end do

   if (imodecnt /= 1) call usage ()
   if (omodecnt <  1) call usage ()

   call GETARG (IARGC()-1, ifname)
   call GETARG (IARGC(),   ofname)

   ! read input file
   if (pmode) status = readpatran   (TRIM(ifname))
   if (imode) status = readideas    (TRIM(ifname))
   if (amode) status = readavs      (TRIM(ifname))

   ! analyze mesh
   status = analyzemesh ()

   ! write output file
   if (gmode) status = writegmv     (TRIM(ofname))
   if (emode) status = writeensight (TRIM(ofname))

contains

   subroutine usage ()
      write (*,*) 'usage: meshtranslate -input -output [-output ...] infile outfileroot'
      write (*,*)
      write (*,*) '       one -input option is required from:'
      write (*,*) '              -patran'
      write (*,*) '              -ideas'
      write (*,*) '              -avs'
      write (*,*) '       at least one -output option is required from:'
      write (*,*) '              -gmv'
      write (*,*) '              -ensight'
      stop
   end subroutine usage

end program meshtranslate

!-------------------------------------------------------------------------------

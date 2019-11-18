!*******************************************************************************
! Copyright Notice
!  + 2010-2012 North Carolina State University
!  + 2010-2012 Pacific Northwest National Laboratory
! 
! This file is part of SCORPIO.
! 
! SCORPIO is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! any later version.
! 
! SCORPIO is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public License
! along with SCORPIO.  If not, see <http://www.gnu.org/licenses/>.
! 
!*******************************************************************************
program test
    implicit none
    include 'mpif.h'
    include 'scorpiof.h'

    integer :: i, ierr
    integer :: mpi_rank, mpi_size
    integer :: fid  ! file handler
    integer :: numgroups ! number of sub-groups to split the global communicator
    integer :: preferredGroupSize ! Each IO group size 
    integer :: mympi_comm
    integer :: mygid ! my IO group ID
    integer :: ndims
    integer,allocatable :: globaldims(:), localdims(:)
    real*8,allocatable :: dvector(:)
    real,allocatable :: svector(:)
    integer,allocatable :: vector(:)
    integer:: vecsize
    integer:: scenario
    integer, parameter :: IOGROUPS=1, ALLNODES=2, SMALLFILE=3
    character(len=80) :: argbuffer

    ! for passing to C, character strings must be terminated this way.
    character(len=80) :: infilename = 'input.h5'//CHAR(0)
    character(len=80) :: outfilename = 'output.h5'//CHAR(0)
    character(len=80) :: celldataset_name = '/Materials/Cell Ids'//CHAR(0)
    character(len=80) :: materialdataset_name = '/Materials/Material Ids'//CHAR(0)
    integer,allocatable :: cellids(:), materialids(:)
    real*8 :: tstart,tend, walltime, rtstart, rtend, wtstart, wtend, otstart, otend, otime, ctime, ctstart, ctend, itstart, itend

    !First argument is scenario number
    call getarg(1, argbuffer)
    read(argbuffer, *) scenario

    call MPI_INIT(ierr)
    tstart = MPI_Wtime()

    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
    mympi_comm = MPI_COMM_WORLD

    if ( scenario == IOGROUPS ) then
        !Second argument is number of IO groups ( only applicable when scenario is IOGROUPS) 
        call getarg(2, argbuffer)
        read(argbuffer, *) numgroups
        !numgroups = mpi_size/1024 ! 96 for 98034 cores, i.e., groups of 1024 cores

        if ( mpi_rank == 0) then
            write (*, '(A, I3,A)') "Scenario : ", scenario ," - IOGROUPS"
            write(*,'(A,I4)') "NumIOgroups = ", numgroups
        end if 
    else if ( scenario == ALLNODES ) then
        numgroups = mpi_size
        if ( mpi_rank == 0) then
            write (*, '(A, I3, A)') "Scenario : ", scenario ," - ALLNODES"
            write(*,'(A,I8)') "NumIOgroups = ", numgroups
        end if 
    else if ( scenario == SMALLFILE) then
        numgroups = 2
        if ( mpi_rank == 0) then
            write (*, '(A, I3, A)') "Scenario : ", scenario ," - SMALLFILE"
            write(*,'(A,I4)') "NumIOgroups = ", numgroups
        end if 
    end if

    write(*,'(A,I10)') "Got MPI Communicator = ", mympi_comm
    ! Set number of IO groups (OR) preferredGroupSize to initialize the parallel
    ! I/O library. Here, we are using numgroups.
    itstart = MPI_Wtime()
    call fscorpio_iogroup_init(numgroups, mympi_comm, mygid, ierr)
    ! call fscorpio_iogroup_init2(preferredGroupSize, numgroups, mympi_comm, mygid, ierr)
    itend = MPI_Wtime()

    otstart = MPI_Wtime()

    call fscorpio_open_file(infilename, mygid, SCORPIO_FILE_READONLY, fid, ierr)

    otend = MPI_Wtime()
    otime = otend - otstart

    call fscorpio_group_exists('/testPhantomGroup'//CHAR(0), fid, mygid, ierr)
    write (*,*) 'Group /testPhantomGroup exists: ', ierr
    call fscorpio_group_exists('Materials'//CHAR(0), fid, mygid, ierr)
    write (*,*) 'Group /Materials exists: ', ierr

    call fscorpio_dataset_exists('Materials'//CHAR(0), fid, mygid, ierr)
    write (*,*) 'Dataset Materials exists: ', ierr
    call fscorpio_dataset_exists('Materials/Cell Ids'//CHAR(0), fid, mygid, ierr)
    write (*,*) 'Dataset Materials/Cell Ids exists: ', ierr

    call fscorpio_get_dataset_ndims(ndims, fid, celldataset_name, mygid, ierr)
    !call fscorpio_get_dataset_ndims(ndims, fid, 'testdata'//CHAR(0), mygid, ierr)
    allocate(globaldims(ndims))
    allocate(localdims(ndims))

    call fscorpio_get_dataset_dims(globaldims, fid, celldataset_name, mygid, ierr)
    !call fscorpio_get_dataset_dims(globaldims, fid, 'testdata'//CHAR(0), mygid, ierr)

    write (*,*) 'Number of dimensions: ', ndims

    do i=1,ndims
        write (*,'(A,I1,A,I10)') "dims(", i, "): ", globaldims(i)
    enddo

    ! equally across all processes
    localdims(1) = globaldims(1)/mpi_size

    vecsize = 1
    do i=1,ndims
        vecsize = vecsize*localdims(i)
    enddo

    allocate(cellids(vecsize))
    allocate(materialids(vecsize))

    !No need to initialize as we are reading from input file and writing to output file
    !do i=1,vecsize
    !    cellids(i) = localdims(1)*mpi_rank + i
    !    materialids(i) = localdims(1)*mpi_rank + i
    !enddo

    rtstart = MPI_Wtime()
    if ( scenario == SMALLFILE) then
        vecsize = 1
        do i=1,ndims
            vecsize = vecsize*globaldims(i)
        enddo
        allocate(cellids(vecsize))
        call fscorpio_read_dataset( cellids, SCORPIO_INTEGER, ndims, globaldims, localdims, & 
                fid, celldataset_name, mygid, SCORPIO_EVERYONE_ENTIRE_DATASET_READ, ierr)
        do i=1,vecsize
            write(*, *) mpi_rank, cellids(i)
        enddo
    else
        call fscorpio_read_dataset( cellids, SCORPIO_INTEGER, ndims, globaldims, localdims, & 
            fid, celldataset_name, mygid, SCORPIO_NONUNIFORM_CONTIGUOUS_READ, ierr)
    end if

    call fscorpio_read_dataset( materialids, SCORPIO_INTEGER, ndims, globaldims, localdims, & 
            fid, materialdataset_name, mygid, SCORPIO_NONUNIFORM_CONTIGUOUS_READ, ierr)
    rtend = MPI_Wtime()

    ctstart = MPI_Wtime()
    call fscorpio_close_file( fid, mygid, ierr)
    ctend = MPI_Wtime()
    ctime = ctend - ctstart

    otstart = MPI_Wtime()
    call fscorpio_open_file(outfilename, mygid, SCORPIO_FILE_CREATE, fid, ierr)
    otend = MPI_Wtime()
    otime = otime + (otend - otstart)

    wtstart = MPI_Wtime()
    call fscorpio_write_dataset( cellids, SCORPIO_INTEGER, ndims, globaldims, localdims, & 
            fid, celldataset_name, mygid, SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE, ierr)

    call fscorpio_write_dataset( materialids, SCORPIO_INTEGER, ndims, globaldims, localdims, & 
            fid, materialdataset_name, mygid, SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE, ierr)
    wtend = MPI_Wtime()

    ctstart = MPI_Wtime()
    call fscorpio_close_file( fid, mygid, ierr)
    ctend = MPI_Wtime()
    ctime = ctime + (ctend - ctstart)

    call fscorpio_iogroup_cleanup( mygid, ierr)

    tend = MPI_Wtime()

    call MPI_FINALIZE(ierr)

    walltime = tend - tstart

    1000 format (A, F10.4, A)
    if ( mpi_rank == 0) then
        write (*,'(A,I8,A)') 'Each process read ', vecsize, ' elements each from two datasets.'
        write (*,1000) 'Init Time: ', (itend-itstart), ' seconds.' 
        write (*,1000) 'Open Time: ', otime, ' seconds.' 
        write (*,1000) 'Close Time: ', ctime, ' seconds.' 
        write (*,1000) 'Read Time: ', (rtend-rtstart), ' seconds.' 
        write (*,1000) 'Write Time: ', (wtend-wtstart), ' seconds.' 
        write (*,1000) 'Wallclock Time: ', walltime, ' seconds.' 
    end if

end program test

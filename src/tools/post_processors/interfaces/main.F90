PROGRAM INTERFACE_POLYGONS
  !=======================================================================
  ! Interfaces_polygon:
  !
  !   A simple post-processing program that reads an interface
  !   dump from TELLURIDE and writes it back out in the form
  !   of triangles.  The triangle output is compatible with either
  !   AVS UCD or GMV format.
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module,       only: ndim, nvc, nec
  use triangle_module,        only: INTERFACE_TRIANGLES
  use triangle_output_module, only: AVS_OUTPUT, GMV_OUTPUT,            &
                                    write_gmv, write_avs
  use user_input_module,      only: dump_end, dump_start, nicells,     &
                                    READ_DUMP, step, time, USER_STDIN, &
                                    int_file

  implicit none

  character(80) :: RCSid

  ! Global Variables
  real(r8), pointer, dimension(:)     :: Rho
  real(r8), pointer, dimension(:,:)   :: Normal
  real(r8), pointer, dimension(:,:,:) :: Xv
  logical :: ortho_mesh

  ! Global Plane arrays.
  integer,  pointer, dimension(:,:)   :: N_Order
  integer,  pointer, dimension(:)     :: Nvrt
  real(r8), pointer, dimension(:,:,:) :: Pv
  real(r8), pointer, dimension(:,:)   :: Mdpnt
  integer,  pointer, dimension(:,:)   :: Perm

  ! Local Variables
  integer :: i, j, k, lst

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Process the user input from stdin.
  call USER_STDIN ()

  ! Loop over dumps
  DUMP_LOOP: do k = 1,dump_end

     ! Read time, cycle, and number of mixed cells
     read (1) time
     read (1) step
     read (1) nicells

     if (k >= dump_start) then
        write(*,1) k
1       format(/,3x,'Processing dump ', i3,' . . .')
        write (*,2) time, step, nicells
2       format(3x,'Time : ',1pe10.2,/,3x,'Cycle: ',i4,/,3x,'Interface cells: ',i5)
     else
        write(*,3) k
3       format(/,3x,'Skipping dump ', i3,' . . .')
     end if

     ! Deallocate input arrays if they've been previously allocated.
     if (ASSOCIATED(Xv)) DEALLOCATE (Xv, Rho, Normal,N_Order, Nvrt, Pv, Perm, Mdpnt)

     ! Allocate the input arrays.
     ALLOCATE (Xv(nvc,ndim,nicells), Rho(nicells), Normal(ndim,nicells),   &
               N_Order(nec,nicells), Nvrt(nicells), Pv(ndim,nec,nicells),  &
               Perm(nicells,nec), Mdpnt(nicells,ndim))

     ! Read this dump.
     call READ_DUMP (k, Xv, Rho, Normal, ortho_mesh)

     ! Go back if we do not process this dump
     if (k < dump_start) cycle DUMP_LOOP

     ! Compute interface polygons (triangles) for this dump.
     call INTERFACE_TRIANGLES (nicells, Xv, Rho, Normal, N_Order, Nvrt, Pv, &
                               Perm, Mdpnt, ortho_mesh)

     ! Write out a UCD dump file readable by AVS
     if (write_avs) call AVS_OUTPUT (nicells, Xv, N_Order, Nvrt, Pv)

     ! Write out a dump file readable by GMV
     if (write_gmv) call GMV_OUTPUT (nicells, Xv, Normal, Nvrt, Pv, &
                           Perm, Mdpnt)

  end do DUMP_LOOP

  ! Close the input file.
  close (1)

  ! Done, inform user.
  write (*,4) TRIM(int_file)
4 format(/,' Finished processing file ',a,/)

    DEALLOCATE(Xv, Rho, Normal, N_Order, Pv, Nvrt, Mdpnt, Perm)

  stop

END PROGRAM INTERFACE_POLYGONS

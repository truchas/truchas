MODULE USER_INPUT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define variables and procedures necessary to process the user stdin.
  !
  ! Public Interface(s):
  !
  !   * call READ_DUMP (Xv, Rho, Normal, ortho_mesh)
  !
  !     Read an interface dump.
  !
  !   * call USER_STDIN ()
  !
  !     Process the user stdin.
  !
  ! Contains: READ_DUMP
  !           USER_STDIN
  !
  ! Author(s): Douglas B. Kothe (LANL Group T-3, dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module,  only: string_len
  implicit none
  private

  ! Public Procedures
  public :: READ_DUMP, USER_STDIN

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  character(string_len), public, save :: int_file
  integer,    public, save :: dump_start, dump_end,  &
                                               nicells, step
  real(r8),      public, save :: time

CONTAINS

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>

  SUBROUTINE READ_DUMP (k, Xv, Rho, Normal, ortho_mesh)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read in the interface data.
    !
    !=======================================================================
    use kind_module,            only: int_kind, real_kind, log_kind
    use parameter_module,       only: ndim, nvc
    use triangle_output_module, only: avs_file, gmv_file

    ! Arguments
    integer,                              intent(IN)  :: k
    real(r8), dimension(nicells),          intent(OUT) :: Rho
    real(r8), dimension(ndim,nicells),     intent(OUT) :: Normal
    real(r8), dimension(nvc,ndim,nicells), intent(OUT) :: Xv
    logical,                              intent(OUT) :: ortho_mesh

    ! Local Variables
    integer                          :: i, v, n
    real(r8), dimension(ndim,nicells) :: X

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

     ! Generate a polygon AVS output file name for this dump.
     avs_file = int_file
     i = LEN_TRIM(int_file)

     ! Append the dump number to the file name.
     if (k >= 0 .and. k < 10) then
        write (avs_file(i+1:), '(a4,i1)') '.000', k
     else if (k > 9 .and. k < 100) then
        write (avs_file(i+1:), '(a3,i2)') '.00', k
     else if (k > 99 .and. k < 1000) then
        write (avs_file(i+1:), '(a2,i3)') '.0', k
     end if

     ! Append a descriptive suffix.
     gmv_file = avs_file
     avs_file = TRIM(avs_file) // '.avs'
     gmv_file = TRIM(gmv_file) // '.gmv'

    ! Read ortho_mesh
       read (1) ortho_mesh

    ! Read the hex vertices
    HEX_VERTICES: do v = 1, nvc
         do n=1,ndim
       read (1) X(n,:)
       Xv(v,n,:) = X(n,:)
         end do
    end do HEX_VERTICES

    ! Read interface plane constant Rho.
    read (1) Rho

    ! Read the plane normal
     do n = 1,ndim
    read (1) Normal(n,:)
     end do

  END SUBROUTINE READ_DUMP

  SUBROUTINE USER_STDIN
    !=======================================================================
    ! Purpose(s):
    !
    !   Process the user stdin.
    !
    !=======================================================================
    use kind_module,            only: log_kind, int_kind
    use parameter_module,       only: string_len
    use triangle_output_module, only: write_avs, write_gmv, write_mesh

    ! Local Variables
    character(string_len) :: line
    logical    :: file_exist
    integer    :: i

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Get the dump file name
    GET_DUMP_FILE: do i = 1,HUGE(i/2)

       write (*, 1, ADVANCE = 'no')
1      format(/' Interface dump file name: ')
       read (*,'(a)') int_file
       int_file = ADJUSTL(int_file)
       inquire (FILE = int_file, EXIST = file_exist)
       if (file_exist) then
          exit GET_DUMP_FILE
       else
          write (*,2) TRIM(int_file)
2         format (/,' File ',a,' not found, please reenter file name')
       end if

    end do GET_DUMP_FILE

    ! Open the interface dump file
    open(1,FILE = int_file, FORM = 'unformatted', STATUS = 'old')
    write(*,3) TRIM(int_file)
3   format (/,' Opening interface dump file ',a)

    ! Get the dump number ranges to process
    DUMP_NUMBERS: do i = 1,HUGE(1/2)

       write(*,4, ADVANCE = 'no')
4      format(/,' Range of dump numbers to process (e.g. 13 25): ')
       read (*,*) dump_start, dump_end
       if (dump_start <= 0) then
          write (*,13) 
13        format (/,' Starting dump number must be positive!')
       else if (dump_end < dump_start) then
          write (*,12) 
12        format (/,' Starting dump number must be less than ending dump number!')
       else
          write (*,5) dump_start, dump_end
5         format (/' Generating interface polygons for dumps ', i3,' to ', i3)
          exit DUMP_NUMBERS
       end if

    end do DUMP_NUMBERS

    ! Check to see if interface cells are to be written out
    write_mesh = .false.
    OUTPUT_CELLS: do i = 1,HUGE(1/2)

       write (*,6, ADVANCE = 'no')
6      format(/' Write out interface cell vertices (y or n)? ')
       read(*,'(a)') line
       line = ADJUSTL(line)
       write_mesh = TRIM(line) == 'y' .or. TRIM(line) == 'yes'
       if ((.not. write_mesh) .and. &
          (TRIM(line) /= 'n' .and. TRIM(line) /= 'no')) then
          write (*,7) TRIM(line)
7         format (/,1x,a,' is an invalid option, try again')
       else
          exit OUTPUT_CELLS
       end if

    end do OUTPUT_CELLS

! let's skip AVS and assume GMV
!     ! Check to see if AVS dumps are to be written
!     write_avs = .false.
!     OUTPUT_AVS: do i = 1,HUGE(1/2)
!
!        write (*,8, ADVANCE = 'no')
! 8      format(/' Write out AVS UCD input files (y or n)? ')
!        read(*,'(a)') line
!        line = ADJUSTL(line)
!        write_avs = TRIM(line) == 'y' .or. TRIM(line) == 'yes'
!        if ((.not. write_avs) .and. &
!           (TRIM(line) /= 'n' .and. TRIM(line) /= 'no')) then
!           write (*,9) TRIM(line)
! 9         format (/,1x,a,' is an invalid option, try again')
!        else
!           exit OUTPUT_AVS
!        end if
!
!     end do OUTPUT_AVS
!
!     ! Check to see if GMV dumps are to be written
!     write_gmv = .false.
!     OUTPUT_GMV: do i = 1,HUGE(1/2)
!
!        write (*,10, ADVANCE = 'no')
! 10     format(/' Write out GMV input files (y or n)? ')
!        read(*,'(a)') line
!        line = ADJUSTL(line)
!        write_gmv = TRIM(line) == 'y' .or. TRIM(line) == 'yes'
!        if ((.not. write_gmv) .and. &
!           (TRIM(line) /= 'n' .and. TRIM(line) /= 'no')) then
!           write (*,11) TRIM(line)
! 11        format (/,1x,a,' is an invalid option, try again')
!        else
!           exit OUTPUT_GMV
!        end if
!
!    end do OUTPUT_GMV

    write_avs = .false.
    write_gmv = .true.

  END SUBROUTINE USER_STDIN

END MODULE USER_INPUT_MODULE

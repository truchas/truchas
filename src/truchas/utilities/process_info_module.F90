!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module process_info_module

  implicit none
  private

  public :: get_process_size, report_memory
  public :: mem_diag_open, mem_diag_write, mem_diag_close

  interface
    subroutine get_process_size (vsize, rsize, dsize) bind(c)
      use,intrinsic :: iso_c_binding, only: c_int
      integer(c_int), intent(out) :: vsize  ! virtual memory size (kB)
      integer(c_int), intent(out) :: rsize  ! resident memory size (kB)
      integer(c_int), intent(out) :: dsize  ! data memory size (kB)
    end subroutine
  end interface

  !! Memory diagnostics file data.
  logical, save :: mem_on = .false.
  integer, save :: mem_lun = -1

Contains

 !!
 !! REPORT_MEMORY -- Logs a summary of process memory usage.
 !!
 !! NNC, Apr 2012: It's questionable how accurate the memory size info is.
 !! I suspect that the data size is more relevant than the virtual size,
 !! but I won't make that change now. Morever, now having multiple meshes,
 !! the words-per-cell metric using the old mesh is losing its relevance.
 !!

  subroutine report_memory

    use,intrinsic :: iso_c_binding, only: c_int
    use parallel_communication, only: nPE, global_all, global_minval, global_maxval, global_sum
    use legacy_mesh_api, only: ncells_tot
    use constants_module, only: FLOATBYTES  ! we ought to compute this here
    use truchas_logging_services

    integer(c_int) :: vsize, rsize, dsize
    character(128) :: string(4)

    call get_process_size (vsize, rsize, dsize)
    if (global_all(vsize /= 0)) then
      if (npe == 1) then
        write(string,1) real(vsize)/1024, int(real(1024*vsize)/(ncells_tot*FLOATBYTES))
        call TLS_info (string(:2))
      else
        write(string,2) global_minval(real(vsize)/1024), global_maxval(real(vsize)/1024), &
            global_sum(real(vsize))/1024, int(global_sum(real(1024*vsize))/(ncells_tot*FLOATBYTES))
        call TLS_info (string(:4))
      end if
    end if

    1 format(9x,'         Process virtual memory used: ',es8.2,' mB',/,&
             9x,'                          words/cell: ',i0)
    2 format(9x,'Smallest process virtual memory used: ',es8.2,' mB',/,&
             9x,' Largest process virtual memory used: ',es8.2,' mB',/,&
             9x,'           Total virtual memory used: ',es8.2,' mB',/,&
             9x,'                          words/cell: ',i0)

  end subroutine report_memory

 !!
 !! MEM_DIAG_OPEN
 !!
 !! Opens the .mem  memory diagnostics output file in the output directory.
 !! Returns the logical argument AVAIL which indicates the availability of
 !! process memory usage information.  If not available, a message to that
 !! effect is written to the file and the file closed; no other MEM_DIAG_*
 !! calls should be made in this case.  This is a collective procedure.
 !!

  subroutine mem_diag_open

    use,intrinsic :: iso_c_binding, only: c_int
    use parallel_communication, only: is_IOP, global_all
    use truchas_env, only: output_file_name

    integer(c_int) :: vsize, rsize, dsize

    !! Check to see if memory usage info is available.
    call get_process_size (vsize, rsize, dsize)
    mem_on = global_all(vsize /= 0)
    !! Open the memory diagnostics file and write header (or regrets).
    if (is_IOP) then
      open(newunit=mem_lun,file=output_file_name('mem'),action='write',status='replace')
      if (mem_on) then
        write(mem_lun,'(a)') '===== MEMORY DIAGNOSTICS ' // repeat('=',55)
        write(mem_lun,'(a)') ' vsize[N]: total virtual memory (mB) used by process N'
        write(mem_lun,'(a)') ' rsize[N]: physical memory (mB) used by process N (resident size)'
        write(mem_lun,'(a)') ' dsize[N]: data memory (mB) used by process N'
        write(mem_lun,'(a)') repeat('=',80)
      else
        write(mem_lun,'(a)') 'Memory diagnostics are not available on this platform.'
        close(mem_lun)
        mem_lun = -1
      end if
    else
      mem_lun = -1
    endif

  end subroutine mem_diag_open

 !!
 !! MEM_DIAG_CLOSE -- closes the memory diagnostics file.
 !!

  subroutine mem_diag_close
    mem_on = .false.
    if (mem_lun /= -1) then
      close(mem_lun)
      mem_lun = -1
    end if
  end subroutine mem_diag_close

 !!
 !! MEM_DIAG_WRITE -- writes the current memory usage info to the memory
 !! diagnostics file.  This is a collective procedure.
 !!

  subroutine mem_diag_write (comment)

    use,intrinsic :: iso_c_binding, only: c_int
    use parallel_communication, only: nPE, is_IOP, gather

    character(*), intent(in) :: comment

    integer :: n
    integer(c_int) :: vsize, rsize, dsize
    real :: vsize_all(nPE), rsize_all(nPE), dsize_all(nPE)

    if (mem_on) then
      !! Gather memory usage from each process.
      call get_process_size (vsize, rsize, dsize)
      call gather (real(vsize)/1024, vsize_all)
      call gather (real(rsize)/1024, rsize_all)
      call gather (real(dsize)/1024, dsize_all)
      !! Write it from the I/O process.
      If (is_IOP) then
        write(mem_lun,'(/,a)') comment
        do n = 1, nPE
          write(mem_lun,'(2x,3(a,i0,a,es9.3,:,a))') &
              'vsize[', n, '] = ', vsize_all(n), '; ', &
              'rsize[', n, '] = ', rsize_all(n), '; ', &
              'dsize[', n, '] = ', dsize_all(n)
        end do
      end if
    end if

  end subroutine mem_diag_write

end module process_info_module

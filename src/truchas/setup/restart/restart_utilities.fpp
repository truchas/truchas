!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! include file for restart_utilities.F90

#ifdef _TYPE_
#undef _TYPE_
#endif

#ifdef _ARRAY_PROCS_
#undef _ARRAY_PROCS_
#endif

#ifdef _INT8_DATA_
#undef _INT8_DATA_
#define _TYPE_ integer(int8)
#define _READ_VAR0_STAT_ read_var_stat_int8_0
#define _READ_VAR0_HALT_ read_var_halt_int8_0
#define _READ_VAR1_STAT_ read_var_stat_int8_1
#define _READ_VAR1_HALT_ read_var_halt_int8_1
#define _ARRAY_PROCS_
#define _READ_ARRAY1_STAT_ read_array_stat_int8_1
#define _READ_ARRAY1_HALT_ read_array_halt_int8_1
#define _READ_ARRAY2_STAT_ read_array_stat_int8_2
#define _READ_ARRAY2_HALT_ read_array_halt_int8_2
#define _READ_ARRAY3_STAT_ read_array_stat_int8_3
#define _READ_ARRAY3_HALT_ read_array_halt_int8_3
#endif

#ifdef _INTEGER_DATA_
#undef _INTEGER_DATA_
#define _TYPE_ integer
#define _READ_VAR0_STAT_ read_var_stat_I0
#define _READ_VAR0_HALT_ read_var_halt_I0
#define _READ_VAR1_STAT_ read_var_stat_I1
#define _READ_VAR1_HALT_ read_var_halt_I1
#define _ARRAY_PROCS_
#define _READ_ARRAY1_STAT_ read_array_stat_I1
#define _READ_ARRAY1_HALT_ read_array_halt_I1
#define _READ_ARRAY2_STAT_ read_array_stat_I2
#define _READ_ARRAY2_HALT_ read_array_halt_I2
#define _READ_ARRAY3_STAT_ read_array_stat_I3
#define _READ_ARRAY3_HALT_ read_array_halt_I3
#endif

#ifdef _REAL_DATA_
#undef _REAL_DATA_
#define _TYPE_ real(r8)
#define _READ_VAR0_STAT_ read_var_stat_R0
#define _READ_VAR0_HALT_ read_var_halt_R0
#define _READ_VAR1_STAT_ read_var_stat_R1
#define _READ_VAR1_HALT_ read_var_halt_R1
#define _ARRAY_PROCS_
#define _READ_ARRAY1_STAT_ read_array_stat_R1
#define _READ_ARRAY1_HALT_ read_array_halt_R1
#define _READ_ARRAY2_STAT_ read_array_stat_R2
#define _READ_ARRAY2_HALT_ read_array_halt_R2
#define _READ_ARRAY3_STAT_ read_array_stat_R3
#define _READ_ARRAY3_HALT_ read_array_halt_R3
#endif

#ifdef _STRING_DATA_
#undef _STRING_DATA_
#define _TYPE_ character(len=*)
#define _READ_VAR0_STAT_ read_var_stat_S0
#define _READ_VAR0_HALT_ read_var_halt_S0
#define _READ_VAR1_STAT_ read_var_stat_S1
#define _READ_VAR1_HALT_ read_var_halt_S1
#endif

  subroutine _READ_VAR0_STAT_ (unit, var, stat)
    use parallel_communication, only: is_IOP
    use pgslib_module, only:  pgslib_bcast
    integer, intent(in) :: unit
    _TYPE_, intent(out) :: var
    integer, intent(out) :: stat
    if (is_IOP) read(unit,iostat=stat) var
    call pgslib_bcast (stat)
    if (stat /= 0) return
    call pgslib_bcast (var)
  end subroutine _READ_VAR0_STAT_

  subroutine _READ_VAR0_HALT_ (unit, var, errmsg)
    integer, intent(in) :: unit
    _TYPE_, intent(out) :: var
    character(len=*), intent(in) :: errmsg
    integer :: stat
    call _READ_VAR0_STAT_ (unit, var, stat)
    if (stat /= 0) call halt (errmsg)
  end subroutine _READ_VAR0_HALT_

  subroutine _READ_VAR1_STAT_ (unit, var, stat)
    use parallel_communication, only: is_IOP
    use pgslib_module, only:  pgslib_bcast
    integer, intent(in) :: unit
    _TYPE_, intent(out) :: var(:)
    integer, intent(out) :: stat
    if (is_IOP) read(unit,iostat=stat) var
    call pgslib_bcast (stat)
    if (stat /= 0) return
    call pgslib_bcast (var)
  end subroutine _READ_VAR1_STAT_

  subroutine _READ_VAR1_HALT_ (unit, var, errmsg)
    integer, intent(in) :: unit
    _TYPE_, intent(out) :: var(:)
    character(len=*), intent(in) :: errmsg
    integer :: stat
    call _READ_VAR1_STAT_ (unit, var, stat)
    if (stat /= 0) call halt (errmsg)
  end subroutine _READ_VAR1_HALT_

#ifdef _ARRAY_PROCS_
#undef _ARRAY_PROCS_
  subroutine _READ_ARRAY1_STAT_ (unit, array, perm, stat)

    use parallel_communication, only: is_IOP
    use pgslib_module, only:  pgslib_bcast, pgslib_dist, pgslib_collate, pgslib_global_sum, pgslib_global_all
    use permutations

    integer, intent(in) :: unit
    _TYPE_, intent(out) :: array(:)
    integer, intent(in), optional :: perm(:)
    integer, intent(out) :: stat

    logical :: coll_perm
    integer :: n
    integer, allocatable :: g_perm(:)
    _TYPE_, allocatable :: g_sect(:)

    n = pgslib_global_sum(size(array,1))
    if (.not.is_IOP) n = 0
    allocate(g_sect(n))

    if (present(perm)) then
      coll_perm = pgslib_global_all(is_IOP .or. size(perm) == 0)
      if (coll_perm) then ! passed a collated permutation; use it.
        ASSERT( size(perm) == n )
        ASSERT( is_perm(perm) .or. n == 0 )
        coll_perm = .true.
      else  ! passed a distributed permutation; collate it.
        ASSERT( size(array,1) == size(perm) )
        allocate(g_perm(n))
        call pgslib_collate (g_perm, perm)
        ASSERT( is_perm(g_perm) .or. n == 0 )
        coll_perm = .false.
      end if
    end if

    if (is_IOP) read(unit,iostat=stat) g_sect
    call pgslib_bcast (stat)
    if (stat == 0) then
      if (is_IOP .and. present(perm)) then
        if (coll_perm) then
          call reorder (g_sect, perm)
        else
          call reorder (g_sect, g_perm)
        end if
      end if
      call pgslib_dist (array, g_sect)
    end if
    deallocate(g_sect)
    if (allocated(g_perm)) deallocate(g_perm)

  end subroutine _READ_ARRAY1_STAT_


  subroutine _READ_ARRAY1_HALT_ (unit, array, perm, errmsg)
    use string_utilities, only: i_to_c
    integer, intent(in) :: unit
    _TYPE_, intent(out) :: array(:)
    integer, intent(in), optional :: perm(:)
    character(len=*), intent(in) :: errmsg
    integer :: stat
    call _READ_ARRAY1_STAT_ (unit, array, perm, stat=stat)
    if (stat /= 0) call halt (errmsg // ': iostat=' // i_to_c(stat))
  end subroutine _READ_ARRAY1_HALT_


  subroutine _READ_ARRAY2_STAT_ (unit, array, perm, stat)

    use parallel_communication, only: is_IOP
    use pgslib_module, only:  pgslib_bcast, pgslib_dist, pgslib_collate, pgslib_global_sum, pgslib_global_all
    use permutations

    integer, intent(in) :: unit
    _TYPE_, intent(out) :: array(:,:)
    integer, intent(in), optional :: perm(:)
    integer, intent(out) :: stat

    logical :: coll_perm
    integer :: n, i
    integer, allocatable :: g_perm(:)
    _TYPE_, allocatable :: g_sect(:)

    n = pgslib_global_sum(size(array,2))
    if (.not.is_IOP) n = 0
    allocate(g_sect(n))

    if (present(perm)) then
      coll_perm = pgslib_global_all(is_IOP .or. size(perm) == 0)
      if (coll_perm) then ! passed a collated permutation; use it.
        ASSERT( size(perm) == n )
        ASSERT( is_perm(perm) .or. n == 0 )
        coll_perm = .true.
      else  ! passed a distributed permutation; collate it.
        ASSERT( size(array,2) == size(perm) )
        allocate(g_perm(n))
        call pgslib_collate (g_perm, perm)
        ASSERT( is_perm(g_perm) .or. n == 0 )
        coll_perm = .false.
      end if
    end if

    INPUT: do i = 1, size(array,1)
      if (is_IOP) read(unit,iostat=stat) g_sect
      call pgslib_bcast (stat)
      if (stat /= 0) exit INPUT
      if (is_IOP .and. present(perm)) then
        if (coll_perm) then
          call reorder (g_sect, perm)
        else
          call reorder (g_sect, g_perm)
        end if
      end if
      call pgslib_dist (array(i,:), g_sect)
    end do INPUT
    deallocate(g_sect)
    if (allocated(g_perm)) deallocate(g_perm)

  end subroutine _READ_ARRAY2_STAT_


  subroutine _READ_ARRAY2_HALT_ (unit, array, perm, errmsg)
    use string_utilities, only: i_to_c
    integer, intent(in) :: unit
    _TYPE_, intent(out) :: array(:,:)
    integer, intent(in), optional :: perm(:)
    character(len=*), intent(in) :: errmsg
    integer :: stat
    call _READ_ARRAY2_STAT_ (unit, array, perm, stat=stat)
    if (stat /= 0) call halt (errmsg // ': iostat=' // i_to_c(stat))
  end subroutine _READ_ARRAY2_HALT_


  subroutine _READ_ARRAY3_STAT_ (unit, array, perm, stat)

    use parallel_communication, only: is_IOP
    use pgslib_module, only:  pgslib_bcast, pgslib_dist, pgslib_collate, pgslib_global_sum, pgslib_global_all
    use permutations

    integer, intent(in) :: unit
    _TYPE_, intent(out) :: array(:,:,:)
    integer, intent(in), optional :: perm(:)
    integer, intent(out) :: stat

    logical :: coll_perm
    integer :: n, i, j
    integer, allocatable :: g_perm(:)
    _TYPE_, allocatable :: g_sect(:)

    n = pgslib_global_sum(size(array,3))
    if (.not.is_IOP) n = 0
    allocate(g_sect(n))

    if (present(perm)) then
      coll_perm = pgslib_global_all(is_IOP .or. size(perm) == 0)
      if (coll_perm) then ! passed a collated permutation; use it.
        ASSERT( size(perm) == n )
        ASSERT( is_perm(perm) .or. n == 0 )
        coll_perm = .true.
      else  ! passed a distributed permutation; collate it.
        ASSERT( size(array,3) == size(perm) )
        allocate(g_perm(n))
        call pgslib_collate (g_perm, perm)
        ASSERT( is_perm(g_perm) .or. n == 0 )
        coll_perm = .false.
      end if
    end if

    INPUT: do j = 1, size(array,2)
      do i = 1, size(array,1)
        if (is_IOP) read(unit,iostat=stat) g_sect
        call pgslib_bcast (stat)
        if (stat /= 0) exit INPUT
        if (is_IOP .and. present(perm)) then
          if (coll_perm) then
            call reorder (g_sect, perm)
          else
            call reorder (g_sect, g_perm)
          end if
        end if
        call pgslib_dist (array(i,j,:), g_sect)
      end do
    end do INPUT
    deallocate(g_sect)
    if (allocated(g_perm)) deallocate(g_perm)

  end subroutine _READ_ARRAY3_STAT_


  subroutine _READ_ARRAY3_HALT_ (unit, array, perm, errmsg)
    use string_utilities, only: i_to_c
    integer, intent(in) :: unit
    _TYPE_, intent(out) :: array(:,:,:)
    integer, intent(in), optional :: perm(:)
    character(len=*), intent(in) :: errmsg
    integer :: stat
    call _READ_ARRAY3_STAT_ (unit, array, perm, stat=stat)
    if (stat /= 0) call halt (errmsg // ': iostat=' // i_to_c(stat))
  end subroutine _READ_ARRAY3_HALT_

#undef _READ_ARRAY1_STAT_
#undef _READ_ARRAY1_HALT_
#undef _READ_ARRAY2_STAT_
#undef _READ_ARRAY2_HALT_
#undef _READ_ARRAY3_STAT_
#undef _READ_ARRAY3_HALT_
#endif

#undef _TYPE_
#undef _READ_SCALAR_STAT_
#undef _READ_SCALAR_HALT_

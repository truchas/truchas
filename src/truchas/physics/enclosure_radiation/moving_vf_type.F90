!!
!! MOVING_VF_TYPE
!!
!! A concrete implementation of the ENCL_VF base class that implements time
!! dependent enclosure view factors that are defined by interpolation on a
!! time sequence of STATIC_VF class objects.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! June 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module moving_vf_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use encl_vf_class
  use static_vf_class
  use parallel_communication
  implicit none
  private

  type, extends(encl_vf), public :: moving_vf
    private
    integer :: n ! current interval index [times(n), times(n+1)]
    real(r8), allocatable, public :: times(:)
    character(:), allocatable :: filenames(:)
    real(r8) :: w1, w2
    class(static_vf), allocatable :: vf1, vf2
    real(r8) :: tlast = -huge(1.0d0)
    ! In parent class:
    ! integer, public :: nface, nface_tot
    ! integer, allocatable :: face_map(:)
    ! real, allocatable, public :: amb_vf(:)
    ! logical :: has_ambient
  contains
    procedure :: init
    procedure :: matvec
    procedure :: set_time
    procedure :: next_interval
    procedure :: filename
  end type

contains

  subroutine init(this, tinit, times, filenames)

    class(moving_vf), intent(out) :: this
    real(r8), intent(in) :: tinit
    real(r8), intent(in) :: times(:)
    character(*), intent(in) :: filenames(:)

    integer :: n
    logical :: flag

    ASSERT(size(times) == size(filenames))
    ASSERT(size(times) > 1)

    this%times = times
    this%filenames = filenames

    !! Verify the TIMES values are strictly increasing
    do n = 1, size(times) - 1
      if (times(n) >= times(n+1)) exit
    end do
    INSIST(n == size(times))

    !! Verify the named files exist
    if (is_IOP) then
      do n = 1, size(times)
        inquire(file=filenames(n), exist=flag)
        if (.not.flag) exit
      end do
    end if
    call broadcast(n)
    ASSERT(n > size(times))

    !! Find N such that TIMES(N) <= TINIT < TIMES(N+1)
    do n = 0, size(times)-1
      if (tinit < times(n+1)) exit
    end do

    !! Load the static view factors at the interval endpoints
    if (n > 0) call alloc_const_vf(this%vf1, filenames(n))
    if (n < size(times)) call alloc_const_vf(this%vf2, filenames(n+1))
    this%n = n

    !! Inherit the characteristics of one of the static VF
    if (allocated(this%vf1)) then
      call copy_base_components(src=this%vf1, dest=this)
    else
      call copy_base_components(src=this%vf2, dest=this)
    end if

    !! Ensure the other static VF has the same characteristics
    if (allocated(this%vf1) .and. allocated(this%vf2)) then
      INSIST(compatible_vf(this%vf1, this%vf2))
    end if

    call set_time(this, tinit)

  end subroutine init

  !! Allocate and initialize the proper STATIC_VF type
  subroutine alloc_const_vf(vf, filename)
    use facet_vf_type
    use patch_vf_type
    use rad_encl_file_type
    use truchas_logging_services
    class(static_vf), allocatable, intent(out) :: vf
    character(*), intent(in) :: filename
    type(rad_encl_file) :: file
    logical :: has_patches
    call TLS_info('    reading enclosure radiation view factors from ' // filename)
    if (is_IOP) call file%open_ro(filename)
    if (is_IOP) has_patches = file%has_patches()
    call broadcast(has_patches)
    if (has_patches) then
      allocate(patch_vf :: vf)
    else
      allocate(facet_vf :: vf)
    end if
    call vf%init(file)
    if (is_IOP) call file%close
  end subroutine

  !! Copy the base ENCL_VF components from SRC to DEST.
  subroutine copy_base_components(src, dest)
    class(encl_vf), intent(in) :: src
    class(encl_vf), intent(inout) :: dest
    dest%nface = src%nface
    dest%nface_tot = src%nface_tot
    dest%face_map = src%face_map
    dest%has_ambient = src%has_ambient
    if (allocated(src%amb_vf)) allocate(dest%amb_vf, mold=src%amb_vf)
  end subroutine

  !! Returns true if the ENCL_VF class objects are compatible
  logical function compatible_vf(a, b)
    class(encl_vf), intent(in) :: a, b
    compatible_vf = .false.
    if (a%nface_tot /= b%nface_tot) return
    if (global_any(a%nface /= b%nface)) return
    if (global_any(a%face_map /= b%face_map)) return
    if (a%has_ambient .neqv. b%has_ambient) return
    compatible_vf = .true.
  end function

  !! Computes the matrix-vector product Y = PHI*X.
  subroutine matvec(this, x, y)
    class(moving_vf), intent(in) :: this
    real(r8), intent(in)  :: x(:)
    real(r8), intent(out) :: y(:)
    real(r8) :: z(this%nface)
    if (this%w1 == 0) then
      call this%vf2%matvec(x, y)
    else if (this%w2 == 0) then
      call this%vf1%matvec(x, y)
    else
      call this%vf1%matvec(x, y)
      y = this%w1 * y
      call this%vf2%matvec(x, z)
      y = y + this%w2 * z
    end if
  end subroutine

  !! Advance view factor data to the next time interval.

  subroutine next_interval(this)
    class(moving_vf), intent(inout) :: this
    INSIST(this%n < size(this%times))
    call move_alloc(this%vf2, this%vf1)
    this%n = this%n + 1
    if (this%n == size(this%times)) return  ! reached the end
    call alloc_const_vf(this%vf2, this%filenames(this%n+1))
    INSIST(compatible_vf(this, this%vf2))
  end subroutine

  !! Set the time T for which MATVEC and AMB_VF will provide values.
  !! The specified time must belong to the current time interval.

  subroutine set_time(this, t)
    class(moving_vf), intent(inout) :: this
    real(r8), intent(in) :: t
    if (t == this%tlast) return ! nothing to do
    if (this%n == 0) then
      INSIST(t <= this%times(1))
      this%w1 = 0
      this%w2 = 1
      this%amb_vf = this%vf2%amb_vf
    else if (this%n == size(this%times)) then
      INSIST(t >= this%times(this%n))
      this%w1 = 1
      this%w2 = 0
      this%amb_vf = this%vf1%amb_vf
    else
      INSIST(t >= this%times(this%n))
      INSIST(t <= this%times(this%n+1))
      this%w1 = (this%times(this%n+1) - t) / (this%times(this%n+1) - this%times(this%n))
      this%w2 = (t - this%times(this%n))   / (this%times(this%n+1) - this%times(this%n))
      this%amb_vf = this%w1 * this%vf1%amb_vf + this%w2 * this%vf2%amb_vf
    end if
    this%tlast = t
  end subroutine

  !! Returns the the filename for one of the current static view factors

  function filename(this)
    class(moving_vf), intent(in) :: this
    character(:), allocatable :: filename
    if (this%n > 0) then
      filename = this%filenames(this%n)
    else
      filename = this%filenames(1)
    end if
  end function

end module moving_vf_type

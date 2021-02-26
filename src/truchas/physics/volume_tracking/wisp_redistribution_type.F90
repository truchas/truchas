!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

! Encapsulates a volume redistribution strategy into a class.  If we decide to
! try out multiple kinds of redistributions, this should be refactored in an
! abstract base class with multiple implementations

module wisp_redistribution_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use truchas_logging_services
  use truchas_timers
  use unstr_mesh_type
  use index_partitioning
  use parallel_communication
  implicit none
  private

  integer, parameter :: &
      remove_all = 0, &
      remove_some = 1, &
      remove_none = 2


  type, public :: wisp_redistribution
    private
    type(unstr_mesh), pointer :: mesh
    integer :: nrealfluid
    logical :: do_redistribution ! use conservative redistribution algortihm to remove wisps
    real(r8) :: wisp_cutoff ! consider fluid volume fractions < wisp_cutoff as wisps
    real(r8) :: wisp_absorption_fraction ! limit amount of wisp material a cell can absorb
    real(r8) :: wisp_neighbor_cutoff !
    real(r8), allocatable :: wisp_donor_volumes(:,:)
    real(r8), allocatable :: wisp_acceptor_volumes(:,:)
    real(r8), allocatable :: wisp_acceptor_fractions(:)
    real(r8), allocatable :: total_wisp_donor_v(:)
    real(r8), allocatable :: total_wisp_acceptor_v(:)
    real(r8), allocatable :: local_wisp_donor_v(:)
    ! needed for statistics
    integer, allocatable :: donor_n(:)
    integer, allocatable :: acceptor_n(:)
    integer, allocatable :: uncorrected_n(:)
  contains
    procedure :: init
    procedure :: redistribute
    procedure :: donor_volumes_view
    procedure :: acceptor_fractions_view
    procedure, private :: wisp_volumes
    procedure, private :: check_acceptor
    procedure, private :: sort_donors
    procedure, private :: remove_wisp_donors
    procedure, private :: fill_wisp_acceptors
    procedure, private :: compute_fractional_field
    procedure, private :: statistics
  end type wisp_redistribution

contains

  function donor_volumes_view(this) result(p)

    class(wisp_redistribution), intent(in), target :: this
    real(r8), pointer :: p(:,:)

    if (this%do_redistribution) then
      p => this%wisp_donor_volumes
    else
      p => null()
    endif
  end function donor_volumes_view

  function acceptor_fractions_view(this) result(p)

    class(wisp_redistribution), intent(in), target :: this
    real(r8), pointer :: p(:)

    if (this%do_redistribution) then
      p => this%wisp_acceptor_fractions
    else
      p => null()
    endif
  end function acceptor_fractions_view

  subroutine init(this, mesh, nrealfluid, params)

    use parameter_list_type

    class(wisp_redistribution), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(in) :: nrealfluid
    type(parameter_list), intent(inout) :: params

    this%mesh => mesh
    this%nrealfluid = nrealfluid

    ! wisp redistribution
    call params%get('wisp_redistribution', this%do_redistribution, default=.false.)
    call params%get('wisp_cutoff', this%wisp_cutoff, default=0.05_r8)
    call params%get('wisp_absorption_fraction', this%wisp_absorption_fraction, default=0.1_r8)
    call params%get('wisp_neighbor_cutoff', this%wisp_neighbor_cutoff, default=0.25_r8)

    if (this%do_redistribution) then
      allocate(this%wisp_donor_volumes(nrealfluid,this%mesh%ncell_onP))
      allocate(this%wisp_acceptor_volumes(nrealfluid,this%mesh%ncell_onP))
      allocate(this%wisp_acceptor_fractions(this%mesh%ncell_onP))
      allocate(this%total_wisp_donor_v(nrealfluid))
      allocate(this%total_wisp_acceptor_v(nrealfluid))
      allocate(this%local_wisp_donor_v(nrealfluid))
      allocate(this%donor_n(nrealfluid))
      allocate(this%acceptor_n(nrealfluid))
      allocate(this%uncorrected_n(nrealfluid))
    end if
  end subroutine init

  subroutine redistribute(this, vof, flux_vol, void)
    class(wisp_redistribution), intent(inout) :: this
    real(r8), intent(inout) :: vof(:,:), flux_vol(:,:)
    integer, intent(in) ::  void
    !-
    real(r8) :: mean_donor_volumes(this%nrealfluid)
    logical :: enough_acceptor(this%nrealfluid)
    integer :: remove_type(this%nrealfluid)
    real(r8) :: remove_volume(this%nrealfluid) ! local wisp volume to remove


    if (.not.this%do_redistribution .or. void /= 1 .or. this%nrealfluid == 0) return ! nothing to do

    call start_timer('redistribution')

    ! Step 1: compute the local donor/acceptor volumes
    call this%wisp_volumes(vof, mean_donor_volumes)

    ! Step 2: determine if we have enough acceptor volume to handle donor volumes
    call this%check_acceptor(enough_acceptor)

    ! Step 3: determine which processors remove (1) all, (2) some, or (3) none of
    !         their donor volume by sorting according to mean_volumes
    call this%sort_donors(mean_donor_volumes, enough_acceptor, remove_type, remove_volume)

    ! Step 4: remove local wisp donors from vof/flux_vol
    call this%remove_wisp_donors(enough_acceptor, remove_type, remove_volume, vof, flux_vol)

    ! Step 5: add wisp acceptor volumes to vof/flux_vol
    call this%fill_wisp_acceptors(vof, flux_vol)

    ! Step 6: redistributing scalar quantities associated with the wisps is
    !         easier if we compresss acceptor_volumes into a fractional field
    call this%compute_fractional_field()

    ! communicate results... I don't think we need to communicate flux_vol changes since
    ! they are strictly local.  May need to be revistied if odd bugs show up
    call gather_boundary(this%mesh%cell_ip, vof)

    call this%statistics()

    call stop_timer('redistribution')

  end subroutine redistribute


  !---
  ! fill out `wisp_donor_volumes` and `wisp_acceptor_volumes` on each process.
  ! compute `mean_donor_volumes` in the same loops.  mean volumes will only be needed
  ! in the event that there is not enough acceptor to go around
  !
  ! global communication to fill out `total_wisp_XXX` arrays
  !---
  subroutine wisp_volumes(this, vof, mean_donor_volumes)
    class(wisp_redistribution), intent(inout) :: this
    real(r8), intent(inout) :: vof(:,:)
    real(r8), intent(out) :: mean_donor_volumes(:)
    !-
    integer :: i, m, void_i
    real(r8) :: neighbor_vof

    this%wisp_donor_volumes = 0.0_r8
    this%wisp_acceptor_volumes = 0.0_r8
    mean_donor_volumes = 0.0_r8
    this%local_wisp_donor_v = 0.0_r8
    this%donor_n = 0
    this%acceptor_n = 0
    void_i = this%nrealfluid + 1

    do i = 1, this%mesh%ncell_onP

      if (vof(void_i,i) > 0.0_r8) then

        associate (cn => this%mesh%cnhbr(this%mesh%xcnhbr(i):this%mesh%xcnhbr(i+1)-1))

          do m = 1, this%nrealfluid

            ! check donor cell criteria. Check for > 0 so we don't bias the mean
            if (vof(m, i) > 0.0_r8 .and. vof(m, i) <= this%wisp_cutoff) then
              ! check neighbor criteria - average vof of neighbors
              neighbor_vof = sum(vof(m, pack(cn, cn > 0))) / size(pack(cn, cn > 0))
              ! check neighbor criteria - max vof of neighbors
              ! neighbor_vof = maxval(vof(m, pack(cn, cn > 0)))
              if (neighbor_vof < this%wisp_neighbor_cutoff) then
                this%wisp_donor_volumes(m, i) = vof(m, i) * this%mesh%volume(i)
                this%donor_n(m) = this%donor_n(m) + 1
                this%local_wisp_donor_v(m) = this%local_wisp_donor_v(m) + this%wisp_donor_volumes(m, i)
              end if
            end if

            ! check acceptor cell criteria
            if (vof(m, i) > this%wisp_neighbor_cutoff .and. vof(void_i, i) > 0.0_r8) then
              if (any(vof(m, pack(cn, cn > 0)) > this%wisp_neighbor_cutoff)) then
                this%wisp_acceptor_volumes(m, i) = this%mesh%volume(i) * &
                    min(vof(void_i, i), vof(m, i) * this%wisp_absorption_fraction)
                this%acceptor_n(m) = this%acceptor_n(m) + 1
              end if
            end if
          end do

        end associate
      end if
    end do
    ! write(*,*) "VOF[:,5048]: ", vof(:,5048)
    ! write(*,*) "WISP DONOR VOLUME[:,5048]: ", this%wisp_donor_volumes(:,5048)/this%mesh%volume(5048)
    ! write(*,*) "WISP ACCEPTOR VOLUME[:,5048]: ", this%wisp_acceptor_volumes(:,5048)

    ! At this point, the sum of acceptor volumes may push the fluid vof out of bounds in
    ! a cell with multiple acceptor fluids.
    if (this%nrealfluid > 1) then

      do i = 1, this%mesh%ncell_onP

        if (count(this%wisp_acceptor_volumes(:, i) > 0.0_r8) > 1) then
          block
            real(r8) :: excess, total

            total = sum(this%wisp_acceptor_volumes(:, i))
            excess = total - (vof(void_i, i) * this%mesh%volume(i))

            if (excess > 0.0_r8) then
              ! proportionally shrink all acceptor volumes
              do m = 1, this%nrealfluid
                this%wisp_acceptor_volumes(m, i) = &
                    this%wisp_acceptor_volumes(m, i) * (1.0_r8 - excess / total)
              end do
            end if
          end block
        end if

      end do
    end if

    ! finish mean calculation
    do m = 1, this%nrealfluid
      if (this%donor_n(m) > 0) then
        mean_donor_volumes(m) = this%local_wisp_donor_v(m) / this%donor_n(m)
      else
        mean_donor_volumes(m) = 0.0_r8
      end if
    end do

    ! parallel communication
    do m = 1, this%nrealfluid
      block
        real(r8) :: tmp

        this%local_wisp_donor_v(m) = sum(this%wisp_donor_volumes(m,:this%mesh%ncell_onP))
        tmp = sum(this%wisp_acceptor_volumes(m,:this%mesh%ncell_onP))

        this%total_wisp_donor_v(m) = global_sum(this%local_wisp_donor_v(m))
        this%total_wisp_acceptor_v(m) = global_sum(tmp)
      end block
    end do

  end subroutine wisp_volumes

  !---
  ! Determine which materials (if any) do not have enough acceptor
  !---
  subroutine check_acceptor(this, enough_acceptor)
    class(wisp_redistribution), intent(inout) :: this
    logical, intent(out) :: enough_acceptor(:)
    !-
    integer :: m

    enough_acceptor = .true.

    do m = 1, this%nrealfluid
      if (this%total_wisp_donor_v(m) > 0.0_r8 .and. &
          this%total_wisp_acceptor_v(m) < this%total_wisp_donor_v(m)) then
        enough_acceptor(m) = .false.
      end if
    end do

    do m = 1, this%nrealfluid
      write(*, '("ENOUGH_ACCEPTOR[",i2,"]: ", L1)') m, enough_acceptor(m)
    end do

  end subroutine check_acceptor

  !---
  !
  !---
  subroutine sort_donors(this, mean_donor_volumes, enough_acceptor, remove_type, remove_vol)

    use sort_utilities, only: quick_sort
#ifdef NO_2008_FINDLOC
    use f08_intrinsics, only: findloc
#endif


    class(wisp_redistribution), intent(inout) :: this
    real(r8), intent(in) :: mean_donor_volumes(:)
    logical, intent(in) :: enough_acceptor(:)
    integer, intent(out) :: remove_type(:)
    real(r8), intent(out) :: remove_vol(:)
    !-
    integer :: m, i, p, tmp(1)
    integer, allocatable :: sorted_ranks(:)
    real(r8), allocatable :: collated_mean(:), donor_volumes(:), sum_volumes(:)

    allocate(collated_mean(nPE))
    allocate(donor_volumes(nPE))
    allocate(sorted_ranks(nPE))
    allocate(sum_volumes(0:nPE))

    do m = 1, this%nrealfluid
      if (enough_acceptor(m)) then
        remove_type(m) = remove_all
        remove_vol(m) = this%local_wisp_donor_v(m)
        cycle
      end if

      ! there is not enough acceptor fluid so first we communicate our means and donor_volumes
      ! After sorting the volumes based on the means
      sorted_ranks = [(i, i=1, nPE)]
      call collate(collated_mean, mean_donor_volumes(m))
      call collate(donor_volumes, this%local_wisp_donor_v(m))

      call broadcast(collated_mean)
      call broadcast(donor_volumes)

      call quick_sort(collated_mean, sorted_ranks, nPE)

      ! position of process in sorted collection
      tmp = findloc(sorted_ranks, this_PE)
      p = tmp(1)
      sum_volumes = 0.0_r8
      do i = 1, p
        sum_volumes(i) = sum_volumes(i-1) + donor_volumes(sorted_ranks(i))
      end do

      ! critera check for `remove_all`
      if (sum_volumes(p) <= this%total_wisp_acceptor_v(m)) then
        remove_type(m) = remove_all
        remove_vol(m) = donor_volumes(this_PE)
        ! criteria for `remove_some`
      elseif (sum_volumes(p-1) <= this%total_wisp_acceptor_v(m)) then
        remove_type(m) = remove_some
        remove_vol(m) = this%total_wisp_acceptor_v(m) - sum_volumes(p-1)
      else
        ! remove_none is the only option left
        remove_type(m) = remove_none
        remove_vol(m) = 0.0_r8
      endif

    end do

  end subroutine sort_donors

  !---
  !
  !---

  subroutine remove_wisp_donors(this, enough_acceptor, remove_type, remove_volume, vof, flux_vol)
    use sort_utilities, only: quick_sort

    class(wisp_redistribution), intent(inout) :: this
    integer, intent(in) :: remove_type(:)
    real(r8), intent(in) :: remove_volume(:)
    real(r8), intent(inout) :: vof(:,:), flux_vol(:,:)
    logical, intent(in) :: enough_acceptor(:)
    !-
    integer :: i, m, ic, void_i, j
    real(r8) :: sum_before, sum_after, remaining, donated, f
    integer, allocatable :: sorted_ic(:)
    real(r8), allocatable :: wisp_v(:)

    allocate(sorted_ic(this%mesh%ncell_onP))
    allocate(wisp_v(this%mesh%ncell_onP))
    void_i = this%nrealfluid + 1


    do m = 1, this%nrealfluid
      this%uncorrected_n(m) = this%donor_n(m)

      if (this%donor_n(m) == 0) cycle

      if (remove_volume(m) == 0.0_r8 .or. remove_type(m) == remove_none) then
        write(*,*) "REMOVE_NONE"
        ! zero the wisp donor volumes so the corrective procedures in 
        ! heat transfer/flow are not triggered
        this%wisp_donor_volumes(m,:) = 0.0_r8
      else if (remove_type(m) == remove_all) then
        write(*,*) "REMOVE_ALL"
        this%uncorrected_n(m) = 0
        do i = 1, this%mesh%ncell_onP
          ! ideally we would like to get rid of this check and just do
          ! vof(m, i) = vof(m, i) - donated.  However, this can violate the invariant
          ! that for any material vof = {0, [cutoff, 1]}
          if (this%wisp_donor_volumes(m, i) > 0.0_r8) then
            donated = this%wisp_donor_volumes(m, i) / this%mesh%volume(i)

            vof(m, i) = 0.0_r8
            vof(void_i, i) = min(vof(void_i, i) + donated, 1.0_r8)
            ! remove flux vol for this cell
            !flux_vol(m, this%mesh%xcface(i):this%mesh%xcface(i+1)-1) = 0.0_r8
          endif
        end do
      else ! remove_type(m) == remove_some
        write(*,*) "REMOVE_SOME"
        sorted_ic = [ (i, i=1, this%mesh%ncell_onP) ]
        wisp_v = this%wisp_donor_volumes(m,:)

        call quick_sort(wisp_v, sorted_ic, this%mesh%ncell_onP)

        sum_before = 0.0_r8
        i = 1
        do while (sum_before < remove_volume(m))

          ic = sorted_ic(i)
          sum_after = sum_before + wisp_v(i)

          if (sum_after <= remove_volume(m)) then
            if (wisp_v(i) > 0.0_r8) then
              this%uncorrected_n(m) = this%uncorrected_n(m) - 1
              ! remove whole volume
              donated = wisp_v(i) / this%mesh%volume(ic)
              ! handle vof
              vof(m, ic) = 0.0_r8
              vof(void_i, ic) = min(vof(void_i, ic) + donated, 1.0_r8)
              ! remove flux vol for this cell
              ! flux_vol(m, this%mesh%xcface(ic):this%mesh%xcface(ic+1)-1) = 0.0_r8
              ! do we need to remove the flux_vol for neighbors fluxing into/out of ic?
            end if
          else
            ! remove a portion of the largest wisp 
            remaining = remove_volume(m) - sum_before
            j = this%mesh%ncell_onP
            ic = sorted_ic(j)

            ASSERT(wisp_v(j) > remaining)

            ! fraction of material volume to be removed
            f = min(1.0_r8, remaining / (this%mesh%volume(ic) * vof(m, ic)))

            ! absolute vof to be removed
            donated = remaining / this%mesh%volume(ic)
            ! handle vof
            ! XXX we should really have some notion of the vof_cutoff here since there's a
            ! chance that this math here may result in 0 < vof < cutoff...
            vof(m, ic) = vof(m, ic) - donated
            vof(void_i, ic) = min(vof(void_i, ic) + donated, 1.0_r8)
            ! shrink flux vol for this cell

            !flux_vol(m, this%mesh%xcface(ic):this%mesh%xcface(ic+1)-1) = &
            !    flux_vol(m, this%mesh%xcface(ic):this%mesh%xcface(ic+1)-1) * (1.0_r8 - f)

            ! remove remaining donor volume in largest wisp cell
            this%wisp_donor_volumes(m, ic) = max(this%wisp_donor_volumes(m, ic) - remaining, 0.0_r8)

            ! remove remaining donor volumes in rest of cells
            do j = i, this%mesh%ncell_onP-1
              ic = sorted_ic(j)
              this%wisp_donor_volumes(m, ic) = 0.0_r8
            end do
            exit
          end if
          ! advance
          sum_before = sum_after
          i = i + 1
        end do
      end if
    end do
    
    do m = 1, this%nrealfluid
      ! update global quantities with adjusted wisp_volumes
      if (.not.enough_acceptor(m)) then
          this%local_wisp_donor_v(m) = sum(this%wisp_donor_volumes(m,:this%mesh%ncell_onP))
          this%total_wisp_donor_v(m) = global_sum(this%local_wisp_donor_v(m))
      end if
    end do

    ! write(*,"(/,a)") "POST WISP DISTRIBUTION"
    ! write(*,*) "VOF[:,5048]: ", vof(:,5048)
    ! write(*,*) "WISP DONOR VOLUME[:,5048]: ", this%wisp_donor_volumes(:,5048)/this%mesh%volume(5048)
    ! write(*,*) "WISP ACCEPTOR VOLUME[:,5048]: ", this%wisp_acceptor_volumes(:,5048)


  end subroutine remove_wisp_donors

  !---
  !
  !---

  subroutine fill_wisp_acceptors(this, vof, flux_vol)
    class(wisp_redistribution), intent(inout) :: this
    real(r8), intent(inout) :: vof(:,:), flux_vol(:,:)
    !-
    integer :: i, m, void_i
    real(r8) :: total_f, local_vol_f

    void_i = this%nrealfluid + 1

    do m = 1, this%nrealfluid
      if (this%total_wisp_donor_v(m) == 0.0_r8 .or. this%total_wisp_acceptor_v(m) == 0.0_r8) cycle

      ! total fraction of wisp_acceptor used, clip at one in case of roundoff errors
      total_f = min(&
          1.0_r8, &
          min(this%total_wisp_donor_v(m), this%total_wisp_acceptor_v(m)) / this%total_wisp_acceptor_v(m))

      do i = 1, this%mesh%ncell_onP
        if (this%wisp_acceptor_volumes(m, i) == 0.0_r8) cycle

        ! fractional increase of material volume - needed for flux_vol adjustment
        local_vol_f = min(&
            this%wisp_absorption_fraction,&
            total_f * this%wisp_acceptor_volumes(m, i) / (this%mesh%volume(i) * vof(m, i)))


        vof(void_i, i) = vof(void_i, i) - (vof(m, i) * local_vol_f)
        vof(m, i) = vof(m, i) * (1.0_r8 + local_vol_f)

        flux_vol(m, this%mesh%xcface(i):this%mesh%xcface(i+1)-1) = &
            flux_vol(m, this%mesh%xcface(i):this%mesh%xcface(i+1)-1) * (1.0_r8 + local_vol_f)

      end do

    end do

  end subroutine fill_wisp_acceptors

  !--
  ! Convert wisp_acceptor_volumes into a fractional field with all materials combined.
  ! Facilitates the redistribution of scalar quantities associated with wisps
  !--
  subroutine compute_fractional_field(this)
    class(wisp_redistribution), intent(inout) :: this

    integer :: i
    real(r8) :: acceptor_total

    acceptor_total = sum(this%total_wisp_acceptor_v)

    if (acceptor_total > 0.0_r8) then
      do i = 1, this%mesh%ncell_onP
        this%wisp_acceptor_fractions(i) = &
          sum(this%wisp_acceptor_volumes(:,i)) / acceptor_total
      end do
    else
      this%wisp_acceptor_fractions = 0.0_r8
    end if

  end subroutine compute_fractional_field
 
  !--
  !  Compute:
  !  ND -> Total number of donor cells
  !  NA -> Total number of acceptor cells
  !  NU -> Total number of uncorrected cells
  !  VD -> DonorVolume
  !  VD/VA -> DonorVolume/AcceptorVolume
  !---
  subroutine statistics(this)
    class(wisp_redistribution), intent(inout) :: this
    character(len=100) :: message(this%nrealfluid+1)
    integer :: nd, na, nu, m
    real(r8) :: vr

    message(1) = ""
    do m = 1, this%nrealfluid
      nd = global_sum(this%donor_n(m))
      na = global_sum(this%acceptor_n(m))
      nu = global_sum(this%uncorrected_n(m))
      if (this%total_wisp_acceptor_v(m) > 0.0_r8) then
        vr = this%total_wisp_donor_v(m) / this%total_wisp_acceptor_v(m)
      else
        vr = 0.0_r8
      endif
      write(message(m+1), 1) m, nd, na, nu, this%total_wisp_donor_v(m), vr
    end do
1   format('RD: F(',i2.2,'), ND:NA:NU=',2(i7.7,:,':'),i7.7, ', VD=',es9.3,', VD/VA=',f9.6)

    call TLS_info(message)
  end subroutine statistics

end module wisp_redistribution_type

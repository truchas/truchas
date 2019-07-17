!!
!! SIMPLE_VOLUME_TRACKER_TYPE
!!
!! Implements a simplified volume tracker, which advects
!! material based on volume fractions in the donor cell
!! without doing a geometric reconstruction. Most useful
!! when only one material is present, or when fluids aren't
!! expected to interact, or when fluid diffusion is acceptable.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module simple_volume_tracker_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use volume_tracker_2d_class
  use unstr_2d_mesh_type
  implicit none
  private

  type, extends(volume_tracker_2d), public :: simple_volume_tracker
    private
    type(unstr_2d_mesh), pointer :: mesh ! unowned reference
    logical :: is_axisym
  contains
    procedure :: init
    procedure :: flux_volumes
    procedure :: set_inflow_material
  end type simple_volume_tracker

contains

  subroutine init(this, mesh, nrealfluid, nfluid, nmat, axisym)

    class(simple_volume_tracker), intent(out) :: this
    type(unstr_2d_mesh), intent(in), target :: mesh
    integer, intent(in) :: nrealfluid, nfluid, nmat
    logical, intent(in) :: axisym

    this%mesh => mesh
    this%is_axisym = axisym

  end subroutine init

  subroutine flux_volumes(this, vel, vof_n, vof, flux_vol, int_normal, fluids, void, dt)

    class(simple_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vof_n(:,:), dt
    real(r8), intent(out) :: flux_vol(:,:), vof(:,:), int_normal(:,:,:)
    integer, intent(in) :: fluids, void

    integer :: i, j, k, n, f0, f1

    ! compute upwind flux volumes for transport
    do i = 1, this%mesh%ncell_onP
      f0 = this%mesh%cstart(i)
      f1 = this%mesh%cstart(i+1)-1

      do j = f0, f1
        k = this%mesh%cface(j)
        n = this%mesh%cnhbr(j)
        if (vel(j) > 0 .or. n == 0) then
          ! if the donator cell is off-process and not a ghost cell, this flux is irrelevant.
          if (this%mesh%fcell(1,k) > this%mesh%ncell .or. this%mesh%fcell(1,k) == 0) cycle
          flux_vol(:,j) = vel(j)*this%mesh%area(k)*dt * vof_n(:,i)
        else
          flux_vol(:,j) = vel(j)*this%mesh%area(k)*dt * vof_n(:,n)
        end if
      end do

      vof(:,i) = vof_n(:,i)
      vof(:,i) = vof(:,i) - sum(flux_vol(:,f0:f1), dim=2) / this%mesh%volume(i)
    end do

    int_normal = 0.0_r8

    call this%mesh%cell_imap%gather_offp(vof)

  end subroutine flux_volumes

  !! Set the inflow material for the given boundary faces. A material index
  !! of 0 will result in materials being fluxed in proportion to the material
  !! fractions present in the cell. This is the preset default. The BC_INDEX
  !! array component is ordered, so we use a binary search to locate the faces
  !! in the the array, so we can set the corresponding inflow material index.
  !! TODO: If FACES is also ordered (likely) the search can be improved further.

  subroutine set_inflow_material(this, mat, faces)
    use truchas_logging_services
    class(simple_volume_tracker), intent(inout) :: this
    integer, intent(in) :: mat  ! material index
    integer, intent(in) :: faces(:) ! face indices
    call TLS_warn('****************************************************')
    call TLS_warn('INFLOW BC NOT IMPLEMENTED FOR SIMPLE VOLUME TRACKING')
    call TLS_warn('****************************************************')
    ! Note that some of the code and data for this would be identical to
    ! what is in geometric_volume_tracker and could probably be pushed
    ! up into the base class.
  end subroutine set_inflow_material

end module simple_volume_tracker_type

!!
!! SIMPLE_VOLUME_TRACKER_TYPE
!!
!! Implements a simplified volume tracker, which advects
!! material based on volume fractions in the donor cell
!! without doing a geometric reconstruction. Most useful
!! when only one material is present, or when fluids aren't
!! expected to interact, or when fluid diffusion is acceptable.

module simple_volume_tracker_type

  use kinds, only: r8
  use volume_tracker_class
  use unstr_mesh_type
  implicit none
  private

  type, extends(volume_tracker), public :: simple_volume_tracker
    private
    type(unstr_mesh), pointer :: mesh ! unowned reference
  contains
    procedure :: init
    procedure :: flux_volumes
  end type simple_volume_tracker

contains

  subroutine init(this, mesh, nrealfluid, nfluid, nmat, liq_matid, params)

    use parameter_list_type

    class(simple_volume_tracker), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    integer, intent(in) :: nrealfluid, nfluid, nmat, liq_matid(:)
    type(parameter_list), intent(inout) :: params

    this%mesh => mesh

  end subroutine init

  subroutine flux_volumes(this, vel, vof_n, vof, flux_vol, fluids, void, dt)

    use index_partitioning, only: gather_boundary

    class(simple_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vof_n(:,:), dt
    real(r8), intent(out) :: flux_vol(:,:), vof(:,:)
    integer, intent(in) :: fluids, void

    integer :: i, j, k, f0, f1

    associate (m => this%mesh)
      ! compute upwind flux volumes for transport
      do i = 1, m%ncell
        f0 = m%xcface(i)
        f1 = m%xcface(i+1)-1
        do j = f0, f1
          k = m%cface(j)
          if (vel(j) > 0 .or. m%fcell(2,k) == 0) then
            ! if the donator cell is off-process and not a ghost cell, this flux is irrelevant.
            if (m%fcell(1,k) > m%ncell .or. m%fcell(1,k) == 0) cycle
            flux_vol(:,j) = vel(j)*m%area(k)*dt * vof_n(:fluids+void,m%fcell(1,k))
          else
            if (m%fcell(2,k) > m%ncell) cycle
            flux_vol(:,j) = vel(j)*m%area(k)*dt * vof_n(:fluids+void,m%fcell(2,k))
          end if
        end do
      end do

      ! update volume fractions
      do i = 1, m%ncell_onP
        f0 = m%xcface(i)
        f1 = m%xcface(i+1)-1
        vof(:,i) = vof_n(:,i)
        vof(:fluids+void,i) = vof(:fluids+void,i) + sum(flux_vol(:,f0:f1), dim=2)
      end do
    end associate

    call gather_boundary(this%mesh%cell_ip, vof)

  end subroutine flux_volumes

end module simple_volume_tracker_type

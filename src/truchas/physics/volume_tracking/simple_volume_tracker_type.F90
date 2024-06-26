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
  use volume_tracker_class
  use unstr_mesh_type
  use vector_func_containers
  implicit none
  private

  type, extends(volume_tracker), public :: simple_volume_tracker
    private
    type(unstr_mesh), pointer :: mesh ! unowned reference
  contains
    procedure :: init
    procedure :: flux_volumes
    procedure :: set_inflow_material
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

  subroutine flux_volumes(this, vel, vof_n, vof, flux_vol, fluids, void, dt, t)

    class(simple_volume_tracker), intent(inout) :: this
    real(r8), intent(in) :: vel(:), vof_n(:,:), dt, t
    real(r8), intent(out) :: flux_vol(:,:), vof(:,:)
    integer, intent(in) :: fluids, void

    integer :: i, j, k, n, f0, f1

    associate (m => this%mesh)
      do i = 1, m%ncell_onP
        f0 = m%xcface(i)
        f1 = m%xcface(i+1)-1

        do j = f0, f1
          k = m%cface(j)
          n = m%cnhbr(m%xcnhbr(i)+j-f0)
          ! If the neighbor cell doesn't exist, I'm a boundary cell and provide
          ! the "upwind" data: more of whatever I have is the inflow volume
          ! fractions.
          if (vel(j) > 0 .or. n == 0) then
            flux_vol(:,j) = vel(j)*m%area(k)*dt * vof_n(:fluids+void,i)
          else
            flux_vol(:,j) = vel(j)*m%area(k)*dt * vof_n(:fluids+void,n)
          end if
        end do

        vof(:,i) = vof_n(:,i)
        vof(:fluids+void,i) = vof(:fluids+void,i) - sum(flux_vol(:,f0:f1), dim=2) / m%volume(i)
      end do
    end associate

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
    class(vector_func), intent(in) :: mat ! material inflow vector function
    integer, intent(in) :: faces(:) ! face indices
    call TLS_warn('****************************************************')
    call TLS_warn('INFLOW BC NOT IMPLEMENTED FOR SIMPLE VOLUME TRACKING')
    call TLS_warn('****************************************************')
    ! Note that some of the code and data for this would be identical to
    ! what is in geometric_volume_tracker and could probably be pushed
    ! up into the base class.
  end subroutine set_inflow_material

end module simple_volume_tracker_type

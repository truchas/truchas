#include "f90_assert.fpp"

module enthalpy_redistributor_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use avg_phase_prop_type
  use parameter_list_type
  implicit none
  private

  type, public :: enthalpy_redistributor
    type(unstr_mesh), pointer :: mesh => null()
    real(r8), pointer :: wisp_donor_volumes(:,:) => null()
    real(r8), pointer :: wisp_acceptor_fractions(:) => null()
    real(r8), pointer :: vtrack_vofs(:,:) => null() ! for debugging only - should be removed
    real(r8), pointer :: vtrack_vofs_old(:,:) => null() ! for debugging only - should be removed
    type(avg_phase_prop), allocatable :: enthalpy
    logical :: do_redistribute
  contains
    procedure :: init
    procedure :: get_redistributed_enthalpy
  end type

contains

  subroutine init(this, mesh)
    use vtrack_driver, only: vtrack_liq_matid_view, vtrack_wisp_donors_view, &
                             vtrack_wisp_acceptor_fractions_view, vtrack_vof_view, vtrack_vof_old_view
    use material_model_driver, only: matl_model
    class(enthalpy_redistributor), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    
    character(:), allocatable :: errmsg
    integer, pointer :: matid(:)

    this%mesh => mesh
    this%wisp_donor_volumes => vtrack_wisp_donors_view()
    this%wisp_acceptor_fractions => vtrack_wisp_acceptor_fractions_view()
    this%vtrack_vofs => vtrack_vof_view()
    this%vtrack_vofs_old => vtrack_vof_old_view()
    this%do_redistribute = &
        associated(this%wisp_acceptor_fractions) .and. associated(this%wisp_donor_volumes)    

    matid => vtrack_liq_matid_view()
    call matl_model%alloc_avg_phase_prop('enthalpy', matid, this%enthalpy, errmsg)
   
  end subroutine init

  subroutine get_redistributed_enthalpy(this, tcell, dq, hlast, h, tot_void_cell)

    use parallel_communication

    class(enthalpy_redistributor), intent(in) :: this
    real(r8), intent(in) :: tcell(:)
    real(r8), intent(in) :: h(:), hlast(:)
    logical, intent(in) :: tot_void_cell(:)
    real(r8), intent(out) :: dq(:)

    integer :: j
    real(r8) :: dq_donor, dq_i, state(1), tmp

    if (.not.this%do_redistribute) then
        dq = 0.0_r8
        return
    end if

    ASSERT(size(dq) == this%mesh%ncell_onP)    
    ASSERT(size(tcell) == this%mesh%ncell_onP)
    INSIST(size(this%vtrack_vofs, DIM=1) == 3)

    ! Compute the loss of enthalpy due to wisp removal from donor volumes
    dq_donor = 0.0_r8
    do j = 1, this%mesh%ncell_onP
        if (.not.tot_void_cell(j) .and. any(this%wisp_donor_volumes(:,j) > 0)) then
            state(1) = tcell(j)            
            call this%enthalpy%compute_value(this%wisp_donor_volumes(:,j), state, dq_i)
            tmp = dq_i / this%mesh%volume(j)
            ! check for strange values
            if (tmp < 0.0_r8) then
              write(*, '("[", i6, "] *** NEGATIVE Q *** TVOID/STATE/HLAST/H/DQ: ", l1, 3es15.5, /, "STATE/DV", 2es15.5)') &
              j, tot_void_cell(j), hlast(j), h(j), tmp, state(1), this%wisp_donor_volumes(1,j)
            else if (tmp > h(j)) then
              write(*, '("[", i6, "] *** DQ > H *** HLAST/H/DQ: ", 3es15.5)') j, hlast(j), h(j), tmp
            end if
            dq(j) = -dq_i
            dq_donor = dq_donor + dq_i
        else
            dq(j) = 0.0_r8
        end if
    end do

    dq_donor = global_sum(dq_donor)
    write(*,*) '>>> WISP ENTHALPY: ', dq_donor
    write(*,*) '>>> ACCEPTOR ENTHALPY: ', dq_donor * sum(this%wisp_acceptor_fractions(:this%mesh%ncell_onP))
    ! Compute enthalpy gains due to wisp additions to acceptor cells
    do j = 1, this%mesh%ncell_onP
        dq(j) = dq(j) + this%wisp_acceptor_fractions(j)*dq_donor
    end do 

  end subroutine get_redistributed_enthalpy

end module enthalpy_redistributor_type

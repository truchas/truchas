module evap_heat_flux_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use tdep_bndry_func_class
  implicit none
  private

  type, extends(tdep_bndry_func), public :: evap_heat_flux
    real(r8), private :: a, beta, b
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_value => compute
    procedure :: compute_deriv => compute
  end type

contains

  subroutine compute(this, t, var)
    class(evap_heat_flux), intent(inout) :: this
    real(r8), intent(in) :: t, var(:)
    integer  :: j
    real(r8) :: temp
    if (this%beta == 0.0_r8) then
      do j = 1, size(this%index)
        temp = var(this%index(j))
        this%value(j) = this%a * exp(-this%b/temp)
        this%deriv(j) = this%value(j) * (this%b/temp**2)
      end do
    else
      do j = 1, size(this%index)
        temp = var(this%index(j))
        this%value(j) = this%a * (temp**this%beta) * exp(-this%b/temp)
        this%deriv(j) = this%value(j) * (this%b/temp + this%beta)/temp
      end do
    end if
  end subroutine compute

  subroutine init(this, mesh, params, stat, errmsg)

    use bitfield_type
    use unstr_mesh_type
    use parameter_list_type
    use parallel_communication, only: global_any

    class(evap_heat_flux), intent(out) :: this
    type(unstr_mesh), intent(in) :: mesh
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, n
    logical, allocatable :: mask(:)
    integer, allocatable :: setids(:)
    type(bitfield) :: bitmask
    real(r8), parameter :: R = 8.3144598_r8 ! Gas constant [J/mol-K]

    call params%get('prefactor', this%a)
    call params%get('temp-exponent', this%beta, default=0.0_r8)
    call params%get('activation-energy', this%b)
    this%b = this%b/R

    call params%get('face-set-ids', setids)
    call mesh%get_face_set_bitmask(setids, bitmask, stat, errmsg)
    if (stat /= 0) return

    !! Identify the faces specified by SETIDS
    mask = (popcnt(iand(bitmask, mesh%face_set_mask)) /= 0)

    !! Verify these are all boundary faces.
    if (global_any(mask .and. .not.btest(mesh%face_set_mask,0))) then
      stat = 1
      errmsg = 'some faces not boundary faces'
      return
    end if

    n = count(mask(:mesh%nface_onP))
    allocate(this%index(n), this%value(n), this%deriv(n))

    n = 0
    do j = 1, mesh%nface_onP
      if (mask(j)) then
        n = n + 1
        this%index(n) = j
      end if
    end do

  end subroutine init

end module evap_heat_flux_type

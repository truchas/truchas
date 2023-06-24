module evap_heat_flux_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_func2_class
  implicit none
  private

  type, extends(bndry_func2), public :: evap_heat_flux
    real(r8), private :: a, b, temp_0
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
    do j = 1, size(this%index)
      temp = var(this%index(j))
      if (temp > 0.0_r8) then
        this%value(j) = (this%a / sqrt(temp)) * &
                        exp(this%b*(1.0_r8/this%temp_0 - 1.0_r8/temp))
        this%deriv(j) = this%value(j) * (this%b/temp - 0.5_r8) / temp
      else ! wonky temperature data -- avoid invalid floating exception
        this%value(j) = 0.0_r8
        this%deriv(j) = 0.0_r8
      end if
    end do
  end subroutine compute

  subroutine init(this, mesh, params, stat, errmsg)

    use bitfield_type
    use unstr_mesh_type
    use parameter_list_type
    use parallel_communication, only: global_any
    use physical_constants, only: R => gas_constant

    class(evap_heat_flux), intent(out) :: this
    type(unstr_mesh), intent(in) :: mesh
    type(parameter_list) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: j, n
    real(r8) :: L, M, p0, lambda
    logical, allocatable :: mask(:)
    integer, allocatable :: setids(:)
    type(bitfield) :: bitmask

    real(r8), parameter :: TWOPI = 6.2831853071795862_r8

    call params%get('vaporization-heat', L, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (L < 0) then
      stat = -1
      errmsg = '"vaporization-heat" is < 0'
      return
    end if

    call params%get('vaporization-temp', this%temp_0, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (this%temp_0 <= 0) then
      stat = -1
      errmsg = '"vaporization-temp" is <= 0'
      return
    end if

    call params%get('molar-mass', M, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (M <= 0) then
      stat = -1
      errmsg = '"molar-mass" is <= 0'
      return
    end if

    ! default is 1 atm in SI units
    call params%get('ambient-pressure', p0, default=1.01325e5_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (p0 <= 0) then
      stat = -1
      errmsg = '"ambient-pressure" is <= 0'
      return
    end if

    call params%get('condensation-factor', lambda, default=0.1_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) return
    if (lambda < 0) then
      stat = -1
      errmsg = '"condensation-factor" is < 0'
      return
    end if

    !! Form the coefficients in the mathematical form of the flux
    this%a = L * lambda * p0 * sqrt(M/(TWOPI*R))
    this%b = M * L / R

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

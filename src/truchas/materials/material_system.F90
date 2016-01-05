!!
!! MATERIAL_SYSTEM
!!
!! Provides a derived type for describing a material system and procedures
!! that operate on instances of that type.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! The derived type MAT_SYSTEM is an opaque structure that encapsulates the
!! data associated with a material system.  This includes material components
!! and phases, phase diagram, and state variables.
!!
!! The module provides a defined assignment operator for instances of the
!! MAT_SYSTEM type that does a deep copy, resulting in a lhs that is
!! completely independent of the rhs.
!!
!!  MS_NUM_PHASE(THIS) returns the number of phases that comprise the material
!!    system THIS.
!!
!!  MS_NUM_COMPONENT(THIS) returns the number of components in the material
!!    system THIS.
!!
!!  MS_TEMP_DEP(THIS) returns true if the material system THIS includes
!!    temperature as one of its state variables; otherwise it returns false.
!!
!!  MS_TEMP_LO(THIS, N) returns the low temperature of the the Nth phase
!!    transformation interval for the material system THIS.  It is an error
!!    if N is not a valid transformation index.
!!
!!  MS_TEMP_HI(THIS, N) returns the high temperature of the the Nth phase
!!    transformation interval for the material system THIS.  It is an error
!!    if N is not a valid transformation index.
!!
!!  MS_LATENT_HEAT(THIS, N) returns the latent_heat of the the Nth phase
!!    transformation for the material system THIS.  It is an error if N is
!!    not a valid transformation index.
!!
!!  CALL MS_GET_PHASE_ID (THIS, ID) returns IDs of the phases that form the
!!    material system THIS in the integer pointer array ID.  This array is
!!    allocated by the subroutine to the proper size, and the caller is
!!    resposible for deallocating it.
!!
!!  CALL MS_PHASE_MIXTURE (THIS, STATE, BETA) returns the material phase
!!    fractions BETA for the material system THIS given the value of the
!!    state variables passed in the rank-1 real array STATE.  If ID is the
!!    array of phase IDs returned by MS_GET_PHASE_ID, then BETA and ID have
!!    the same shape and BETA(j) is the fraction of phase ID(j) present at
!!    the specified state.  We have BETA >= 0 and SUM(BETA) = 1.  This is a
!!    query of the material system's phase diagram, and as such should only
!!    be called with multi-phase material systems (necessarily temperature
!!    dependent).  The first element of STATE is assumed to be temperature,
!!    and the remaining elements the component concentrations.  Note that
!!    only single-component material systems (unary phase diagrams) are
!!    currently implemented.
!!
!!  CALL MS_CREATE (THIS, NUM_COMP, TEMP_DEP, PHASE_ID, TEMP_LO, TEMP_HI,
!!    LATENT_HEAT, EPS, REF_TEMP, REF_ENTHALPY) creates a material system
!!    THIS having the attributes specified by the remaining arguments:
!!
!!    NUM_COMP -- the number of components in the material system
!!    TEMP_DEP -- true if the material system has temperature as a state var
!!    PHASE_ID -- array of the IDs of the phases that form the material system.
!!      For single-component material systems, which have a unary phase diagram,
!!      the order of the IDs is significant and are assumed to be in
!!    TEMP_LO/TEMP_HI -- specifies the (smeared) unarray phase diagram for
!!      single-component, multi-phase material systems.  The size of the
!!      arrays are 1 less than the number of phases.  (TEMP_LO(j),TEMP_HI(j))
!!      is the temperature interval of the jth phase transformation (smeared
!!      out from the proper isothermal transformation of a unary phase diagram)
!!      from PHASE_ID(j) to PHASE_ID(j+1), and the transformations are required
!!      to be non-overlapping and ordered with increasing temperature.
!!    EPS -- specifies the smoothing radius at the lo/hi corners of the phase
!!      transformation intervals.  Within those intervals linear interpolation
!!      nominally defines the phase fractions, which results in non-smooth
!!      behavior at the end-points, but this is smoothed using a polynomial
!!      interpolant in a EPS(j) neighborhood of the end-points of the jth
!!      transformation interval.  EPS must be non-negative and no larger than
!!      half the transformation width.
!!    LATENT_HEAT -- array of latent heats associated with each of the phase
!!      transformation intervals.
!!    REF_TEMP/REF_ENTHALPY -- specifies the enthalpy of the material system
!!      at a reference temperature.  Used to uniquely define the enthalpy
!!      function from the user-specified heat capacity function.  The reference
!!      temperature must lie in the domain of heat capacity function of the
!!      first (or only) phase.
!!
!!  CALL DESTROY (THIS) deallocates all storage associated with the material
!!      system THIS, returning THIS to its default initialization state.
!!

#include "f90_assert.fpp"

module material_system

  use kinds
  implicit none
  private

  public :: ms_create, destroy
  public :: ms_num_phase, ms_num_component, ms_temp_dep, ms_get_phase_id
  public :: ms_temp_lo, ms_temp_hi, ms_latent_heat, ms_ref_temp, ms_ref_enthalpy
  public :: ms_phase_mixture

  type, public :: mat_system
    private
    integer :: num_phase = 0
    integer :: num_component = 0
    logical :: temp_dep = .true.
    integer, pointer :: phase_id(:) => null()
    !! Unary phase diagram data (hardwired structures)
    real(r8), pointer :: temp_hi(:) => null()
    real(r8), pointer :: temp_lo(:) => null()
    real(r8), pointer :: eps(:) => null()  ! smoothing radii
    !! Used to define the enthalpy function only (don't really belong here)
    real(r8), pointer :: latent_heat(:) => null()
    real(r8) :: ref_temp = 0.0_r8, ref_enthalpy = 0.0_r8
  end type mat_system

  public :: assignment(=)
  interface assignment(=)
    module procedure copy_mat_system
  end interface

  interface destroy
    module procedure destroy_mat_system
  end interface

contains

  integer function ms_num_phase (this)
    type(mat_system), intent(in) :: this
    ms_num_phase = this%num_phase
  end function ms_num_phase

  integer function ms_num_component (this)
    type(mat_system), intent(in) :: this
    ms_num_component = this%num_component
  end function ms_num_component

  logical function ms_temp_dep (this)
    type(mat_system), intent(in) :: this
    ms_temp_dep = this%temp_dep
  end function ms_temp_dep

  subroutine ms_get_phase_id (this, phase_id)
    type(mat_system), intent(in) :: this
    integer, pointer :: phase_id(:)
    ASSERT(associated(this%phase_id))
    allocate(phase_id(size(this%phase_id)))
    phase_id = this%phase_id
  end subroutine ms_get_phase_id
  
  function ms_temp_lo (this, n) result (t)
    type(mat_system), intent(in) :: this
    integer, intent(in) :: n
    real(r8) :: t
    INSIST(n > 0 .and. n <= size(this%temp_lo))
    t = this%temp_lo(n)
  end function ms_temp_lo

  function ms_temp_hi (this, n) result (t)
    type(mat_system), intent(in) :: this
    integer, intent(in) :: n
    real(r8) :: t
    INSIST(n > 0 .and. n <= size(this%temp_hi))
    t = this%temp_hi(n)
  end function ms_temp_hi
  
  function ms_latent_heat (this, n) result (l)
    type(mat_system), intent(in) :: this
    integer, intent(in) :: n
    real(r8) :: l
    INSIST(n > 0 .and. n <= size(this%latent_heat))
    l = this%latent_heat(n)
  end function ms_latent_heat

  function ms_ref_temp (this) result (t)
    type(mat_system), intent(in) :: this
    real(r8) :: t
    t = this%ref_temp
  end function
  
  function ms_ref_enthalpy (this) result (h)
    type(mat_system), intent(in) :: this
    real(r8) :: h
    h = this%ref_enthalpy
  end function

  subroutine ms_phase_mixture (this, state, beta)

    type(mat_system) :: this
    real(r8), intent(in)  :: state(:)
    real(r8), intent(out) :: beta(:)

    integer :: p
    logical :: two_phase
    real(r8) :: alpha, temp

    ASSERT(size(beta) == this%num_phase)

    if (this%num_phase == 1) then ! we shouldn't have been called, but ...

      beta(1) = 1.0_r8

    else

      ASSERT(this%temp_dep)
      ASSERT(size(state) >= 1)
      temp = state(1)

      !! Find p such that TEMP_HI(p-1) + EPS(p-1) <= TEMP < TEMP_HI(p) + EPS(p);
      !! this implies phase p is present in this region.
      p = 1
      do while (temp >= this%temp_hi(p) + this%eps(p))
        p = p + 1
        if (p > this%num_phase-1) exit
      end do

      !! Determine whether this is a single-phase region with phase p only,
      !! or a two-phase region with phases p and p+1.
      two_phase = .false.
      if (p <= this%num_phase-1) then
        if (temp > this%temp_lo(p) - this%eps(p)) two_phase = .true.
      end if

      if (two_phase) then

        !! Compute the fraction ALPHA of phase p+1.  Normally this would just
        !! be the linear interpolant from 0 to 1 over the lo to hi temperature
        !! interval.  Instead we smooth in an EPS-neighborhood of the corners
        !! at the lo and hi temperatures using a polynomial with matched values,
        !! first and second derivatives to the piecewise-linear ALPHA.

#define C2_SMOOTHING
#ifdef C2_SMOOTHING
        alpha = C2_ramp(this%temp_lo(p), this%temp_hi(p), this%eps(p), temp)
#else
        alpha = C1_ramp(this%temp_lo(p), this%temp_hi(p), this%eps(p), temp)
#endif

        ASSERT( 0.0_r8 <= alpha .and. alpha <= 1.0_r8 )

        beta = 0.0_r8
        beta(p) = 1.0_r8 - alpha
        beta(p+1) = alpha

      else ! single-phase region with phase p

        beta = 0.0_r8
        beta(p) = 1.0_r8

      end if

    end if

  contains

    function C1_ramp (x1, x2, eps, x) result (alpha)

      real(r8), intent(in) :: x1, x2, eps, x
      real(r8) :: alpha, z

      if (x < x1 + eps) then
        z = (x - x1) / eps
        alpha = (eps/(x2-x1))*(1.0d0+z)**2/4.0d0
      else if (x > x2 - eps) then
        z = (x - x2) / eps
        alpha = 1.0_r8 - (eps/(x2-x1))*(1.0d0-z)**2/4.0d0
      else ! the usual linear form
        alpha = (x - x1) / (x2 - x1)
      end if

    end function C1_ramp

    function C2_ramp (x1, x2, eps, x) result (alpha)

      real(r8), intent(in) :: x1, x2, eps, x
      real(r8) :: alpha, z

      if (x < x1 + eps) then
        z = (x - x1) / eps
        alpha = (eps/(x2-x1))*(1.0d0+z)**3*(3.0d0-z)/16.0d0
      else if (x > x2 - eps) then
        z = (x - x2) / eps
        alpha = 1.0_r8 - (eps/(x2-x1))*(1.0d0-z)**3*(3.0d0+z)/16.0d0
      else ! the usual linear form
        alpha = (x - x1) / (x2 - x1)
      end if

    end function C2_ramp

  end subroutine ms_phase_mixture

  subroutine ms_create (this, num_component, temp_dep, phase_id, &
        temp_lo, temp_hi, latent_heat, eps, ref_temp, ref_enthalpy)

    use phase_property_table, only: ppt_valid_phase

    type(mat_system), intent(out) :: this
    integer,  intent(in) :: num_component
    logical,  intent(in) :: temp_dep
    integer,  intent(in) :: phase_id(:)
    real(r8), intent(in) :: temp_lo(:), temp_hi(:), latent_heat(:), eps(:)
    real(r8), intent(in) :: ref_temp, ref_enthalpy

    integer :: n

    INSIST( num_component > 0 )
    INSIST( size(phase_id) > 0 )

    this%num_component = num_component
    this%num_phase = size(phase_id)

    INSIST(temp_dep .or. (num_component>1 .and. this%num_phase==1))
    this%temp_dep = temp_dep

    INSIST( distinct(phase_id) )
    do n = 1, size(phase_id)
      INSIST( ppt_valid_phase(phase_id(n)) )
    end do

    allocate(this%phase_id(this%num_phase))
    this%phase_id = phase_id

    if (this%num_phase > 1) then  ! have a phase diagram

      if (this%num_component == 1) then ! unary phase diagram

        n = this%num_phase - 1

        INSIST( size(temp_lo) == n )
        INSIST( valid_unary_diagram(temp_lo, temp_hi, eps) )
        INSIST( size(latent_heat) == n )
        INSIST( all(latent_heat >= 0.0_r8) )

        allocate(this%temp_lo(n), this%temp_hi(n), this%latent_heat(n), this%eps(n))
        this%temp_lo = temp_lo
        this%temp_hi = temp_hi
        this%latent_heat = latent_heat
        this%eps = eps

      else

        !! Multi-component phase diagrams are not yet implemented.
        INSIST( .false. )

      end if

    else  ! single-phase -- no phase diagram

      allocate(this%temp_lo(0), this%temp_hi(0), this%eps(0), this%latent_heat(0))

    end if

    this%ref_temp = ref_temp
    this%ref_enthalpy = ref_enthalpy

  contains

    logical function distinct (array)
      integer, intent(in) :: array(:)
      integer :: i, j
      distinct = .false.
      do i = 1, size(array)-1
        do j = i+1, size(array)
          if (array(i) == array(j)) return
        end do
      end do
      distinct = .true.
    end function

    logical function valid_unary_diagram (temp_lo, temp_hi, eps)

      real(r8), intent(in) :: temp_lo(:), temp_hi(:), eps(:)

      integer :: n, j

      valid_unary_diagram = .false.

      n = size(temp_lo)
      if (size(temp_hi) /= n .or. size(eps) /= n) return

      !! Check individual transformation intervals for correctness.
      if (any(temp_hi <= temp_lo)) return

      !! Check smoothing radii for correctness.
      if (any(eps < 0.0_r8)) return
      if (any(eps > (temp_hi-temp_lo)/2)) return

      !! Check for properly ordered transformation intervals
      do j = 1, n-1
        if (temp_hi(j) >= temp_lo(j+1)) return
      end do

      !! Check that smoothed transformation intervals remain separated.
      do j = 1, n-1
        if (temp_hi(j) + eps(j) > temp_lo(j+1) - eps(j+1)) return
      end do

      valid_unary_diagram = .true.

    end function

  end subroutine ms_create

  subroutine destroy_mat_system (this)
    type(mat_system), intent(inout) :: this
    type(mat_system) :: default
    if (associated(this%phase_id)) deallocate(this%phase_id)
    if (associated(this%temp_lo)) deallocate(this%temp_lo)
    if (associated(this%temp_hi)) deallocate(this%temp_hi)
    if (associated(this%latent_heat)) deallocate(this%latent_heat)
    if (associated(this%eps)) deallocate(this%eps)
    this = default  ! assign default initialization values
  end subroutine destroy_mat_system

  subroutine copy_mat_system (dest, src)
    type(mat_system), intent(out) :: dest
    type(mat_system), intent(in)  :: src
    dest%num_phase = src%num_phase
    dest%num_component = src%num_component
    dest%temp_dep = src%temp_dep
    !! fix up the pointer components
    if (associated(src%phase_id)) then
      allocate(dest%phase_id(size(src%phase_id)))
      dest%phase_id = src%phase_id
    end if
    if (associated(src%temp_lo)) then
      allocate(dest%temp_lo(size(src%temp_lo)))
      dest%temp_lo = src%temp_lo
    end if
    if (associated(src%temp_hi)) then
      allocate(dest%temp_hi(size(src%temp_hi)))
      dest%temp_hi = src%temp_hi
    end if
    if (associated(src%latent_heat)) then
      allocate(dest%latent_heat(size(src%latent_heat)))
      dest%latent_heat = src%latent_heat
    end if
    if (associated(src%eps)) then
      allocate(dest%eps(size(src%eps)))
      dest%eps = src%eps
    end if
    dest%ref_temp = src%ref_temp
    dest%ref_enthalpy = src%ref_enthalpy
  end subroutine copy_mat_system

end module material_system

!!
!! USTRUC_LDRD
!!
!! A concrete implementation of USTRUC_COMP/USTRUC_ANALYSIS that adds the
!! identification of microstructure type and characteristics.
!!
!! This implements an initial model from Seth Imhoff (MST-6) that takes as
!! input the thermal gradient magnitude (G) and solidification front speed (V)
!! at a certain point during soldification, and does a lookup for the type of
!! microstructure (planar, cellular, dendritic) and for the dendritic case,
!! it computes the primary and secondary arm spacing using explicit formulas.
!! This work was funded by LDRD 20140639ER "Solute and Microstructure
!! Prediction during Processing".
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! September 2014; updated July 2015
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  See the comments that accompany the source for the USTRUC_ANALYSIS and
!!  USTRUC_COMP classes which define the interface to this object.  The
!!  only public entity provided by the module is the following function.
!!
!!  NEW_USTRUC_LDRD(COMP, PARAMS) returns a pointer to a new USTRUC_ANALYSIS
!!    class object whose dynamic type is USTRUC_LDRD.  The USTRUC_LDRD type is
!!    itself private. COMP is a USTRUC_ANALYSIS pointer whose target is being
!!    wrapped by this new analysis component specified by PARAMS. The new
!!    object takes ownership of the target. PARAMS is a PARAMETER_LIST object;
!!    the relevant parameters in PARAMS are:
!!
!!    theta1 -- low solid fraction threshold; solidification is deemed
!!        to have started when the solid fraction crosses this threshold.
!!    theta2 -- high solid fraction threshold; solidification is deemed
!!       to have finished when the solid fraction crosses this threshold.
!!    theta1p -- if the solid fraction drops below this threshold while
!!       solidifying, it is deemed to have returned to a liquid state;
!!       theta1p <= theta1, optional, default theta1p = theta1.
!!    theta2p -- if the solid fraction drops below this threshold while
!!       considered solid, it is considered to have remelted and the
!!       previously computed solidification time erased.
!!    theta-gv -- the solid fraction at which G and V are taken
!!
!!    The remaining parameters are specific to the model (see other
!!    documentation for their meaning): liquidus-slope, solute-conc,
!!    partition-coef, liq-sol-delta-T, diffusivity, gibbs-thomson-coef,
!!    instability-coef, and coarsening-coef.
!!
!!  Objects of this type respond to the following data names in the generic
!!  GET subroutine: 'ldrd-invalid-GV', 'ldrd-t_sol', 'ldrd-type',
!!  'ldrd-lambda1', 'ldrd-lambda2', 'ldrd-G', and 'ldrd-V'.
!!

#include "f90_assert.fpp"

module ustruc_ldrd_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use ustruc_comp_class
  implicit none
  private

  public :: new_ustruc_ldrd

  !! The microstructure selection model
  type :: gv
    real(r8) :: theta     ! solid fraction threshold (for picking off G and V)
    real(r8) :: m_liq     ! liquidus slope (should be negative)
    real(r8) :: c_0       ! initial concentration
    real(r8) :: k_part    ! partition coefficient (should be > 1)
    real(r8) :: delta_T   ! liquidus-solidus temperature difference
    real(r8) :: d         ! diffusion coefficient
    real(r8) :: gamma     ! Gibbs-Thomson coefficient
    real(r8) :: alpha     ! instability proportionality constant
    real(r8) :: k_coarse  ! coarsening constant
  contains
    procedure :: init   => gv_init
    procedure :: start  => gv_start
    procedure :: update => gv_update
    procedure :: finish => gv_finish
    procedure :: reset  => gv_reset
    procedure :: lambda1 => gv_lambda1
    procedure :: lambdaCD => gv_lambdaCD
  end type gv

  !! GV_STATE ustruc values
  integer, parameter :: GV_INVALID   = -1
  integer, parameter :: GV_UNDEFINED = 0
  integer, parameter :: GV_PLANAR    = 1
  integer, parameter :: GV_CELLULAR  = 2
  integer, parameter :: GV_DENDRITIC = 3

  !! Microstructure state data
  type :: gv_state
    real(r8) :: G = 0.0_r8        ! thermal gradient magnitude at onset
    real(r8) :: V = 0.0_r8        ! solidification front speed at onset
    real(r8) :: sfrac = 0.0_r8    ! current solid fraction
    real(r8) :: lambda1 = 0.0_r8  ! primary arm spacing (dendritic only)
    real(r8) :: lambda2 = 0.0_r8  ! secondary arm spacing (dendritic only)
    integer  :: ustruc = GV_INVALID
  end type

  type, extends(ustruc_comp) :: ustruc_ldrd
    real(r8) :: f1, f1p, f2, f2p
    integer,  allocatable :: state(:)
    real(r8), allocatable :: dt(:)
    type(gv) :: ustruc_model
    type(gv_state), allocatable :: ustruc_state(:)
  contains
    procedure :: set_state
    procedure :: update_state
    procedure :: get_comp_list
    procedure :: has
    procedure :: getl1
    procedure :: getr1
    procedure :: serialize
    procedure :: deserialize
  end type ustruc_ldrd

  integer, parameter :: STATE_INVALID   = 0
  integer, parameter :: STATE_UNDEFINED = 1
  integer, parameter :: STATE_LIQUID    = 2
  integer, parameter :: STATE_MUSHY     = 3
  integer, parameter :: STATE_SOLID     = 4

#ifndef INTEL_BUG20210619
  !! Number of bytes (per cell) of internal state for serialization/deserialization
  type(ustruc_ldrd), allocatable :: dummy  ! only use is in the following parameter declaration
  integer, parameter :: NBYTES = storage_size(dummy%ustruc_state%G)/8 + &
                                 storage_size(dummy%ustruc_state%V)/8 + &
                                 storage_size(dummy%ustruc_state%sfrac)/8 +  &
                                 storage_size(dummy%ustruc_state%lambda1)/8 + &
                                 storage_size(dummy%ustruc_state%lambda2)/8 + &
                                 storage_size(dummy%ustruc_state%ustruc)/8 + &
                                 storage_size(dummy%state)/8 + &
                                 storage_size(dummy%dt)/8
#endif

contains

  subroutine gv_init(this, params)
    use parameter_list_type
    use truchas_logging_services
    class(gv), intent(out) :: this
    type(parameter_list) :: params
    integer :: stat
    character(:), allocatable :: errmsg
    call params%get('theta-gv', this%theta, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
    call params%get('liquidus-slope', this%m_liq, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
    if(this%m_liq >= 0.0_r8) call TLS_fatal('liquidus-slope must be < 0')
    call params%get('solute-conc', this%c_0, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
    call params%get('partition-coef', this%k_part, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
    call params%get('liq-sol-delta-T', this%delta_T, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
    call params%get('diffusivity', this%d, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
    call params%get('gibbs-thomson-coef', this%gamma, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
    call params%get('instability-coef', this%alpha, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
    call params%get('coarsening-coef', this%k_coarse, stat, errmsg)
    if (stat /= 0) call TLS_fatal(errmsg)
  end subroutine

  !! This method is called when transitioning from the liquid to mushy state.
  subroutine gv_start(this, sfrac, G, L, state)
    class(gv), intent(in) :: this
    real(r8),  intent(in) :: sfrac, G, L
    type(gv_state), intent(out) :: state
    state%ustruc = GV_UNDEFINED
    call gv_update(this, sfrac, G, L, state)
  end subroutine

  !! This method is called when continuing in the mushy state.
  subroutine gv_update(this, sfrac, G, L, state)

    class(gv), intent(in) :: this
    real(r8),  intent(in) :: sfrac, G, L
    type(gv_state), intent(inout) :: state

    real(r8) :: c, lambda1

    if (state%ustruc == GV_UNDEFINED) then
      if (sfrac >= this%theta) then
        if (L >= huge(1.0d0)*min(1.0d0, G)) then ! velocity overflows
          state%ustruc = GV_INVALID
        else
          state%G = G
          state%V = L/G
          c = -this%m_liq * this%c_0 * (1 - this%k_part) / (this%k_part * this%d)
          if (state%G > c*state%V) then
            state%ustruc = GV_PLANAR
          else
            lambda1 = this%lambda1(state%G,state%V)
            if (lambda1 < this%lambdaCD(state%G,state%V)) then
              state%ustruc = GV_CELLULAR
            else
              state%ustruc = GV_DENDRITIC
              state%lambda1 = lambda1
            end if
          end if
        end if
      end if
    end if

  end subroutine gv_update

  !! This method is called when making the transition from mushy to solid.
  subroutine gv_finish(this, dt, state)
    class(gv), intent(in) :: this
    real(r8), intent(in) :: dt
    type(gv_state), intent(inout) :: state
    if (state%ustruc == GV_UNDEFINED) then
      state%ustruc = GV_INVALID
    elseif (state%ustruc == GV_DENDRITIC) then
      state%lambda2 = this%k_coarse*sqrt(dt)
    end if
  end subroutine

  !! This method is called when making a transition other than above.
  subroutine gv_reset(this, state)
    class(gv), intent(in) :: this
    type(gv_state), intent(out) :: state  ! default initialized
  end subroutine

  pure function gv_lambda1(this, G, V) result(lambda1)
    class(gv), intent(in) :: this
    real(r8), intent(in) :: G, V
    real(r8) :: lambda1, t1, t2
    t1 = (6*this%delta_T/((1-this%k_part)*G))*((this%d/V) - (this%delta_T*this%k_part/G))
    t1 = sqrt(max(0.0_r8, t1))
    t2 = 4.3_r8*sqrt(sqrt(this%d*this%gamma/(this%delta_T*this%k_part*V))*this%delta_T/G)
    lambda1 = max(t1,t2)
  end function

  pure function gv_lambdaCD(this, G, V) result(lambdaCD)
    class(gv), intent(in) :: this
    real(r8), intent(in) :: G, V
    real(r8) :: lambdaCD
    lambdaCD = (this%alpha/this%c_0**0.25_r8)*(this%d*this%gamma/(G*V))**(1.0_r8/3.0_r8)
  end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function new_ustruc_ldrd(comp, params) result(this)

    use parameter_list_type
    use truchas_logging_services

    class(ustruc_analysis), pointer, intent(in) :: comp
    type(parameter_list) :: params
    type(ustruc_ldrd), pointer :: this
    integer :: stat
    character(:), allocatable :: errmsg

    allocate(this)
    call this%init(comp)

    call params%get('theta1', this%f1, stat, errmsg)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if(this%f1 <= 0.0 .or. this%f1 >= 1.0) then
      call TLS_fatal('theta1 must be > 0.0 and < 1.0')
    end if

    call params%get('theta2', this%f2, stat, errmsg)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%f2 <= 0.0 .or. this%f2 >= 1.0) then
      call TLS_fatal('theta2 must be > 0.0 and < 1.0')
    end if
    if (this%f2 <= this%f1) call TLS_fatal('theta2 <= theta1')

    call params%get('theta1p', this%f1p, stat, errmsg, default=this%f1)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%f1p > this%f1 .or. this%f1p < 0.0) then
      call TLS_fatal('theta1p must be >= 0.0 and <= theta1')
    end if

    call params%get('theta2p', this%f2p, stat, errmsg, default=this%f2)
    if (stat /= 0) then
      call TLS_fatal(errmsg)
    else if (this%f1p > this%f1 .or. this%f1p < 0.0) then
      call TLS_fatal('theta2p must be >= theta1 and <= theta2')
    end if

    call this%ustruc_model%init(params)

    if (this%ustruc_model%theta < this%f1 .or. this%ustruc_model%theta > this%f2) then
      call TLS_fatal('theta-gv must be >= theta1 and <= theta2')
    end if

    allocate(this%state(this%n), this%dt(this%n), this%ustruc_state(this%n))

  end function new_ustruc_ldrd

  subroutine set_state(this, t, temp, temp_grad, frac, invalid)

    class(ustruc_ldrd), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    integer :: j

    !! Set the initial state for the next analysis component in the chain.
    !! This will define the CORE component state arrays.
    call this%ustruc_comp%set_state(t, temp, temp_grad, frac, invalid)

    do j = 1, this%n
      if (invalid(j)) then
        this%state(j) = STATE_INVALID
      else if (frac(j) <= this%f1) then
        this%state(j) = STATE_LIQUID
      else
        this%state(j) = STATE_UNDEFINED
      end if
    end do

    !! Assign dummy values to the remaining arrays.
    this%dt = 0.0_r8

  end subroutine set_state

  subroutine update_state(this, t, temp, temp_grad, frac, invalid)

    class(ustruc_ldrd), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:)
    logical,  intent(in) :: invalid(:)

    integer  :: j
    real(r8) :: G, L

    !! Update the state of this analysis component.
    associate (prev_t => this%core%t, prev_temp => this%core%temp, prev_frac => this%core%frac)
      do j = 1, this%n
        if (invalid(j)) then
          this%state(j) = STATE_INVALID
          call this%ustruc_model%reset(this%ustruc_state(j))
        else
          select case (this%state(j))
          case (STATE_LIQUID)
            if (frac(j) >= this%f1) then
              this%dt(j) = interp_frac(prev_t, prev_frac(j), t, frac(j), this%f1)
              this%state(j) = STATE_MUSHY
            end if
            if (frac(j) > this%f2) then
              this%dt(j) = interp_frac(prev_t, prev_frac(j), t, frac(j), this%f2) - this%dt(j)
              this%state(j) = STATE_SOLID
            end if
            if (this%state(j) == STATE_MUSHY) then
              G = this%vector_magnitude(temp_grad(:,j))
              L = (prev_temp(j) - temp(j))/(t - prev_t)
              call this%ustruc_model%start(frac(j), G, L, this%ustruc_state(j))
            end if
          case (STATE_MUSHY)
            if (frac(j) < this%f1p) then
              this%state(j) = STATE_LIQUID
              call this%ustruc_model%reset(this%ustruc_state(j))
            else if (frac(j) > this%f2) then
              this%dt(j) = interp_frac(prev_t, prev_frac(j), t, frac(j), this%f2) - this%dt(j)
              this%state(j) = STATE_SOLID
              call this%ustruc_model%finish(this%dt(j), this%ustruc_state(j))
            else
              G = this%vector_magnitude(temp_grad(:,j))
              L = (prev_temp(j) - temp(j))/(t - prev_t)
              call this%ustruc_model%update(frac(j), G, L, this%ustruc_state(j))
            end if
          case (STATE_SOLID)
            if (frac(j) < this%f1) then
              this%state(j) = STATE_LIQUID
            else if (frac(j) <= this%f2p) then
              this%state(j) = STATE_UNDEFINED
            end if
            if (this%state(j) /= STATE_SOLID) call this%ustruc_model%reset(this%ustruc_state(j))
          case (STATE_UNDEFINED)
            if (frac(j) < this%f1) this%state(j) = STATE_LIQUID
          case (STATE_INVALID)
            if (frac(j) < this%f1) then
              this%state(j) = STATE_LIQUID
            else
              ! might be able to sensibly identify microstructure in this case too
              this%dt(j) = interp_frac(prev_t, 0.0_r8, t, frac(j), this%f1)
              this%state(j) = STATE_MUSHY
              if (frac(j) > this%f2) then
                this%dt(j) = interp_frac(prev_t, 0.0_r8, t, frac(j), this%f2) - this%dt(j)
                this%state(j) = STATE_SOLID
              end if
            end if
          case default
            INSIST(.false.)
          end select
        end if
      end do
    end associate

    !! Update the next analysis component in the chain.
    !! Note that the core analysis component always ends the chain.
    call this%ustruc_comp%update_state(t, temp, temp_grad, frac, invalid)

  contains

    pure function interp_frac(t1, f1, t2, f2, f) result(t)
      real(r8), intent(in) :: t1, f1, t2, f2, f
      real(r8) :: t
      t = t1*((f2-f)/(f2-f1)) + t2*((f-f1)/(f2-f1))
    end function

  end subroutine update_state

  subroutine get_comp_list(this, list)
    class(ustruc_ldrd), intent(in) :: this
    integer, allocatable, intent(out) :: list(:)
    integer, allocatable :: rest(:)
    call this%ustruc_comp%get_comp_list(rest)
    allocate(list(size(rest)+1))
    list(1) = USTRUC_LDRD_ID
    list(2:) = rest
  end subroutine

  logical function has(this, name)
    class(ustruc_ldrd), intent(in) :: this
    character(*), intent(in) :: name
    select case (name)
    case ('ldrd-invalid-GV', 'ldrd-t_sol', 'ldrd-type', 'ldrd-lambda1', 'ldrd-lambda2', 'ldrd-G', 'ldrd-V')
      has = .true.
    case default
      has = this%ustruc_comp%has(name)
    end select
  end function

  subroutine getl1(this, name, array)
    class(ustruc_ldrd), intent(in) :: this
    character(*), intent(in) :: name
    logical, intent(out) :: array(:)
    select case (name)
    case ('ldrd-invalid-GV')
      ASSERT(size(array) == this%n)
      array = (this%state /= STATE_SOLID)
    case default
      call this%ustruc_comp%get(name, array)
    end select
  end subroutine

  subroutine getr1(this, name, array, invalid)
    class(ustruc_ldrd), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    select case (name)
    case ('ldrd-t_sol')
      ASSERT(size(array) == this%n)
      where (this%state == STATE_SOLID)
        array = this%dt
      else where
        array = 0.0_r8
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%state /= STATE_SOLID)
      end if
    case ('ldrd-type')
      ASSERT(size(array) == this%n)
      array = this%ustruc_state%ustruc
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%state /= STATE_SOLID .or. this%ustruc_state%ustruc == GV_UNDEFINED &
                                             .or. this%ustruc_state%ustruc == GV_INVALID)
      end if
    case ('ldrd-lambda1')
      ASSERT(size(array) == this%n)
      where (this%state == STATE_SOLID .and. this%ustruc_state%ustruc == GV_DENDRITIC)
        array = this%ustruc_state%lambda1
      else where
        array = 0.0_r8
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%state /= STATE_SOLID .or. this%ustruc_state%ustruc /= GV_DENDRITIC)
      end if
    case ('ldrd-lambda2')
      ASSERT(size(array) == this%n)
      where (this%state == STATE_SOLID .and. this%ustruc_state%ustruc == GV_DENDRITIC)
        array = this%ustruc_state%lambda2
      else where
        array = 0.0_r8
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%state /= STATE_SOLID .or. this%ustruc_state%ustruc /= GV_DENDRITIC)
      end if
    case ('ldrd-G')
      ASSERT(size(array) == this%n)
      where (this%ustruc_state%ustruc > GV_UNDEFINED)
        array = this%ustruc_state%G
      else where
        array = 0.0_r8
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%ustruc_state%ustruc <= GV_UNDEFINED)
      end if
    case ('ldrd-V')
      ASSERT(size(array) == this%n)
      where (this%ustruc_state%ustruc > GV_UNDEFINED)
        array = this%ustruc_state%V
      else where
        array = 0.0_r8
      end where
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%ustruc_state%ustruc <= GV_UNDEFINED)
      end if
    case default
      call this%ustruc_comp%get(name, array, invalid)
    end select
  end subroutine getr1

  subroutine serialize(this, cid, array)

    use,intrinsic :: iso_fortran_env, only: int8
    use serialization_tools, only: copy_to_bytes

    class(ustruc_ldrd), intent(in) :: this
    integer, intent(in) :: cid
    integer(int8), allocatable, intent(out) :: array(:,:)

    integer :: j, offset
#ifdef INTEL_BUG20210619
    integer :: NBYTES
    NBYTES = storage_size(this%ustruc_state%G)/8 + &
             storage_size(this%ustruc_state%V)/8 + &
             storage_size(this%ustruc_state%sfrac)/8 +  &
             storage_size(this%ustruc_state%lambda1)/8 + &
             storage_size(this%ustruc_state%lambda2)/8 + &
             storage_size(this%ustruc_state%ustruc)/8 + &
             storage_size(this%state)/8 + &
             storage_size(this%dt)/8
#endif

    if (cid == USTRUC_LDRD_ID) then
      allocate(array(NBYTES,this%n))
      do j = 1, this%n
        offset = 0
        call copy_to_bytes(this%ustruc_state(j)%G, array(:,j), offset)
        call copy_to_bytes(this%ustruc_state(j)%V, array(:,j), offset)
        call copy_to_bytes(this%ustruc_state(j)%sfrac, array(:,j), offset)
        call copy_to_bytes(this%ustruc_state(j)%lambda1, array(:,j), offset)
        call copy_to_bytes(this%ustruc_state(j)%lambda2, array(:,j), offset)
        call copy_to_bytes(this%ustruc_state(j)%ustruc, array(:,j), offset)
        call copy_to_bytes(this%state(j), array(:,j), offset)
        call copy_to_bytes(this%dt(j), array(:,j), offset)
      end do
    else
      call this%ustruc_comp%serialize(cid, array)
    end if

  end subroutine serialize

  subroutine deserialize(this, cid, array)

    use,intrinsic :: iso_fortran_env, only: int8
    use serialization_tools, only: copy_from_bytes

    class(ustruc_ldrd), intent(inout) :: this
    integer, intent(in) :: cid
    integer(int8), intent(in) :: array(:,:)

    integer :: j, offset
#ifdef INTEL_BUG20210619
    integer :: NBYTES
    NBYTES = storage_size(this%ustruc_state%G)/8 + &
             storage_size(this%ustruc_state%V)/8 + &
             storage_size(this%ustruc_state%sfrac)/8 +  &
             storage_size(this%ustruc_state%lambda1)/8 + &
             storage_size(this%ustruc_state%lambda2)/8 + &
             storage_size(this%ustruc_state%ustruc)/8 + &
             storage_size(this%state)/8 + &
             storage_size(this%dt)/8
#endif

    if (cid == USTRUC_LDRD_ID) then
      INSIST(size(array,1) == NBYTES)
      INSIST(size(array,2) == this%n)
      do j = 1, this%n
        offset = 0
        call copy_from_bytes(array(:,j), offset, this%ustruc_state(j)%G)
        call copy_from_bytes(array(:,j), offset, this%ustruc_state(j)%V)
        call copy_from_bytes(array(:,j), offset, this%ustruc_state(j)%sfrac)
        call copy_from_bytes(array(:,j), offset, this%ustruc_state(j)%lambda1)
        call copy_from_bytes(array(:,j), offset, this%ustruc_state(j)%lambda2)
        call copy_from_bytes(array(:,j), offset, this%ustruc_state(j)%ustruc)
        call copy_from_bytes(array(:,j), offset, this%state(j))
        call copy_from_bytes(array(:,j), offset, this%dt(j))
      end do
    else
      call this%ustruc_comp%deserialize(cid, array)
    end if

  end subroutine deserialize

end module ustruc_ldrd_type

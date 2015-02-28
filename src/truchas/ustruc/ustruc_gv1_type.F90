!!
!! USTRUC_GV1
!!
!! A concrete implementation of USTRUC_PLUGIN/USTRUC_COMP that adds the
!! identification of microstructure type and characteristics (first model).
!!
!! NB: THIS IS AN INCOMPLETE STUB.  ONLY THE THERMAL GRADIENT MAGNITUDE (G)
!! AND SOLIDIFICATION FRONT SPEED (V) ARE COMPUTED.  THESE WILL BE THE INPUTS
!! TO A TABLE LOOKUP FOR THE TYPE OF MICROSTRUCTURE (DENDRITIC, PLANAR, ...)
!! AND THE COMPUTATION OF CHARACTERISTICS (PRIMARY AND SECONDARY ARM SPACING)
!! AWAITING THIS INFO FROM SETH IMHOFF AND PAUL GIBBS (MST-6)
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! September 2014
!!
!! PROGRAMMING INTERFACE
!!
!!  See the comments that accompany the source for the USTRUC_COMP and
!!  USTRUC_PLUGIN classes which define the interface to this object.  The
!!  only public entity provided by the module is the following function.
!!
!!  NEW_USTRUC_GV1(COMP, PARAMS) returns a pointer to a new USTRUC_COMP class
!!    object whose dynamic type is USTRUC_GV1.  The USTRUC_GV1 type is itself
!!    private.  COMP is a USTRUC_COMP pointer whose target is being wrapped
!!    by this new analysis component specified by PARAMS.  The new object
!!    takes ownership of the target.  PARAMS is a PARAMETER_LIST object;
!!    the relevant parameters in PARAMS are:
!!
!!      'theta1' -- low solid fraction threshold; solidification is deemed
!!          to have started when the solid fraction crosses this threshold.
!!      'theta2' -- high solid fraction threshold; solidification is deemed
!!          to have finished when the solid fraction crosses this threshold.
!!      'theta1p' -- if the solid fraction drops below this threshold while
!!          solidifying, it is deemed to have returned to a liquid state;
!!          theta1p <= theta1, optional, default theta1p = theta1.
!!      'theta2p' -- if the solid fraction drops below this threshold while
!!          considered solid, it is considered to have remelted and the
!!          previously computed solidification time erased.
!!
!!  Objects of this type respond to the following data names in the generic
!!  GET subroutine: 'solid-time', 'g', and 'v'.
!!

#include "f90_assert.fpp"

module ustruc_gv1_type

  use kinds, only: r8
  use ustruc_plugin_class
  implicit none
  private

  public :: new_ustruc_gv1

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

  type, extends(ustruc_plugin) :: ustruc_gv1
    real(r8) :: f1, f1p, f2, f2p
    integer,  allocatable :: state(:)
    real(r8), allocatable :: dt(:)
    type(gv) :: ustruc_model
    type(gv_state), allocatable :: ustruc_state(:)
  contains
    procedure :: set_state
    procedure :: update_state
    procedure :: getl1
    procedure :: getr1
  end type ustruc_gv1

  integer, parameter :: STATE_INVALID   = 0
  integer, parameter :: STATE_UNDEFINED = 1
  integer, parameter :: STATE_LIQUID    = 2
  integer, parameter :: STATE_MUSHY     = 3
  integer, parameter :: STATE_SOLID     = 4

contains

  subroutine gv_init (this, params)
    use parameter_list_type
    class(gv), intent(out) :: this
    type(parameter_list) :: params
    call params%get ('liquidus-slope', this%m_liq)
    call params%get ('solute-conc', this%c_0)
    call params%get ('partition-coef', this%k_part)
    call params%get ('liq-sol-delta-T', this%delta_T)
    call params%get ('diffusivity', this%d)
    call params%get ('gibbs-thomson-coef', this%gamma)
    call params%get ('instability-coef', this%alpha)
    call params%get ('coarsening-coef', this%k_coarse)
  end subroutine gv_init
  
  !! This method is called when transitioning from the liquid to mushy state.
  subroutine gv_start (this, sfrac, G, V, invalid_V, state)
    class(gv), intent(in) :: this
    real(r8),  intent(in) :: sfrac, G, V
    logical,   intent(in) :: invalid_V
    type(gv_state), intent(out) :: state
    state%ustruc = GV_UNDEFINED
    call gv_update (this, sfrac, G, V, invalid_V, state)
  end subroutine gv_start
  
  !! This method is called when continuing in the mushy state.
  subroutine gv_update (this, sfrac, G, V, invalid_V, state)
  
    class(gv), intent(in) :: this
    real(r8),  intent(in) :: sfrac, G, V
    logical,   intent(in) :: invalid_V
    type(gv_state), intent(inout) :: state
    
    real(r8) :: c, lambda1
    
    if (state%ustruc == GV_UNDEFINED) then
      if (sfrac >= this%theta) then
        if (invalid_V) then
          state%ustruc = GV_INVALID ! cannot handle this
        else
          state%G = G
          state%V = V
          c = -this%m_liq * this%c_0 * (1 - this%k_part) / (this%k_part * this%d)
          if (state%G < c*state%V) then
            state%ustruc = GV_PLANAR
          else
            lambda1 = this%lambda1(G,V)
            if (lambda1 < this%lambdaCD(G,V)) then
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
  subroutine gv_finish (this, dt, state)
    class(gv), intent(in) :: this
    real(r8), intent(in) :: dt
    type(gv_state), intent(inout) :: state
    if (state%ustruc == GV_UNDEFINED) then
      state%ustruc = GV_INVALID
    elseif (state%ustruc == GV_DENDRITIC) then
      state%lambda2 = this%k_coarse*sqrt(dt)
    else
      INSIST(.false.)
    end if
  end subroutine gv_finish
  
  !! This method is called when making a transition other than above.
  subroutine gv_reset (this, state)
    class(gv), intent(in) :: this
    type(gv_state), intent(out) :: state  ! default initialized
  end subroutine gv_reset
  
  pure function gv_lambda1 (this, G, V) result (lambda1)
    class(gv), intent(in) :: this
    real(r8), intent(in) :: G, V
    real(r8) :: lambda1, t1, t2
    t1 = (6*this%delta_T/((1-this%k_part)*G))*((this%d/V) - (this%delta_T*this%k_part/G))
    t1 = sqrt(max(0.0_r8, t1))
    t2 = 4.3_r8*sqrt(sqrt(this%d*this%gamma/(this%delta_T*this%k_part*V))*this%delta_T/G)
    lambda1 = max(t1,t2)
  end function gv_lambda1
  
  pure function gv_lambdaCD (this, G, V) result (lambdaCD)
    class(gv), intent(in) :: this
    real(r8), intent(in) :: G, V
    real(r8) :: lambdaCD
    lambdaCD = (this%alpha/this%c_0**0.25_r8)*(this%d*this%gamma/(G*V))**(1.0_r8/3.0_r8)
  end function gv_lambdaCD
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function new_ustruc_gv1 (comp, params) result (this)

    use parameter_list_type
    use parameter_list_json

    class(ustruc_comp), pointer, intent(in) :: comp
    type(parameter_list) :: params
    type(ustruc_gv1), pointer :: this
    
    character(:), allocatable :: filename, errmsg
    type(parameter_list), pointer :: plist
    integer :: lun

    allocate(this)
    call this%init (comp)

    call params%get ('theta1',  this%f1)
    INSIST(this%f1 >= 0.0_r8)
    call params%get ('theta2',  this%f2)
    INSIST(this%f2 <= 1.0_r8)
    INSIST(this%f1 <= this%f2)
    call params%get ('theta1p', this%f1p, default=this%f1)
    INSIST(this%f1p >= 0.0_r8 .and. this%f1p <= this%f1)
    call params%get ('theta2p', this%f2p, default=this%f2)
    INSIST(this%f2p >= this%f1 .and. this%f2p <= this%f2)
    
    call params%get ('parameter-file', filename)
    
    !FIXME: every process is reading -- better to read on one and broadcast(?)
    open(newunit=lun,file=filename,action='read',access='stream',form='unformatted')
    call parameter_list_from_json_stream (lun, plist, errmsg)
    INSIST(associated(plist))
    
    call this%ustruc_model%init (plist)
    deallocate(plist)

    allocate(this%state(this%n), this%dt(this%n), this%ustruc_state(this%n))

  end function new_ustruc_gv1

  subroutine set_state (this, t, temp, temp_grad, frac, frac_grad, invalid)

    class(ustruc_gv1), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:), frac_grad(:,:)
    logical,  intent(in) :: invalid(:)

    integer :: j

    !! Set the initial state for the next analysis component in the chain.
    !! This will define the CORE component state arrays.
    call this%ustruc_plugin%set_state (t, temp, temp_grad, frac, frac_grad, invalid)

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

  subroutine update_state (this, t, temp, temp_grad, frac, frac_grad, invalid)

    class(ustruc_gv1), intent(inout) :: this
    real(r8), intent(in) :: t, temp(:), temp_grad(:,:), frac(:), frac_grad(:,:)
    logical,  intent(in) :: invalid(:)

    integer  :: j
    real(r8) :: prev_t, prev_frac(this%n), G, V

    prev_t = this%core%t
    prev_frac = this%core%frac

    !! Update the next analysis component in the chain.
    !! Note that the core analysis component always ends the chain.
    call this%ustruc_plugin%update_state (t, temp, temp_grad, frac, frac_grad, invalid)

    !! Update the state of this analysis component.
    do j = 1, this%n
      if (invalid(j)) then
        this%state(j) = STATE_INVALID
        call this%ustruc_model%reset (this%ustruc_state(j))
      else
        associate (temp_grad => this%core%temp_grad, velocity => this%core%velocity, &
                   invalid_velocity => this%core%invalid_velocity)
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
              V = this%vector_magnitude(velocity(:,j))
              call this%ustruc_model%start (frac(j), G, V, invalid_velocity(j), this%ustruc_state(j))
            end if
          case (STATE_MUSHY)
            if (frac(j) < this%f1p) then
              this%state(j) = STATE_LIQUID
              call this%ustruc_model%reset (this%ustruc_state(j))
            else if (frac(j) > this%f2) then
              this%dt(j) = interp_frac(prev_t, prev_frac(j), t, frac(j), this%f2) - this%dt(j)
              this%state(j) = STATE_SOLID
              call this%ustruc_model%finish (this%dt(j), this%ustruc_state(j))
            else
              call this%ustruc_model%update (frac(j), G, V, invalid_velocity(j), this%ustruc_state(j))
            end if
          case (STATE_SOLID)
            if (frac(j) < this%f1) then
              this%state(j) = STATE_LIQUID
            else if (frac(j) <= this%f2p) then
              this%state(j) = STATE_UNDEFINED
            end if
            if (this%state(j) /= STATE_SOLID) call this%ustruc_model%reset (this%ustruc_state(j))
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
        end associate
      end if
    end do

  contains

    pure function interp_frac (t1, f1, t2, f2, f) result (t)
      real(r8), intent(in) :: t1, f1, t2, f2, f
      real(r8) :: t
      t = t1*((f-f1)/(f2-f1)) + t2*((f2-f)/(f2-f1))
    end function

  end subroutine update_state

  subroutine getl1 (this, name, array)
    class(ustruc_gv1), intent(in) :: this
    character(*), intent(in) :: name
    logical, intent(out) :: array(:)
    select case (name)
    case ('invalid-gv')
      ASSERT(size(array) == this%n)
      array = (this%state /= STATE_SOLID)
    case default
      call this%ustruc_plugin%get (name, array)
    end select
  end subroutine getl1

  subroutine getr1 (this, name, array, invalid)
    class(ustruc_gv1), intent(in) :: this
    character(*), intent(in) :: name
    real(r8), intent(out) :: array(:)
    logical, intent(out), optional :: invalid(:)
    select case (name)
    case ('solid-time')
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
    case ('ustruc')
      ASSERT(size(array) == this%n)
      array = this%ustruc_state%ustruc
      if (present(invalid)) then
        ASSERT(size(invalid) == this%n)
        invalid = (this%state /= STATE_SOLID .or. this%ustruc_state%ustruc == GV_UNDEFINED &
                                             .or. this%ustruc_state%ustruc == GV_INVALID)
      end if
    case ('lambda1')
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
    case ('lambda2')
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
    case default
      call this%ustruc_plugin%get (name, array, invalid)
    end select
  end subroutine getr1

end module ustruc_gv1_type

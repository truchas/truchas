!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This module defines the type AM_STATE.
!! Objects of type AM_STATE are responsible for 
!! for holding and manipulating the state data: temperatures,
!! enthalpies, velocities, denisties, volume fractions, and
!! material information such as all material indices corresponding
!! to fluid materials
!!
!! PROGRAMMING INTERFACE
!!
!!  The type AM_STATE has the following type bound procedures:
!!     ALLOC_AND_INIT(THIS, VOL_FRACS, TEMPS, &
!!                    DENSITIES, VELOCITIES_X, &
!!                    VELOCITIES_Y, VELOCITIES_Z)
!!      Allocates arrays for holding e.g. volume fractions, temperatures, 
!!      enthalpies, velocities, and densities, as well as an array
!!      for storing which material indices correspond to fluids, and initializes these arrays
!!      with those from the argument list. Also computes and stores material indices correspond to fluids
!!
!!     SET_STATE(THIS, VOL_FRACS, TEMPS, &
!!                    DENSITIES, VELOCITIES_X, &
!!                    VELOCITIES_Y, VELOCITIES_Z)
!!      Sets volume fractions, temperatures, enthalpies, velocities, and densities; also
!!      calculates the enthalpy densities from the temperatures. NOTE: It can happen that
!!      temperatures that are passed in are negative, since the diffusion solver does can
!!      return negative temperatures (invariably in cells that are mostly void but 
!!      contain tiny fluid fragments); for such cells n, THIS%ignore_cells(n) = .true., 
!!      which has the consequence that the source term dQ(:) formed for the diffusion solver is
!!      set to zero in this cell (dQ(n) = 0). In addition, this cell is ignored for the purposes
!!      of energy deposition.
!!
!!     COMPUTE_ENTHALPY_DENSITY (THIS, N, H)
!!      Computes enthalpy denisty H in cell n from THIS%temps(n) in cell n; 
!!      uses volume fractions THIS%vol_fracs(1:nmat,n) for this computation
!!
!!     COMPUTE_TEMP (this, n, H, T)
!!      Computes temperature T in cell n from the enthalpy
!       density H; uses volume fractions THIS%vol_fracs(1:nmat,n)
!!      for this computation (hence the cell argument n). NOTE: this computation
!!      requires solving h(T) = H, where h(T) gives the enthalpy denisty as a function
!!      of the temperature; we use the object THIS%root_solver of type ROOT_FINDER for this purpose
!!
!!      FIXUP_ENERGY_IN_NEARLY_VOID_CELLS(THIS, IS_NEARLY_VOID)
!!       In cells n that are nearly void (is_nearly_void(n) = .true.), set the temperature
!!       of the cell to be the average of its neighbors (only positive neighboring temperatures
!!       contribute to the average); set the corresponding enthalpy density so that it is
!!       consistent with its new temperature value. NOTE: in rare cases, a cell can have a 
!!       negative temperature and be surrounded by all void cells (invariably a cell with a
!!       tiny amount of fluid); currently this routine sets the temperarure in this cell to
!!       be 300 (arbitrarily)---need to change this to something more reasonable, e.g. setting
!!       the temperature to its previous value.

module am_state_type

  use kinds, only: r8
  use root_finder_type
  use truchas_logging_services


  type, public:: am_state
    integer :: void_index
    !!!! TODO: CHANGE TO READ THIS IN !!!!!
    real(r8) :: one_tol = 10.0_r8**(-8)
    real(r8) :: am_tol = 10.0_r8**(-3)
    !!!! TODO: CHANGE TO READ THIS IN !!!!
    integer :: nfaces, ncells, nmat, nfluids
    integer, allocatable :: fluid_indices(:)
    real(r8), allocatable :: vol_fracs(:,:) ! material volume fractions
    real(r8), allocatable :: void_vol_fracs(:) ! void volume fractions
    real(r8), allocatable :: void_vol_fracs_old(:) ! previous time step
    real(r8), allocatable :: temps(:)
    ! holds the enthalpy per unit volume for each cell
    real(r8), allocatable :: enthalpies(:)
    real(r8), allocatable :: enthalpies_old(:), temps_old(:) ! from previous time step
    real(r8), allocatable :: densities(:)
    real(r8), allocatable :: velocities_x(:), velocities_y(:), velocities_z(:)
    ! if a non-void cell has negative temperature (from the diffusion solve), 
    ! then ignore this
    logical, allocatable :: ignore_cell(:)
    ! root_solver is used in  compute_temp_from_enthalpy_density to
    ! solve for H(T) = H_0, where T = temp and H = enthalpy density
    type(root_finder) :: root_solver
  contains
    procedure :: alloc_and_init
    procedure :: set_state
    procedure :: compute_enthalpy_density
    procedure :: compute_temp => compute_temp_from_enthalpy_density
    procedure :: fixup_energy_in_nearly_void_cells
  end type am_state

  contains


    subroutine alloc_and_init(this, vol_fracs, temps, &
                              densities, velocities_x, &
                              velocities_y, velocities_z)

      use parameter_module, only: nmat, ncells, nfc
      use fluid_data_module, only: Void_Material_Index, isImmobile

      class(am_state), intent(inout) :: this
      real(r8), intent(in) :: vol_fracs(:,:), temps(:), &
                              densities(:), velocities_x(:), velocities_y(:), &
                              velocities_z(:)
    
      integer :: status, m, m_fluid, n
      real(r8) :: H

      ! Assume there is only one void material; extract the void material index
      this%void_index = Void_Material_Index(1)

      ! Set the number of faces, number of cells, and number of materials
      ! TODO: note that these are imported from parameter_module, pass as parameters instead
      this%ncells = ncells
      this%nfaces = nfc
      this%nmat = nmat


      ! Get the number of fluid materials
      m_fluid = 0
      do  m = 1, nmat
        if ( .not. isImmobile(m) .and. m /= this%void_index ) then
          m_fluid = m_fluid + 1
        end if
      end do
      this%nfluids = m_fluid

      ! Allocate array to hold indices of interface cells
      allocate(   this%vol_fracs(nmat,ncells), &
                  this%void_vol_fracs(ncells), &
                  this%void_vol_fracs_old(ncells), &
                  this%temps(ncells), &
                  this%temps_old(ncells), &
                  this%enthalpies(ncells), &
                  this%enthalpies_old(ncells), &
                  this%densities(ncells), &
                  this%velocities_x(ncells), &
                  this%velocities_y(ncells), &
                  this%velocities_z(ncells), &
                  this%fluid_indices(m_fluid), &
                  this%ignore_cell(ncells), &
                  STAT = status   )
      if (status /= 0) call TLS_fatal ('alloc_am_state_vars: allocation failed') 


      ! Initialize the array holding the material indices that are fluid
      m_fluid = 0
      do  m = 1, nmat
        if ( .not. isImmobile(m) .and. m /= this%void_index ) then
          m_fluid = m_fluid + 1
          this%fluid_indices(m_fluid) = m
        end if
      end do


      ! Initialize the root solver object (this is used for
      ! computing the enthalpy density from the temperature and 
      ! volume fractions
      !!!!!! TODO: CHANGE---READ THIS IN FROM INPUT !!!!!!!!
      call this%root_solver%init(100, 100, 100, 10.0d0**(-12))
      !!!!!! CHANGE !!!!!!!!     

      ! Initialize  the state variables
      this%vol_fracs = vol_fracs
      this%void_vol_fracs = vol_fracs(this%void_index, :)
      this%temps = temps
      this%densities = densities
      this%velocities_x = velocities_x
      this%velocities_y = velocities_y
      this%velocities_z = velocities_z
      ! Initially don't ignore any cells
      this%ignore_cell = .false.

      ! Compute the enthalpy density corresponding to the current temperature
      ! if the cell is not void
      do n = 1, ncells
        if ( abs(this%void_vol_fracs(n) - 1) < this%one_tol  ) then
          this%enthalpies(n) = 0
        else
          call this%compute_enthalpy_density(n, H)
          this%enthalpies(n) = H
        end if
      end do


      ! Set old enthalpy densities and volume fractions
      this%enthalpies_old = this%enthalpies
      this%void_vol_fracs_old = this%void_vol_fracs


    end subroutine alloc_and_init


    subroutine set_state(this, vol_fracs, temps, &
                         densities, velocities_x, &
                         velocities_y, velocities_z)

      class(am_state), intent(inout) :: this
      real(r8), intent(in) :: vol_fracs(:,:), temps(:), &
                              densities(:), velocities_x(:), &
                              velocities_y(:), velocities_z(:)

      integer :: n, ncells
      real(r8) :: H, one_tol


      ncells = this%ncells
      one_tol = this%one_tol

      ! initialize state vars
      this%vol_fracs = vol_fracs
      this%void_vol_fracs = vol_fracs(this%void_index, :)
      this%temps = temps
      this%densities = densities
      this%velocities_x = velocities_x
      this%velocities_y = velocities_y
      this%velocities_z = velocities_z


      ! Compute the enthalpy density corresponding to the current temperature
      ! if the cell is not void and was not void on the previous time step 
      ! (if void on the previous time step, the enthalpy will be updated based
      ! on the energy fluxed into the cell from the advection step)
      this%ignore_cell = .false.
      do n = 1, ncells

        if ( abs(this%void_vol_fracs(n) - 1) < one_tol .or. &
             abs(this%void_vol_fracs_old(n) - 1) < one_tol  ) cycle
        ! Compute the enthalpy density corresponding to the current temperature
        if ( this%temps(n) <= 0 ) then
          ! Negative temperature from the diffusion solve---ignore for energy
          ! deposition purposes.
          ! Ignore this cell in forming the source term for the diffusion solver
          this%ignore_cell(n) = .true.
        else
          call this%compute_enthalpy_density(n, H)
          this%enthalpies(n) = H
        end if

      end do



    end subroutine set_state



    subroutine compute_enthalpy_density(this, n, H)

      use material_interop,     only: ds_enthalpy_density

      class(am_state), intent(inout) :: this
      integer, intent(in) :: n
      real(r8), intent(out) :: H

      integer :: m, nmat
      real(r8) :: state(1)

      nmat = this%nmat
      state(1) = this%temps(n)
      H = 0
      do m = 1, nmat
        !! note: need to multiply by the cell volume to get the total cell energy
        H = H + this%vol_fracs(m,n) * ds_enthalpy_density(m,state) 
      end do
    

    end subroutine compute_enthalpy_density



    subroutine compute_temp_from_enthalpy_density(this, n, H, T)

      class(am_state) :: this
      integer, intent(in) :: n ! cell number
      real(r8), intent(in) :: H ! enthalpy density
      real(r8), intent(out) :: T ! temperature

      integer :: nmat
      real(r8) :: Tmin, Tmax
      real(r8), allocatable :: params(:)

      nmat = this%nmat
      allocate(params(1:nmat+1))

      ! Require non-negative enthalpy
      if ( H <= 0  ) then
        call TLS_fatal ('compute_temp_from_enthalpy_density: negative enthalpy')
      end if

      ! Compute new temp. that corresponds to the enthalpy density H
      params(1) = H
      params(2:nmat+1) = this%vol_fracs(1:nmat,n)

      this%root_solver%found_root = .false.
      this%root_solver%is_bracketed = .false.
      call this%root_solver%bracket_root(residual_temp, params, Tmin, Tmax)
      if ( this%root_solver%is_bracketed ) then
        call this%root_solver%compute_root(residual_temp, params, Tmin, Tmax)
      else
        print *, "compute_updated_temperature: couldn't bracket root: "
      end if

    !  error = residual_for_temp(T, params)
      if ( .not. this%root_solver%found_root  ) then
        print *, "compute_temp_from_enthalpy_density: couldn't find root"
      else
        T = this%root_solver%root
      end if

      deallocate(params)

    end subroutine compute_temp_from_enthalpy_density


    function residual_temp(T, params)

      use parameter_module,     only: nmat
      use material_interop,     only: ds_enthalpy_density

      real(r8), intent(in) :: T, params(:) 
      real(r8) :: residual_temp

      integer :: m
      real(r8) :: H
      real(r8) :: state(1), H_target, Vofm(1:nmat)

      ! extract parameters for evaluation of the residual
      H_target = params(1)
      Vofm(1:nmat) = params(2:nmat+1)

      ! compute the enthalpy density H(T) 
      state(1) = T
      H = 0
      do m = 1, nmat
        H = H + Vofm(m) * ds_enthalpy_density(m,state) 
      end do

      residual_temp = (H - H_target)/max(H_target,1.0_r8)


    end function residual_temp


    subroutine fixup_energy_in_nearly_void_cells(this, is_nearly_void)

      use gs_module,     only: EE_GATHER

      class(am_state) :: this
      logical, intent(in) :: is_nearly_void(:)

      integer :: n, ncells, f, nfc, k, status
      real(r8) :: one_tol, am_tol, H, T
      real(r8), allocatable :: temps_nghbrs(:,:), void_vol_fracs_nghbrs(:,:)

      nfc = this%nfaces
      one_tol = this%one_tol
      am_tol = this%am_tol
      ncells = this%ncells

      allocate( temps_nghbrs(nfc,ncells), &
                void_vol_fracs_nghbrs(nfc,ncells), &
                STAT = status )
      if (status /= 0) call TLS_panic ('fixup_energy_in_nearly_void_cells: cannot allocate arrays') 

      call EE_GATHER(temps_nghbrs, this%temps)
      call EE_GATHER(void_vol_fracs_nghbrs, this%void_vol_fracs)

      do n = 1, ncells
        if ( this%void_vol_fracs(n) >  (1 - this%one_tol) ) cycle
        if ( .not. is_nearly_void(n) ) cycle
        if ( this%ignore_cell(n) ) cycle
        ! average over faces with nonzero temps
        T = 0
        k = 0
        do f = 1, nfc
          if ( temps_nghbrs(f,n) > 0  ) then
            T = T + temps_nghbrs(f,n)
            k = k + 1
          end if
        end do     
        if ( k > 0 ) then
          this%temps(n) = T/k
        else
          ! Can't fixup energy since all neighboring cells have zero or negative temperature;
          ! thus, ignore this cell as far as forming a source dQ for the diffusion solver
       !   this%ignore_cell(n) = .true.
          !!!! TODO: NEED TO CHANGE HOW THIS IS HANDLED !!!!
          print *, 'fixup_energy_in_nearly_void_cells: cannot fix cell energy'
          this%temps(n) = 300
       !   this%temps(n) = this%temps_old(n)
          !!!! TODO: !!!!
        !  call TLS_fatal ('fixup_energy_in_nearly_void_cells: cannot fix cell energy')
        end if
        ! Compute the associated enthalpy density that is consistent with this temperature
        call this%compute_enthalpy_density(n, H)
        this%enthalpies(n) = H
      end do

      deallocate(temps_nghbrs, void_vol_fracs_nghbrs)

    end subroutine fixup_energy_in_nearly_void_cells





end module am_state_type

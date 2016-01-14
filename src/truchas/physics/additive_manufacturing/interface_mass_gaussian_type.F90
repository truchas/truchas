!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module interface_mass_gaussian_type

  use kinds, only: r8
  use parameter_list_type
  use interface_mass_class
  use am_coord_system_class
  use am_state_type
  implicit none

  private

  type, extends(interface_mass), public :: interface_mass_gaussian
    real(r8) :: rate, sigma, vof_cutoff
    real(r8) :: zero_tol
    logical, allocatable :: interface_cells(:), void_interface_cells(:)
    real(r8), allocatable :: powder_absorp(:)
    contains
      ! These procedures implement the abstract interface
      ! of the base class interface_mass
      procedure :: init => init_mass_depos
      procedure :: alloc => alloc_mass_depos
      procedure :: deposit_mass => deposit_mass_gaussian
      ! These procedures are specific to the derived class
      procedure :: redistribute_overflow
      procedure :: mass_source_gaussian
      procedure :: find_and_init_interface_cells
      procedure :: update_density
  end type interface_mass_gaussian


  contains


    subroutine alloc_mass_depos(this)

      use parameter_module, only: ncells

      class(interface_mass_gaussian), intent(inout) :: this
 
      integer :: status

      ! Allocate array to hold fraction of heat to be deposited in cell

      allocate( this%interface_cells(ncells), &
                this%void_interface_cells(ncells), &
                this%powder_absorp(ncells), &
                STAT = status )

      if (status /= 0) call TLS_fatal ('alloc_mass_gaussian: allocation failed')      

    end subroutine alloc_mass_depos



    subroutine init_mass_depos(this, params)

      use parameter_list_type

      class(interface_mass_gaussian), intent(inout) :: this
      type(parameter_list), intent(inout) :: params

      call params%get('powder_rate', this%rate)
      call params%get('powder_sigma', this%sigma)  
  !    call params%get('vof_cutoff', this%vof_cutoff)  

      !!!! CHANGE !!!!
      this%zero_tol = 10.0_r8**(-8)
      this%vof_cutoff = 10.0_r8**(-12)
      !!!! CHANGE !!!!

    end subroutine init_mass_depos

   

    subroutine deposit_mass_gaussian(this, state, am_geom, dt)

      use mesh_module,          only: Cell

      class(interface_mass_gaussian), intent(inout) :: this
      type(am_state), intent(inout) :: state
      class(am_coord_system), intent(in) :: am_geom
      real(r8), intent(in) :: dt

      integer :: m, n, nfluid, m_void, m_fluid
      real(r8) :: dx, dy, delta_vof

      nfluid = size(state%fluid_indices)
      m_void = state%void_index

      ! Find all interface cells located in the mushy zone. 
      ! Note: for interface cells that are pure void, the temperature
      ! and velocity of the cell are set to that of the cell below; 
      ! in particular, the state object, state, is modified from this 
      ! function call
      call this%find_and_init_interface_cells(state)

      do n = 1, state%ncells
       if ( .not. this%interface_cells(n) ) cycle
        ! Compute the powder volume added to the cell over a time step dt
        dx = Cell(n)%Centroid(1) - am_geom%am_coords(1)
        dy = Cell(n)%Centroid(2) - am_geom%am_coords(2)
        delta_vof = this%mass_source_gaussian(dx, dy, dt)
        ! Convert this to a volume fraction, and multiply by fraction
        ! of powder that sticks to the cell
        delta_vof = this%powder_absorp(n) * delta_vof/Cell(n)%Volume
        ! Add volume fraction to fluid materials by splitting into equal amounts
        do m = 1, nfluid
          m_fluid = state%fluid_indices(m)
          state%vol_fracs(m_fluid,n) = state%vol_fracs(m_fluid,n) + delta_vof/nfluid
        end do 
        ! Subtract volume fraction from the void (this can go negative here---if it does, 
        ! it will be fixed in  the fluid redistribution phase)
        state%vol_fracs(m_void,n) = state%vol_fracs(m_void,n) - delta_vof
        state%void_vol_fracs(n) = state%vol_fracs(m_void,n)
      end do

      ! redistribute volume fractions in cels that have "overflowed"
      call this%redistribute_overflow(state)
      ! change the densities to be consistent with the new volume fractions
      call this%update_density(state)

    end subroutine deposit_mass_gaussian



    subroutine redistribute_overflow(this, state)

      use gs_module,     only: EE_GATHER
      use mesh_module,      only: Cell

      class(interface_mass_gaussian), intent(inout) :: this
      type(am_state), intent(inout) :: state

      integer :: m_fluid, n, m, nfc, ncells, f, f_up(1), f_down(1), status, m_void
      real(r8) :: volume_void, volume_void_nghbr, upward_vec(3), downward_vec(3)
      real(r8), allocatable :: temps_nghbrs(:,:), velocities_x_nghbrs(:,:), &
                               velocities_y_nghbrs(:,:), velocities_z_nghbrs(:,:), &
                               void_vol_fracs_nghbrs(:,:), proj_onto_down_vec(:)

      nfc = state%nfaces
      ncells = state%ncells
      m_void = state%void_index
      downward_vec = [0,0,-1]

      allocate( void_vol_fracs_nghbrs(nfc,ncells), &
                temps_nghbrs(nfc,ncells), &
                velocities_x_nghbrs(nfc,ncells), &
                velocities_y_nghbrs(nfc,ncells), &
                velocities_z_nghbrs(nfc,ncells), &
                proj_onto_down_vec(nfc), &
                STAT = status )
      if (status /= 0) call TLS_panic ('redistribute_overflow: allocation failed')
                               
      ! Gather neighboring cell-centered velocities, temperatures, and
      ! enthalpies in order to update initially void that cells that have
      ! get filled with incoming fluid from neighboring, "overflowing" cells
      call EE_GATHER(void_vol_fracs_nghbrs, state%void_vol_fracs)
      call EE_GATHER(temps_nghbrs, state%temps)
      call EE_GATHER(velocities_x_nghbrs, state%velocities_x)
      call EE_GATHER(velocities_y_nghbrs, state%velocities_y)
      call EE_GATHER(velocities_z_nghbrs, state%velocities_z)

      ! Accept all overflowing fluid from the cell below
      do n = 1, state%ncells
        volume_void = state%void_vol_fracs(n)
        ! Find the cell below 
        do f = 1, nfc
          proj_onto_down_vec(f) = &
                   sum( Cell(n)%Face_Normal(1:3,f) * downward_vec(1:3) )
        end do
        f_down(1:1) = maxloc(proj_onto_down_vec(1:nfc)) 
        ! If the cell below has a negative void volume fraction, 
        ! then get the oveflowing fluid
        volume_void_nghbr = void_vol_fracs_nghbrs(f_down(1),n)
        if ( volume_void_nghbr < 0 ) then
          do m = 1, state%nfluids
            m_fluid = state%fluid_indices(m)
            ! Split the amount of overflowing fluid in the cell below 
            state%vol_fracs(m_fluid,n) = state%vol_fracs(m_fluid,n) + abs(volume_void_nghbr)/state%nfluids
          end do
          state%vol_fracs(m_void,n) = state%vol_fracs(m_void,n) - abs(volume_void_nghbr)
          state%void_vol_fracs(n) = state%vol_fracs(m_void,n)
          ! If the accepting cell was pure void, set the temperatures and velocities
          ! to be those of the donating cell below
          if ( volume_void > (1-state%one_tol) ) then
            state%temps(n) = temps_nghbrs(f_down(1),n)
            state%velocities_x(n) = velocities_x_nghbrs(f_down(1),n)
            state%velocities_y(n) = velocities_y_nghbrs(f_down(1),n)
            state%velocities_z(n) = velocities_z_nghbrs(f_down(1),n) 
          end if
        end if
      end do

      ! Donate all overflowing fluid to cell the above
      do n = 1, state%ncells
        volume_void = state%void_vol_fracs(n)
        if ( volume_void < 0 ) then
          ! Split the amount of outgoing fluid equally among the (non-void) fluids
          do m = 1, state%nfluids
            m_fluid = state%fluid_indices(m)
            ! Subtract off the amount of fluid material m_fluid leaving the cell
            state%vol_fracs(m_fluid,n) = state%vol_fracs(m_fluid,n) - abs(volume_void)/state%nfluids
          end do
          state%vol_fracs(m_void,n) = 0
          state%void_vol_fracs(n) = 0
        end if
      end do

      deallocate( void_vol_fracs_nghbrs, temps_nghbrs, &
                  velocities_x_nghbrs, velocities_y_nghbrs, &
                  velocities_z_nghbrs, proj_onto_down_vec  )

    end subroutine redistribute_overflow



    function mass_source_gaussian(this, dx, dy, dt) result(q)

      class(interface_mass_gaussian), intent(in) :: this
      real(r8), intent(in) :: dx, dy, dt
      real(r8) :: q

      real(r8) :: r, two_sigma2
      real(r8), parameter :: PI = 3.141592653589793_r8

      two_sigma2 = this%sigma**2
      r = 2 * (dx**2 + dy**2) / two_sigma2
   !   q = dt * this%rate * exp(-r) / (PI*two_sigma2)
      q = dt * this%rate * exp(-r)


    end function mass_source_gaussian



    subroutine find_and_init_interface_cells(this, state)

      use gs_module,     only: EE_GATHER

      use parameter_module, only: ndim
      use mesh_module,      only: Cell

      class(interface_mass_gaussian), intent(inout) :: this
      type(am_state), intent(inout) :: state

      integer :: m, n, f, f_down(1), ncells, nfc, &
                 n_fluid, m_fluid, m_void, nmat, status
      real(r8) :: downward_unit_vec(1:3), volume_void_neighbor, &
                  vof_cutoff, zero_tol, &
                  volume_fluid_neighbor, fluid_vol_tot
      logical :: is_pure_void_over_mushy
      real(r8), allocatable :: normals_dotted_with_down_vec(:), &
                               temps_nghbrs(:,:), &
                               velocities_x_nghbrs(:,:), &
                               velocities_y_nghbrs(:,:), &
                               velocities_z_nghbrs(:,:), &
                               vol_fracs_nghbrs(:,:,:)


      nfc = state%nfaces
      ncells = state%ncells
      nmat = state%nmat
      n_fluid = size(state%fluid_indices,1)
      m_void = state%void_index
      vof_cutoff = this%vof_cutoff
      zero_tol = this%zero_tol

      allocate( vol_fracs_nghbrs(1:nmat,1:nfc,1:ncells), &
                temps_nghbrs(1:nfc,1:ncells), &
                velocities_x_nghbrs(1:nfc,1:ncells), &
                velocities_y_nghbrs(1:nfc,1:ncells), &
                velocities_z_nghbrs(1:nfc,1:ncells), &
                normals_dotted_with_down_vec(1:nfc), &
                STAT = status )
      if (status /= 0) call TLS_panic ('fixup_energy_in_cells_with_tiny_fluid_volumes') 

      do m = 1,nmat
        call EE_GATHER(vol_fracs_nghbrs(m,:,:), state%vol_fracs(m,:))
      end do
      call EE_GATHER(temps_nghbrs, state%temps)
      call EE_GATHER(velocities_x_nghbrs, state%velocities_x)
      call EE_GATHER(velocities_y_nghbrs, state%velocities_y)
      call EE_GATHER(velocities_z_nghbrs, state%velocities_z)
    
      downward_unit_vec = (/0,0,-1/) 
      this%interface_cells = .false.
      this%powder_absorp = 0
      do n = 1, ncells

        ! Find all pure void cells that are directly on top of mushy cells without any void
        if ( state%void_vol_fracs(n) > (1 - zero_tol) ) then  ! is a pure void cell

          ! find the face number of the cell directly under the cell n
          ! (the normal vector dotted with the downward direction is largest)
          ! NOTE: THIS ASSUMES THAT THE CELLS ARE HEXES
          do f = 1, nfc
            normals_dotted_with_down_vec(f) = &
                   sum( Cell(n)%Face_Normal(1:ndim,f) * downward_unit_vec(1:ndim) )
          end do
          f_down(1:1) = maxloc( normals_dotted_with_down_vec(1:nfc) )
          ! Find the void and (non-void) fluid volumes of the cell below
          volume_void_neighbor = vol_fracs_nghbrs(m_void,f_down(1),n)
          volume_fluid_neighbor = 0
          do m = 1, n_fluid
            m_fluid = state%fluid_indices(m)
            volume_fluid_neighbor = volume_fluid_neighbor + &
                                      vol_fracs_nghbrs(m_fluid,f_down(1),n)
          end do 
          ! If the cell is pure void and the fluid volume in the cell below is larger than
          ! the cutoff (i.e., this is a mushy cell), then this is an interface cell
          this%interface_cells(n) = &
              (volume_void_neighbor < vof_cutoff) .and. (volume_fluid_neighbor > vof_cutoff)
          is_pure_void_over_mushy = this%interface_cells(n)
          ! The fraction of powder particles that stick to the cell is equal to the volume
          ! fraction of fluid in the cell below
          this%powder_absorp(n) = volume_fluid_neighbor
          ! Set the velocities and temperatures of the void cell to that of the cell below if this is a pure
          ! void cell over a mushy cell
          if ( is_pure_void_over_mushy ) then
            state%velocities_x(n) = velocities_x_nghbrs(f_down(1),n)
            state%velocities_y(n) = velocities_y_nghbrs(f_down(1),n)
            state%velocities_z(n) = velocities_z_nghbrs(f_down(1),n)
            state%temps(n) = temps_nghbrs(f_down(1),n)
          end if

        end if


        if ( state%void_vol_fracs(n) > vof_cutoff .and. &
             state%void_vol_fracs(n) < 1-vof_cutoff ) then
          ! Check if the cell is in the mushy zone (otherwise powder particles bounce off)
          fluid_vol_tot = 0
          ! Find the total (non-void) fluid volume
          do m = 1, n_fluid
            m_fluid = state%fluid_indices(m)
            fluid_vol_tot = fluid_vol_tot + state%vol_fracs(m_fluid,n)
          end do 
          ! If the fluid volume is larger than the cutoff, this is a mushy cell
          this%interface_cells(n) = ( fluid_vol_tot > vof_cutoff )
          ! The fraction of powder particles that stick to the cell is equal to the volume
          ! fraction of fluid in the cell
          this%powder_absorp(n) = fluid_vol_tot
        end if

      end do


      deallocate( vol_fracs_nghbrs, temps_nghbrs, &
                  velocities_x_nghbrs, velocities_y_nghbrs, &
                  velocities_z_nghbrs, normals_dotted_with_down_vec )


    end subroutine find_and_init_interface_cells



    subroutine update_density(this, state)

      use property_module, only: DENSITY_MATERIAL

      class(interface_mass_gaussian), intent(inout) :: this
      type(am_state), intent(inout) :: state

      integer :: m
      real(r8) :: m_density

      state%densities = 0.0_r8
      do m = 1, state%nmat
         m_density = DENSITY_MATERIAL(m)
         where (state%vol_fracs(m,:) > 0.0_r8) &
              state%densities = state%densities + state%vol_fracs(m,:)*m_density
      end do



    end subroutine update_density




    subroutine print_obj(this)

      class(interface_mass_gaussian), intent(inout) :: this

      print *, "powder rate: ", this%rate
      print *, "powder sigma: ", this%sigma

    end subroutine print_obj





end module interface_mass_gaussian_type


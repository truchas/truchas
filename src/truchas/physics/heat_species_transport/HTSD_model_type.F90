!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module HTSD_model_type

  use kinds
  use unstr_mesh_type
  use mfd_disc_type
  use data_layout_type
  use property_mesh_function
  use source_mesh_function
  use boundary_data
  use bndry_func_class
  use interface_data
  use rad_problem_type
  use index_partitioning
  use timing_tree
  implicit none
  private
  
  type, public :: HT_model
    !! Equation parameters
    type(prop_mf) :: conductivity ! thermal conductivity
    type(prop_mf) :: H_of_T       ! enthalpy as a function of temperature
    type(source_mf) :: source     ! external heat source
    !! Boundary condition data
    type(bd_data) :: bc_dir  ! Dirichlet
    type(bd_data) :: bc_flux ! simple flux
    type(bd_data) :: bc_htc  ! external HTC (coef, ref temp)
    type(bd_data) :: bc_rad  ! simple radiation (eps, amb temp)
    type(if_data) :: ic_htc  ! internal HTC
    type(if_data) :: ic_rad  ! internal gap radiation
    class(bndry_func), allocatable :: evap_flux
    real(r8) :: sbconst, abszero ! Stefan-Boltzmann constant and absolute zero for radiation BC
    !! Enclosure radiation problems
    type(rad_problem), pointer :: vf_rad_prob(:) => null()
    real(r8), allocatable :: face_flux(:), face_temp(:)
  end type HT_model
  
  type, public :: SD_model
    !! Equation parameters
    type(prop_mf) :: diffusivity
    type(source_mf) :: source
    type(prop_mf), pointer :: soret => null()
    !! Boundary condition data
    type(bd_data) :: bc_dir   ! Dirichlet
    type(bd_data) :: bc_flux  ! simple flux
  end type SD_model
  
  type, public :: HTSD_model
    integer :: num_comp = 0
    type(HT_model), pointer :: ht => null()
    type(SD_model), pointer :: sd(:) => null()
    type(mfd_disc),  pointer :: disc => null()
    type(unstr_mesh), pointer :: mesh => null()
    logical, pointer :: void_cell(:) => null(), void_face(:) => null()
    real(r8) :: void_temp = 0.0_r8
    type(data_layout) :: layout
    integer :: cell_heat_segid, cell_temp_segid, face_temp_segid
    integer, pointer :: rad_segid(:) => null()
    integer, pointer :: cell_conc_segid(:) => null(), face_conc_segid(:) => null()
  end type HTSD_model
  
  public :: HTSD_model_init
  public :: HTSD_model_delete
  public :: HTSD_model_size
  public :: HTSD_model_compute_f
  public :: HTSD_model_new_state_array
  
  public :: HTSD_model_get_cell_heat_view, HTSD_model_get_cell_heat_copy
  public :: HTSD_model_get_cell_temp_view, HTSD_model_get_cell_temp_copy
  public :: HTSD_model_get_face_temp_view, HTSD_model_get_face_temp_copy
  public :: HTSD_model_get_radiosity_view, HTSD_model_get_radiosity_copy
  public :: HTSD_model_get_cell_conc_view, HTSD_model_get_cell_conc_copy
  public :: HTSD_model_get_face_conc_view, HTSD_model_get_face_conc_copy
  
  public :: HTSD_model_set_cell_heat
  public :: HTSD_model_set_cell_temp
  public :: HTSD_model_set_face_temp
  public :: HTSD_model_set_radiosity
  public :: HTSD_model_set_cell_conc
  public :: HTSD_model_set_face_conc
  
contains
  
  subroutine HTSD_model_init (this, disc, htmodel, sdmodel)
  
    type(HTSD_model), intent(out) :: this
    type(mfd_disc), intent(in), target :: disc
    type(HT_model), pointer :: htmodel
    type(SD_model), pointer :: sdmodel(:)
    
    integer :: n
    
    this%disc => disc
    this%mesh => disc%mesh
    
    this%ht => htmodel  ! take ownership
    this%sd => sdmodel  ! take ownership
    
    ASSERT(associated(this%ht) .or. associated(this%sd))
    
    !! Create the packed layout of the model variables.
    if (associated(this%ht)) then
      this%cell_heat_segid = alloc_segment(this%layout, this%mesh%ncell_onP)
      this%cell_temp_segid = alloc_segment(this%layout, this%mesh%ncell_onP)
      this%face_temp_segid = alloc_segment(this%layout, this%mesh%nface_onP)
      if (associated(this%ht%vf_rad_prob)) then
        ASSERT(size(this%ht%vf_rad_prob) > 0)
        allocate(this%rad_segid(size(this%ht%vf_rad_prob)))
        do n = 1, size(this%ht%vf_rad_prob)
          this%rad_segid(n) = alloc_segment(this%layout, size(this%ht%vf_rad_prob(n)%faces))
        end do
      end if
      allocate (this%ht%face_flux(this%mesh%nface))
      allocate (this%ht%face_temp(this%mesh%nface))      
    end if
    if (associated(this%sd)) then
      ASSERT(size(this%sd) > 0)
      this%num_comp = size(this%sd)
      allocate(this%cell_conc_segid(this%num_comp), this%face_conc_segid(this%num_comp))
      do n = 1, this%num_comp
        this%cell_conc_segid(n) = alloc_segment(this%layout, this%mesh%ncell_onP)
        this%face_conc_segid(n) = alloc_segment(this%layout, this%mesh%nface_onP)
      end do
    end if
    call alloc_complete (this%layout)
    
  end subroutine HTSD_model_init
  
  subroutine HTSD_model_delete (this)
    type(HTSD_model), intent(inout) :: this
    integer :: n
    call delete_layout (this%layout)
    if (associated(this%rad_segid)) deallocate(this%rad_segid)
    if (associated(this%cell_conc_segid)) deallocate(this%cell_conc_segid)
    if (associated(this%face_conc_segid)) deallocate(this%face_conc_segid)
    if (associated(this%ht)) then
      call HT_model_delete (this%ht)
      deallocate(this%ht)
    end if
    if (associated(this%sd)) then
      do n = 1, size(this%sd)
        call SD_model_delete (this%sd(n))
      end do
      deallocate(this%sd)
    end if
  end subroutine HTSD_model_delete
  
  subroutine HT_model_delete (this)
    type(HT_model), intent(inout) :: this
    integer :: n
    call destroy (this%conductivity)
    call destroy (this%H_of_T)
    call smf_destroy (this%source)
    call bd_data_destroy (this%bc_dir)
    call bd_data_destroy (this%bc_flux)
    call bd_data_destroy (this%bc_htc)
    call bd_data_destroy (this%bc_rad)
    call if_data_destroy (this%ic_htc)
    call if_data_destroy (this%ic_rad)
    if (associated(this%vf_rad_prob)) deallocate(this%vf_rad_prob)
    if (allocated(this%face_flux)) deallocate(this%face_flux)
    if (allocated(this%face_temp)) deallocate(this%face_temp)
  end subroutine HT_model_delete
  
  subroutine SD_model_delete (this)
    type(SD_model), intent(inout) :: this
    call destroy (this%diffusivity)
    if (associated(this%soret)) then
      call destroy (this%soret)
      deallocate(this%soret)
    end if
    call smf_destroy (this%source)
    call bd_data_destroy (this%bc_dir)
    call bd_data_destroy (this%bc_flux)
  end subroutine SD_model_delete
  
  function HTSD_model_new_state_array (this, u) result (state)
    type(HTSD_model), intent(in) :: this
    real(r8), intent(in) :: u(:)
    real(r8), pointer :: state(:,:)
    integer :: n
    n = 1
    if (associated(this%ht)) n = 0
    allocate(state(this%mesh%ncell,n:this%num_comp))
    if (associated(this%ht)) then
      call HTSD_model_get_cell_temp_copy (this, u, state(:,0))
      call gather_boundary (this%mesh%cell_ip, state(:,0))
    end if
    if (associated(this%sd)) then
      do n = 1, this%num_comp
        call HTSD_model_get_cell_conc_copy (this, n, u, state(:,n))
        call gather_boundary (this%mesh%cell_ip, state(:,n))
      end do
    end if
  end function HTSD_model_new_state_array

  subroutine HTSD_model_compute_f (this, t, u, udot, f)
  
    use mfd_disc_type
  
    type(HTSD_model), intent(inout) :: this
    real(r8), intent(in)  :: t, u(:), udot(:)
    real(r8), intent(out) :: f(:)
    target :: u, udot, f
    
    integer :: n
    !real(r8), allocatable, target :: state(:,:)
    real(r8), pointer :: state(:,:)
    
    call start_timer ('HTSD function')
    
    state => HTSD_model_new_state_array(this, u)
    
    !! HT residual.
    if (associated(this%ht)) then
      call start_timer ('HT function')
      call HT_model_compute_f
      call stop_timer ('HT function')
    end if
    
    !! SD residual.
    if (associated(this%sd)) then
      call start_timer ('SD function')
      do n = 1, this%num_comp
        call SD_model_compute_f (n)
      end do
      call stop_timer ('SD function')
    end if
    
    deallocate(state)
    
    call stop_timer ('HTSD function')

  contains
  
    subroutine HT_model_compute_f ()

      integer :: j, n, n1, n2
      real(r8) :: term
      real(r8) :: bcflux
      real(r8), pointer :: uptr(:), fptr(:), qrad(:)
      real(r8), dimension(this%mesh%ncell) :: Fcell, Tcell, Hdot, value
      real(r8), dimension(this%mesh%nface) :: Fface, Tface
      real(r8), allocatable :: Fdir(:), flux(:)
      integer, pointer :: faces(:)
      logical, allocatable :: void_link(:)

    !!!! RESIDUAL OF THE ALGEBRAIC ENTHALPY-TEMPERATURE RELATION !!!!!!!!!!!!!!!!!

      call HTSD_model_get_cell_heat_view (this, u, uptr)
      call HTSD_model_get_cell_heat_view (this, f, fptr)
      call pmf_eval (this%ht%H_of_T, state, value)
      fptr = uptr - value(:this%mesh%ncell_onP)

      !! Overwrite function value on void cells with dummy equation H=0.
      if (associated(this%void_cell)) where (this%void_cell(:this%mesh%ncell_onP)) fptr = uptr

    !!!! RESIDUAL OF THE HEAT EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !! Off-process-extended cell and face temperatures.
      call HTSD_model_get_cell_temp_copy (this, u, Tcell)
      call gather_boundary (this%mesh%cell_ip, Tcell)
      call HTSD_model_get_face_temp_copy (this, u, Tface)
      call gather_boundary (this%mesh%face_ip, Tface)

      !! Off-process-extended cell enthalpy time derivative.
      call HTSD_model_get_cell_heat_copy (this, udot, Hdot)
      call gather_boundary (this%mesh%cell_ip, Hdot)

      !! Pre-compute the Dirichlet condition residual and
      !! impose the Dirichlet data on the face temperature.
      call bd_data_eval (this%ht%bc_dir, t)
      allocate(Fdir(size(this%ht%bc_dir%faces)))
      do j = 1, size(this%ht%bc_dir%faces)
        Fdir(j) = Tface(this%ht%bc_dir%faces(j)) - this%ht%bc_dir%values(1,j)
        Tface(this%ht%bc_dir%faces(j)) = this%ht%bc_dir%values(1,j)
      end do

      !! Compute the generic heat equation residual.
      call pmf_eval (this%ht%conductivity, state, value)
      if (associated(this%void_cell)) where (this%void_cell) value = 0.0_r8
      call this%disc%apply_diff (value, Tcell, Tface, Fcell, Fface,this%ht%face_flux)
      call smf_eval (this%ht%source, t, value)
      Fcell = Fcell + this%mesh%volume*(Hdot - value)

      !! Dirichlet condition residuals.
      do j = 1, size(this%ht%bc_dir%faces)
        n = this%ht%bc_dir%faces(j)
        Fface(n) = Fdir(j)  !! overwrite with pre-computed values
      end do
      deallocate(Fdir)

      !! Simple flux BC contribution.
      call bd_data_eval (this%ht%bc_flux, t)
      do j = 1, size(this%ht%bc_flux%faces)
        n = this%ht%bc_flux%faces(j)
        bcflux = this%ht%bc_flux%values(1,j) 
        Fface(n) = Fface(n) + this%mesh%area(n) * bcflux
        this%ht%face_flux(n) = bcflux
      end do

      !! External HTC flux contribution.
      call bd_data_eval (this%ht%bc_htc, t)
      do j = 1, size(this%ht%bc_htc%faces)
        n = this%ht%bc_htc%faces(j)
        bcflux = this%ht%bc_htc%values(1,j) * &
                                  (Tface(n) - this%ht%bc_htc%values(2,j))
        Fface(n) = Fface(n) + this%mesh%area(n) * bcflux
!        this%ht%face_flux(n) = this%ht%face_flux(n) + bcflux      
      end do

      !! Internal HTC flux contribution.
      call if_data_eval (this%ht%ic_htc, t)
      allocate(void_link(size(this%ht%ic_htc%faces,2)))
      if (associated(this%void_face)) then
        do j = 1, size(void_link)
          void_link(j) = any(this%void_face(this%ht%ic_htc%faces(:,j)))
        end do
      else
        void_link = .false.
      end if
      do j = 1, size(this%ht%ic_htc%faces,2)
        if (void_link(j)) cycle
        n1 = this%ht%ic_htc%faces(1,j)
        n2 = this%ht%ic_htc%faces(2,j)
        bcflux = this%ht%ic_htc%values(1,j) * (Tface(n1) - Tface(n2))
        term = this%mesh%area(n1) * bcflux
        Fface(n1) = Fface(n1) + term
        Fface(n2) = Fface(n2) - term
!        this%ht%face_flux(n1) = this%ht%face_flux(n1) + bcflux
!        this%ht%face_flux(n2) = this%ht%face_flux(n2) - bcflux        
      end do
      if (allocated(void_link)) deallocate(void_link)
      
      !! Internal gap radiation contribution.
      call if_data_eval (this%ht%ic_rad, t)
      allocate(void_link(size(this%ht%ic_rad%faces,2)))
      if (associated(this%void_face)) then
        do j = 1, size(void_link)
          void_link(j) = any(this%void_face(this%ht%ic_rad%faces(:,j)))
        end do
      else
        void_link = .false.
      end if
      do j = 1, size(this%ht%ic_rad%faces,2)
        if (void_link(j)) cycle
        n1 = this%ht%ic_rad%faces(1,j)
        n2 = this%ht%ic_rad%faces(2,j)
        bcflux = this%ht%ic_rad%values(1,j) * this%ht%sbconst * &
               ((Tface(n1)-this%ht%abszero)**4 - (Tface(n2)-this%ht%abszero)**4)
        term = this%mesh%area(n1) * bcflux
        Fface(n1) = Fface(n1) + term
        Fface(n2) = Fface(n2) - term
!        this%ht%face_flux(n1) = this%ht%face_flux(n1) + bcflux
!        this%ht%face_flux(n2) = this%ht%face_flux(n2) - bcflux        
      end do
      if (allocated(void_link)) deallocate(void_link)

      !! Ambient radiation BC flux contribution.
      call bd_data_eval (this%ht%bc_rad, t)
      do j = 1, size(this%ht%bc_rad%faces)
        n = this%ht%bc_rad%faces(j)
        bcflux = this%ht%bc_rad%values(1,j) * &
                                this%ht%sbconst * ((Tface(n) - this%ht%abszero)**4 - &
                                           (this%ht%bc_rad%values(2,j)-this%ht%abszero)**4)
        Fface(n) = Fface(n) + this%mesh%area(n) * bcflux
        this%ht%face_flux(n) = this%ht%face_flux(n) + bcflux
      end do

      !! Experimental evaporation heat flux
      if (allocated(this%ht%evap_flux)) then
        call this%ht%evap_flux%compute_value(t, Tface)
        do j = 1, size(this%ht%evap_flux%index)
          n = this%ht%evap_flux%index(j)
          bcflux = this%ht%evap_flux%value(j)
          Fface(n) = Fface(n) + this%mesh%area(n)*bcflux
          this%ht%face_flux(n) = this%ht%face_flux(n) + bcflux
        end do
      end if

      !! Overwrite function value on void cells and faces with dummy equation T=0.
      if (associated(this%void_cell)) where (this%void_cell) Fcell = Tcell - this%void_temp
      if (associated(this%void_face)) where (this%void_face) Fface = Tface - this%void_temp

    !!!! RESIDUALS OF THE ENCLOSURE RADIATION SYSTEMS !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (associated(this%ht%vf_rad_prob)) then
        do n = 1, size(this%ht%vf_rad_prob)
          call HTSD_model_get_radiosity_view (this, n, u, qrad)
          faces => this%ht%vf_rad_prob(n)%faces
          !! Radiative heat flux contribution to the heat conduction face residual.
          allocate(flux(size(faces)))
          call this%ht%vf_rad_prob(n)%heat_flux (t, qrad, Tface(faces), flux)
          do j = 1, size(faces)
            Fface(faces(j)) = Fface(faces(j)) + this%mesh%area(faces(j)) * flux(j)
          end do
          deallocate(flux)
          !! Residual of the algebraic radiosity system.
          call HTSD_model_get_radiosity_view (this, n, f, fptr)
          call this%ht%vf_rad_prob(n)%residual (t, qrad, Tface(faces), fptr)
          fptr = -fptr
        end do
      end if

      !! Return the on-process part of the heat conduction residuals (the rest is junk)
      call HTSD_model_set_cell_temp (this, Fcell, f)
      call HTSD_model_set_face_temp (this, Fface, f)

      this%ht%face_temp(:) = Tface(:)

    end subroutine HT_model_compute_f

    subroutine SD_model_compute_f (index)
    
      integer, intent(in) :: index

      integer :: j, n
      real(r8), dimension(this%mesh%ncell) :: Fcell, Ccell, Cdot, D, value
      real(r8), dimension(this%mesh%nface) :: Fface, Cface
      real(r8), allocatable :: Fdir(:)
      real(r8), pointer :: fview(:), Tcell(:)
      real(r8), allocatable :: Tface(:)

      !! Off-process extended cell and face concentrations.
      call HTSD_model_get_cell_conc_copy (this, index, u, Ccell)
      call gather_boundary (this%mesh%cell_ip, Ccell)
      call HTSD_model_get_face_conc_copy (this, index, u, Cface)
      call gather_boundary (this%mesh%face_ip, Cface)

      !! Pre-compute the Dirichlet condition constraint and
      !! impose the Dirichlet data on the face concentrations.
      call bd_data_eval (this%sd(index)%bc_dir, t)
      allocate(Fdir(size(this%sd(index)%bc_dir%faces)))
      do j = 1, size(this%sd(index)%bc_dir%faces)
        Fdir(j) = Cface(this%sd(index)%bc_dir%faces(j)) - this%sd(index)%bc_dir%values(1,j)
        Cface(this%sd(index)%bc_dir%faces(j)) = this%sd(index)%bc_dir%values(1,j)
      end do

      !! Diffusion operator function value.
      call pmf_eval (this%sd(index)%diffusivity, state, D)
      if (associated(this%void_cell)) where (this%void_cell) D = 0.0_r8
      call this%disc%apply_diff (D, Ccell, Cface, Fcell, Fface)

      !! Time derivative and source contribution.
      call HTSD_model_get_cell_conc_copy (this, index, udot, Cdot)
      call gather_boundary (this%mesh%cell_ip, Cdot)
      call smf_eval (this%sd(index)%source, t, value)
      Fcell = Fcell + this%mesh%volume*(Cdot - value)

      !! Dirichlet condition residuals.
      do j = 1, size(this%sd(index)%bc_dir%faces)
        Fface(this%sd(index)%bc_dir%faces(j)) = Fdir(j)  ! overwrite with the pre-computed values
      end do
      deallocate(Fdir)

      !! Simple flux BC contribution.
      call bd_data_eval (this%sd(index)%bc_flux, t)
      do j = 1, size(this%sd(index)%bc_flux%faces)
        n = this%sd(index)%bc_flux%faces(j)
        Fface(n) = Fface(n) + this%mesh%area(n) * this%sd(index)%bc_flux%values(1,j)
      end do

      !! Overwrite function value on void cells and faces with dummy equation C=0.
      if (associated(this%void_cell)) where (this%void_cell) Fcell = Ccell
      if (associated(this%void_face)) where (this%void_face) Fface = Cface

      !! Return the on-process part of the function values (the rest is junk).
      call HTSD_model_set_cell_conc (this, index, Fcell, f)
      call HTSD_model_set_face_conc (this, index, Fface, f)

      !! Soret coupling contribution.
      if (associated(this%sd(index)%soret)) then
        Tcell => state(:,0)
        !! Get face temperatures overwritten with Dirichlet values.
        allocate(Tface(this%mesh%nface))
        call HTSD_model_get_face_temp_copy (this, u, Tface)
        call gather_boundary (this%mesh%face_ip, Tface)
        call bd_data_eval (this%ht%bc_dir, t)
        Tface(this%ht%bc_dir%faces) = this%ht%bc_dir%values(1,:)
        call pmf_eval (this%sd(index)%soret, state, value)
        value = D*value
        call this%disc%apply_diff (value, Tcell, Tface, Fcell, Fface)
        call bd_data_eval (this%sd(index)%bc_dir, t)
        Fface(this%sd(index)%bc_dir%faces) = 0.0_r8
        !! Fcell and Fface should already be 0 at all void cells and faces as required.
        call HTSD_model_get_cell_conc_view (this, index, f, fview)
        fview = fview + Fcell(:size(fview))
        call HTSD_model_get_face_conc_view (this, index, f, fview)
        fview = fview + Fface(:size(fview))
        deallocate(Tface)
      end if

    end subroutine SD_model_compute_f

  end subroutine HTSD_model_compute_f

  
  pure integer function HTSD_model_size (this)
    type(HTSD_model), intent(in) :: this
    HTSD_model_size = layout_size(this%layout)
  end function HTSD_model_size
  
  subroutine HTSD_model_get_cell_heat_view (this, array, view)
    type(HTSD_model), intent(in) :: this
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    ASSERT(associated(this%ht))
    call get_segment_view (this%layout, array, this%cell_heat_segid, view)
  end subroutine HTSD_model_get_cell_heat_view
  
  subroutine HTSD_model_get_cell_temp_view (this, array, view)
    type(HTSD_model), intent(in) :: this
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    ASSERT(associated(this%ht))
    call get_segment_view (this%layout, array, this%cell_temp_segid, view)
  end subroutine HTSD_model_get_cell_temp_view
  
  subroutine HTSD_model_get_face_temp_view (this, array, view)
    type(HTSD_model), intent(in) :: this
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    ASSERT(associated(this%ht))
    call get_segment_view (this%layout, array, this%face_temp_segid, view)
  end subroutine HTSD_model_get_face_temp_view
  
  subroutine HTSD_model_get_radiosity_view (this, index, array, view)
    type(HTSD_model), intent(in) :: this
    integer, intent(in) :: index
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    ASSERT(associated(this%rad_segid))
    ASSERT(index > 0 .and. index <= size(this%rad_segid))
    call get_segment_view (this%layout, array, this%rad_segid(index), view)
  end subroutine HTSD_model_get_radiosity_view
  
  subroutine HTSD_model_get_cell_heat_copy (this, array, copy)
    type(HTSD_model), intent(in) :: this
    real(r8), intent(in), target :: array(:)
    real(r8), intent(inout) :: copy(:)
    ASSERT(associated(this%ht))
    call get_segment_copy (this%layout, array, this%cell_heat_segid, copy)
  end subroutine HTSD_model_get_cell_heat_copy

  subroutine HTSD_model_get_cell_temp_copy (this, array, copy)
    type(HTSD_model), intent(in) :: this
    real(r8), intent(in), target :: array(:)
    real(r8), intent(inout) :: copy(:)
    ASSERT(associated(this%ht))
    call get_segment_copy (this%layout, array, this%cell_temp_segid, copy)
  end subroutine HTSD_model_get_cell_temp_copy

  subroutine HTSD_model_get_face_temp_copy (this, array, copy)
    type(HTSD_model), intent(in) :: this
    real(r8), intent(in), target :: array(:)
    real(r8), intent(inout) :: copy(:)
    ASSERT(associated(this%ht))
    call get_segment_copy (this%layout, array, this%face_temp_segid, copy)
  end subroutine HTSD_model_get_face_temp_copy

  subroutine HTSD_model_get_radiosity_copy (this, index, array, copy)
    type(HTSD_model), intent(in) :: this
    integer, intent(in) :: index
    real(r8), intent(in), target :: array(:)
    real(r8), intent(inout) :: copy(:)
    ASSERT(associated(this%rad_segid))
    ASSERT(index > 0 .and. index <= size(this%rad_segid))
    call get_segment_copy (this%layout, array, this%rad_segid(index), copy)
  end subroutine HTSD_model_get_radiosity_copy
  
  subroutine HTSD_model_get_cell_conc_view (this, index, array, view)
    type(HTSD_model), intent(in) :: this
    integer, intent(in) :: index
    real(r8), intent(in), target :: array(:)
    real(r8), pointer :: view(:)
    ASSERT(associated(this%cell_conc_segid))
    ASSERT(index > 0 .and. index <= size(this%cell_conc_segid))
    call get_segment_view (this%layout, array, this%cell_conc_segid(index), view)
  end subroutine HTSD_model_get_cell_conc_view
  
  subroutine HTSD_model_get_face_conc_view (this, index, array, view)
    type(HTSD_model), intent(in) :: this
    integer, intent(in) :: index
    real(r8), intent(in), target :: array(:)
    real(r8), pointer :: view(:)
    ASSERT(associated(this%face_conc_segid))
    ASSERT(index > 0 .and. index <= size(this%face_conc_segid))
    call get_segment_view (this%layout, array, this%face_conc_segid(index), view)
  end subroutine HTSD_model_get_face_conc_view
  
  subroutine HTSD_model_get_cell_conc_copy (this, index, array, copy)
    type(HTSD_model), intent(in) :: this
    integer, intent(in) :: index
    real(r8), intent(in), target :: array(:)
    real(r8), intent(inout) :: copy(:)
    ASSERT(associated(this%cell_conc_segid))
    ASSERT(index > 0 .and. index <= size(this%cell_conc_segid))
    call get_segment_copy (this%layout, array, this%cell_conc_segid(index), copy)
  end subroutine HTSD_model_get_cell_conc_copy
  
  subroutine HTSD_model_get_face_conc_copy (this, index, array, copy)
    type(HTSD_model), intent(in) :: this
    integer, intent(in) :: index
    real(r8), intent(in), target :: array(:)
    real(r8), intent(inout) :: copy(:)
    ASSERT(associated(this%face_conc_segid))
    ASSERT(index > 0 .and. index <= size(this%face_conc_segid))
    call get_segment_copy (this%layout, array, this%face_conc_segid(index), copy)
  end subroutine HTSD_model_get_face_conc_copy
  
  subroutine HTSD_model_set_cell_heat (this, source, array)
    type(HTSD_model), intent(in) :: this
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout) :: array(:)
    INSIST(associated(this%ht))
    call set_segment (this%layout, source, array, this%cell_heat_segid)
  end subroutine HTSD_model_set_cell_heat

  subroutine HTSD_model_set_cell_temp (this, source, array)
    type(HTSD_model), intent(in) :: this
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout) :: array(:)
    INSIST(associated(this%ht))
    call set_segment (this%layout, source, array, this%cell_temp_segid)
  end subroutine HTSD_model_set_cell_temp

  subroutine HTSD_model_set_face_temp (this, source, array)
    type(HTSD_model), intent(in) :: this
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout) :: array(:)
    INSIST(associated(this%ht))
    call set_segment (this%layout, source, array, this%face_temp_segid)
  end subroutine HTSD_model_set_face_temp

  subroutine HTSD_model_set_radiosity (this, index, source, array)
    type(HTSD_model), intent(in) :: this
    integer, intent(in) :: index
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout) :: array(:)
    ASSERT(associated(this%rad_segid))
    ASSERT(index > 0 .and. index <= size(this%rad_segid))
    call set_segment (this%layout, source, array, this%rad_segid(index))
  end subroutine HTSD_model_set_radiosity

  subroutine HTSD_model_set_cell_conc (this, index, source, array)
    type(HTSD_model), intent(in) :: this
    integer, intent(in) :: index
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout) :: array(:)
    ASSERT(associated(this%cell_conc_segid))
    ASSERT(index > 0 .and. index <= size(this%cell_conc_segid))
    call set_segment (this%layout, source, array, this%cell_conc_segid(index))
  end subroutine HTSD_model_set_cell_conc

  subroutine HTSD_model_set_face_conc (this, index, source, array)
    type(HTSD_model), intent(in) :: this
    integer, intent(in) :: index
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout) :: array(:)
    ASSERT(associated(this%face_conc_segid))
    ASSERT(index > 0 .and. index <= size(this%face_conc_segid))
    call set_segment (this%layout, source, array, this%face_conc_segid(index))
  end subroutine HTSD_model_set_face_conc
  
end module HTSD_model_type

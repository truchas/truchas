#include "f90_assert.fpp"

module HTSD_model_type

  use kinds
  use dist_mesh_type
  use mfd_disc_type
  use data_layout_type
  use property_mesh_function
  use source_mesh_function
  use boundary_data
  use interface_data
  use ER_driver
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
    real(r8) :: sbconst, abszero ! Stefan-Boltzmann constant and absolute zero for radiation BC
    !! Enclosure radiation problems
    type(ERD_problem), pointer :: vf_rad_prob(:) => null()
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
    type(dist_mesh), pointer :: mesh => null()
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
  public :: HTSD_model_compute_udot
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
    if (associated(this%vf_rad_prob)) then
      do n = 1, size(this%vf_rad_prob)
        call ERD_problem_destroy (this%vf_rad_prob(n))
      end do
      deallocate(this%vf_rad_prob)
    end if
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
      call this%disc%apply_diff (value, Tcell, Tface, Fcell, Fface)
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
        Fface(n) = Fface(n) + this%mesh%area(n) * this%ht%bc_flux%values(1,j)
      end do

      !! External HTC flux contribution.
      call bd_data_eval (this%ht%bc_htc, t)
      do j = 1, size(this%ht%bc_htc%faces)
        n = this%ht%bc_htc%faces(j)
        Fface(n) = Fface(n) + this%mesh%area(n) * this%ht%bc_htc%values(1,j) * &
                                  (Tface(n) - this%ht%bc_htc%values(2,j))
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
        term = this%mesh%area(n1) * this%ht%ic_htc%values(1,j) * (Tface(n1) - Tface(n2))
        Fface(n1) = Fface(n1) + term
        Fface(n2) = Fface(n2) - term
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
        term = this%mesh%area(n1) * this%ht%ic_rad%values(1,j) * this%ht%sbconst * &
               ((Tface(n1)-this%ht%abszero)**4 - (Tface(n2)-this%ht%abszero)**4)
        Fface(n1) = Fface(n1) + term
        Fface(n2) = Fface(n2) - term
      end do
      if (allocated(void_link)) deallocate(void_link)

      !! Ambient radiation BC flux contribution.
      call bd_data_eval (this%ht%bc_rad, t)
      do j = 1, size(this%ht%bc_rad%faces)
        n = this%ht%bc_rad%faces(j)
        Fface(n) = Fface(n) + this%mesh%area(n) * this%ht%bc_rad%values(1,j) * &
                                this%ht%sbconst * ((Tface(n) - this%ht%abszero)**4 - &
                                           (this%ht%bc_rad%values(2,j)-this%ht%abszero)**4)
      end do

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
          call ERD_compute_heat_flux (this%ht%vf_rad_prob(n), t, qrad, Tface(faces), flux)
          do j = 1, size(faces)
            Fface(faces(j)) = Fface(faces(j)) + this%mesh%area(faces(j)) * flux(j)
          end do
          deallocate(flux)
          !! Residual of the algebraic radiosity system.
          call HTSD_model_get_radiosity_view (this, n, f, fptr)
          call ERD_compute_residual (this%ht%vf_rad_prob(n), t, qrad, Tface(faces), fptr)
          fptr = -fptr
        end do
      end if

      !! Return the on-process part of the heat conduction residuals (the rest is junk)
      call HTSD_model_set_cell_temp (this, Fcell, f)
      call HTSD_model_set_face_temp (this, Fface, f)

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

!TODO!  This procedure needs to be reimplemented to exactly compute the
!TODO!  the initial state and state time derivative instead of the simple
!TODO!  approximation computed here.  Because the exact solution must
!TODO!  involve a solver, this new procedure probably doesn't belong here,
!TODO!  but should become a new class, like HT_precon, that builds on
!TODO!  HT_model.

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! EVAL_UDOT
 !!
 !! Given the values of the cell unknowns, this routine computes values for
 !! the face unknowns and time derivatives for all the unknowns using the
 !! DAE system.  However, instead of solving the linear system for the face
 !! unknowns as we ought, we use simple averaging of adjacent cell unknowns
 !! to get a cheap and decent approximation for the face unknowns.  Likewise
 !! the time derivative of the face unknowns is approximated by a simple
 !! averaging of the computed time derivatives of the cell unknowns.  These
 !! values are intended only to initialize the DAE integrator, which uses
 !! them to extrapolate an initial guess for the first time step, and so the
 !! approximations made are not of great consequence.
 !!

  subroutine HTSD_model_compute_udot (this, t, temp, conc, u, udot)

    type(HTSD_model), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), pointer :: temp(:), conc(:,:)
    real(r8), intent(out), target :: u(:), udot(:)

    integer :: n
    integer, pointer :: faces(:)
    real(r8), allocatable, target :: f(:)
    real(r8), allocatable :: Tface(:), Ccell(:), Cface(:), Tdot(:), dHdT(:)
    real(r8), pointer :: state(:,:), var(:), Hcell(:), Hdot(:), Cdot(:), Fcell(:)
    
    ASSERT(size(u) == layout_size(this%layout))
    ASSERT(size(u) == size(udot))
    
    !! Set cell temperatures.
    if (associated(this%ht)) then
      ASSERT(associated(temp))
      ASSERT(size(temp) == this%mesh%ncell_onP)
      call HTSD_model_set_cell_temp (this, temp, u)
      call HTSD_model_get_cell_temp_view (this, u, var)
      if (associated(this%void_cell)) where (this%void_cell(:size(var))) var = this%void_temp
    end if
    
    !! Set cell concentrations.
    if (associated(this%sd)) then
      ASSERT(associated(conc))
      ASSERT(size(conc,dim=1) == this%mesh%ncell_onP)
      ASSERT(size(conc,dim=2) == this%num_comp)
      do n = 1, this%num_comp
        call HTSD_model_set_cell_conc (this, n, conc(:,n), u)
        call HTSD_model_get_cell_conc_view (this, n, u, var)
        if (associated(this%void_cell)) where (this%void_cell(:size(var))) var = 0.0_r8
      end do
    end if

    !! Here we ought to be solving for the algebraically-coupled face
    !! temperatures and radiosities. Instead we cheat and just approximate
    !! the face temperature by averaging the adjacent cell temperatures, and
    !! compute the radiosities accordingly. The end result is face temperatures
    !! and radiosities that do not give exact heat flux matching at the faces,
    !! and ignores heat flux due to view factor radiation.

    state => HTSD_model_new_state_array(this, u)
    
    if (associated(this%ht)) then
      allocate(Tface(this%mesh%nface))
      call eval_face_averages (state(:,0), Tface)
      !! Overwrite face temperatures with Dirichlet boundary data.
      call bd_data_eval (this%ht%bc_dir, t)
      Tface(this%ht%bc_dir%faces) = this%ht%bc_dir%values(1,:)
      !! Overwrite void face components with dummy value.
      if (associated(this%void_face)) where (this%void_face) Tface = 0.0_r8
      call HTSD_model_set_face_temp (this, Tface, u)

      if (associated(this%ht%vf_rad_prob)) then
        do n = 1, size(this%ht%vf_rad_prob)
          call HTSD_model_get_radiosity_view (this, n, u, var)
          faces => this%ht%vf_rad_prob(n)%faces
          var = 0.0_r8
          call ERD_solve_radiosity (this%ht%vf_rad_prob(n), t, Tface(faces), var)
        end do
      end if
      deallocate(Tface)
    end if
    
    !! Here we ought to be solving the algebraic constraints for the face
    !! concentrations to complete the definition of U.  Instead we cheat and
    !! just approximate them by averaging adjacent cell concentrations.  The
    !! end result is a (somewhat) inconsistent U that does not give exact
    !! concentration flux matching at the faces.

    if (associated(this%sd)) then
      allocate(Cface(this%mesh%nface))
      do n = 1, this%num_comp
        call eval_face_averages (state(:,n), Cface)
        !! Overwrite face concentrations with Dirichlet boundary data.
        call bd_data_eval (this%sd(n)%bc_dir, t)
        Cface(this%sd(n)%bc_dir%faces) = this%sd(n)%bc_dir%values(1,:)
        !! Overwrite void face components with dummy value.
        if (associated(this%void_face)) where (this%void_face) Cface = 0.0_r8
        call HTSD_model_set_face_conc (this, n, Cface, u)
      end do
      deallocate(Cface)
    end if
    
    !! Compute the cell heat density, its time derivative, and the time derivative
    !! of the cell concentrations.  We use the model's F function with appropriate
    !! inputs to do this.  The calculated heat density is exact, but the time
    !! derivatives are not because U is only approximately consistent (see above).
    
    if (associated(this%ht)) then
      call HTSD_model_get_cell_heat_view (this, u, Hcell)
      Hcell = 0.0_r8
    end if
    udot = 0.0_r8
    allocate(f(size(u)))
    call HTSD_model_compute_f (this, t, u, udot, f)
    if (associated(this%ht)) then
      !! Extract the heat densities from F.
      call HTSD_model_get_cell_heat_view (this, u, Hcell)
      call HTSD_model_get_cell_heat_view (this, f, Fcell)
      Hcell = -Fcell
      if (associated(this%void_cell)) where (this%void_cell(:size(Hcell))) Hcell = 0.0_r8
      !! Extract the heat density derivatives from F.
      call HTSD_model_get_cell_heat_view (this, udot, Hdot)
      call HTSD_model_get_cell_temp_view (this, f, Fcell)
      Hdot = -Fcell / this%mesh%volume(:size(Hdot))
      if (associated(this%void_cell)) where (this%void_cell(:size(Hdot))) Hdot = 0.0_r8
    end if
    if (associated(this%sd)) then
      !! Extract the cell concentration derivatives from F.
      do n = 1, this%num_comp
        call HTSD_model_get_cell_conc_view (this, n, udot, Cdot)
        call HTSD_model_get_cell_conc_view (this, n, f, Fcell)
        Cdot = -Fcell / this%mesh%volume(:size(Cdot))
        if (associated(this%void_cell)) where (this%void_cell(:size(Cdot))) Cdot = 0.0_r8
      end do
    end if
    deallocate(f)

    if (associated(this%ht)) then
      !! Time derivative of the cell temperatures.
      allocate(Tdot(this%mesh%ncell), dHdT(this%mesh%ncell))
      call pmf_eval_deriv (this%ht%H_of_T, state, 1, dHdT)
      if (associated(this%void_cell)) where (this%void_cell) dHdT = 1.0_r8
      Tdot(:size(Hdot)) = Hdot / dHdT(:size(Hdot))
      call gather_boundary (this%mesh%cell_ip, Tdot)
      call HTSD_model_set_cell_temp (this, Tdot, udot)
      deallocate(dHdT, state)

      !! Approximate face temperature derivative by average of adjacent cell derivatives.
      allocate(Tface(this%mesh%nface))
      call eval_face_averages (Tdot, Tface)
      call HTSD_model_set_face_temp (this, Tface, udot)
      deallocate(Tdot, Tface)

      !! Here we should compute/approximate the radiosity derivatives.
      !! For now we just leave them 0 (set when udot was set to 0 above).
    end if

    if (associated(this%sd)) then
      !! Approximate face concentration derivative by average of adjacent cell derivatives.
      allocate(Ccell(this%mesh%ncell), Cface(this%mesh%nface))
      do n = 1, this%num_comp
        call HTSD_model_get_cell_conc_copy (this, n, udot, Ccell)
        call gather_boundary (this%mesh%cell_ip, Ccell)
        call eval_face_averages (Ccell, Cface)
        call HTSD_model_set_face_conc (this, n, Cface, udot)
      end do
      deallocate(Ccell, Cface)
    end if

  contains

    subroutine eval_face_averages (ucell, uface)
    
      use bitfield_type
      use index_partitioning

      real(r8), intent(in)  :: ucell(:)
      real(r8), intent(out) :: uface(:)

      integer :: j
      integer :: scale(size(uface))

      ASSERT(size(ucell) == this%mesh%ncell)
      ASSERT(size(uface) == this%mesh%nface)

      uface = 0.0_r8
      scale = 0.0_r8
      do j = 1, this%mesh%ncell
        if (associated(this%void_cell)) then
          if (this%void_cell(j)) cycle
        end if
        uface(this%mesh%cface(:,j)) = uface(this%mesh%cface(:,j)) + ucell(j)
        scale(this%mesh%cface(:,j)) = scale(this%mesh%cface(:,j)) + 1
      end do
      call gather_boundary (this%mesh%face_ip, uface)
      call gather_boundary (this%mesh%face_ip, scale)
      
      where (scale == 0)
        uface = 0.0_r8
      elsewhere
        uface = uface / scale
      end where

    end subroutine eval_face_averages

  end subroutine HTSD_model_compute_udot

  
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

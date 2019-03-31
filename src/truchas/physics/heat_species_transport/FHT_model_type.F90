!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module FHT_model_type

  use kinds
  use mfd_disc_type
  use unstr_mesh_type
  use data_layout_type
  use property_mesh_function
  use source_mesh_function
  use bndry_func1_class
  use bndry_func2_class
  use intfc_func2_class
  use rad_problem_type
  use index_partitioning
  use truchas_timers
  implicit none
  private

  type, public :: FHT_model
    type(mfd_disc),   pointer :: disc => null()
    type(unstr_mesh), pointer :: mesh => null()
    type(data_layout) :: layout
    integer :: cell_temp_segid, face_temp_segid
    integer, pointer :: rad_segid(:) => null()
    logical, pointer :: void_cell(:) => null(), void_face(:) => null()
    !real(r8), pointer :: vfrac(:) => null()
    !! The remaining components must be defined BEFORE calling the init function.
    !! Heat equation parameters
    type(prop_mf), pointer :: conductivity => null()   ! thermal conductivity
    type(prop_mf), pointer :: H_of_T => null()          ! enthalpy as a function of temperature
    type(source_mf), pointer :: q => null()    ! external heat source
    !! Boundary condition data
    class(bndry_func1), allocatable :: bc_dir  ! Dirichlet
    class(bndry_func1), allocatable :: bc_flux ! simple flux
    class(bndry_func2), allocatable :: bc_htc  ! external HTC (coef, ref temp)
    class(bndry_func2), allocatable :: bc_rad  ! simple radiation (eps, amb temp)
    class(intfc_func2), allocatable :: ic_htc  ! internal HTC
    class(intfc_func2), allocatable :: ic_rad  ! internal gap radiation
    !! Enclosure radiation problems
    type(rad_problem), pointer :: vf_rad_prob(:) => null()
  end type FHT_model

  public :: FHT_model_init
  public :: FHT_model_delete
  public :: FHT_model_compute_f
  public :: FHT_model_compute_udot

  public :: FHT_model_size
  public :: FHT_model_get_cell_temp_view
  public :: FHT_model_get_face_temp_view
  public :: FHT_model_get_radiosity_view

  public :: FHT_model_get_cell_temp_copy
  public :: FHT_model_get_face_temp_copy
  public :: FHT_model_get_radiosity_copy

  public :: FHT_model_set_cell_temp
  public :: FHT_model_set_face_temp
  public :: FHT_model_set_radiosity

  public :: FHT_model_cell_temp_index
  public :: FHT_model_face_temp_index

contains

  subroutine FHT_model_init (this, disc)

    type(FHT_model), intent(inout) :: this ! must be inout because some components set beforehand
    type(mfd_disc), intent(in), target :: disc

    integer :: n

    this%disc => disc
    this%mesh => disc%mesh

    !! Create the packed layout of the model variables.
    this%cell_temp_segid = alloc_segment(this%layout, this%mesh%ncell_onP)
    this%face_temp_segid = alloc_segment(this%layout, this%mesh%nface_onP)
    if (associated(this%vf_rad_prob)) then
      ASSERT(size(this%vf_rad_prob) > 0)
      allocate(this%rad_segid(size(this%vf_rad_prob)))
      do n = 1, size(this%vf_rad_prob)
        this%rad_segid(n) = alloc_segment(this%layout, size(this%vf_rad_prob(n)%faces))
      end do
    end if
    call alloc_complete (this%layout)

  end subroutine FHT_model_init
  
  subroutine FHT_model_delete (this)
    type(FHT_model), intent(inout) :: this
    integer :: n
    this%disc => null()
    this%mesh => null()
    call delete_layout (this%layout)
    if (associated(this%rad_segid)) deallocate(this%rad_segid)
    this%void_cell => null()
    this%void_face => null()
    if (associated(this%conductivity)) then
      call destroy (this%conductivity)
      deallocate(this%conductivity)
    end if
    if (associated(this%H_of_T)) then
      call destroy (this%H_of_T)
      deallocate(this%H_of_T)
    end if
    if (associated(this%q)) then
      call smf_destroy (this%q)
      deallocate(this%q)
    end if
    if (associated(this%vf_rad_prob)) deallocate(this%vf_rad_prob)
  end subroutine FHT_model_delete

  subroutine FHT_model_compute_f (this, t, u, hdot, f)

    type(FHT_model), intent(inout) :: this
    real(r8), intent(in)  :: t, u(:), hdot(:)
    real(r8), intent(out) :: f(:)
    target :: u, f

    integer :: j, n, n1, n2
    real(r8) :: term
    real(r8), target :: state(this%mesh%ncell,1)
    real(r8), pointer :: rptr(:), qrad(:), Tcell(:)
    real(r8), dimension(this%mesh%ncell) :: Fcell, value
    real(r8), dimension(this%mesh%nface) :: Fface, Tface
    real(r8), allocatable :: Fdir(:), flux(:)
    integer, pointer :: faces(:)
    logical, allocatable :: void_link(:)

    call start_timer ('FHT function')

    !! Initialize the STATE array.
    call FHT_model_get_cell_temp_copy (this, u, state(:,1))
    call gather_boundary (this%mesh%cell_ip, state(:,1))

  !!!! RESIDUAL OF THE HEAT EQUATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! Off-process-extended cell and face temperatures.
    Tcell => state(:,1)
    call FHT_model_get_face_temp_copy (this, u, Tface)
    call gather_boundary (this%mesh%face_ip, Tface)

    !! Pre-compute the Dirichlet condition residual and
    !! impose the Dirichlet data on the face temperature.
    if (allocated(this%bc_dir)) then
      call this%bc_dir%compute(t)
      allocate(Fdir(size(this%bc_dir%index)))
      do j = 1, size(this%bc_dir%index)
        Fdir(j) = Tface(this%bc_dir%index(j)) - this%bc_dir%value(j)
        Tface(this%bc_dir%index(j)) = this%bc_dir%value(j)
      end do
    end if

    !! Compute the generic heat equation residual.
    call pmf_eval (this%conductivity, state, value)
    where (this%void_cell) value = 0.0_r8
    call this%disc%apply_diff (value, Tcell, Tface, Fcell, Fface)

    !! Add the source and time deriviative contribution.
    !! The result is complete on on-process cells only, but the
    !! off-process values are not used hereafter.
    call smf_eval (this%q, t, value)
    do j = 1, this%mesh%ncell_onP
      Fcell(j) = Fcell(j) + this%mesh%volume(j)*(hdot(j) - value(j))
    end do
    !call gather_boundary (this%mesh%cell_ip, Fcell) ! off-process not used below

    !! Dirichlet condition residuals.
    if (allocated(this%bc_dir)) then
      do j = 1, size(this%bc_dir%index)
        n = this%bc_dir%index(j)
        Fface(n) = Fdir(j)  !! overwrite with pre-computed values
      end do
      deallocate(Fdir)
    end if

    !! Simple flux BC contribution.
    if (allocated(this%bc_flux)) then
      call this%bc_flux%compute(t)
      do j = 1, size(this%bc_flux%index)
        n = this%bc_flux%index(j)
        Fface(n) = Fface(n) + this%mesh%area(n) * this%bc_flux%value(j)
      end do
    end if

    !! External HTC flux contribution.
    if (allocated(this%bc_htc)) then
    call this%bc_htc%compute(t, Tface)
      do j = 1, size(this%bc_htc%index)
        n = this%bc_htc%index(j)
        Fface(n) = Fface(n) + this%bc_htc%value(j)
      end do
    end if

    !! Internal HTC flux contribution.
    if (allocated(this%ic_htc)) then
      call this%ic_htc%compute(t, Tface)
      allocate(void_link(size(this%ic_htc%index,2)))
      do j = 1, size(void_link)
        void_link(j) = any(this%void_face(this%ic_htc%index(:,j)))
      end do
      do j = 1, size(this%ic_htc%index,2)
        if (void_link(j)) cycle
        n1 = this%ic_htc%index(1,j)
        n2 = this%ic_htc%index(2,j)
        Fface(n1) = Fface(n1) + this%ic_htc%value(j)
        Fface(n2) = Fface(n2) - this%ic_htc%value(j)
      end do
      deallocate(void_link)
    end if
      
    !! Internal gap radiation contribution.
    if (allocated(this%ic_rad)) then
      call this%ic_rad%compute(t, Tface)
      allocate(void_link(size(this%ic_rad%index,2)))
      do j = 1, size(void_link)
        void_link(j) = any(this%void_face(this%ic_rad%index(:,j)))
      end do
      do j = 1, size(this%ic_rad%index,2)
        if (void_link(j)) cycle
        n1 = this%ic_rad%index(1,j)
        n2 = this%ic_rad%index(2,j)
        Fface(n1) = Fface(n1) + this%ic_rad%value(j)
        Fface(n2) = Fface(n2) - this%ic_rad%value(j)
      end do
      deallocate(void_link)
    end if

    !! Ambient radiation BC flux contribution.
    if (allocated(this%bc_rad)) then
      call this%bc_rad%compute(t, Tface)
      do j = 1, size(this%bc_rad%index)
        n = this%bc_rad%index(j)
        Fface(n) = Fface(n) + this%bc_rad%value(j)
      end do
    end if

    !! Overwrite function value on void cells and faces with dummy equation T=0.
    where (this%void_cell) Fcell = Tcell
    where (this%void_face) Fface = Tface

  !!!! RESIDUALS OF THE ENCLOSURE RADIATION SYSTEMS !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (associated(this%vf_rad_prob)) then
      do n = 1, size(this%vf_rad_prob)
        call FHT_model_get_radiosity_view (this, n, u, qrad)
        faces => this%vf_rad_prob(n)%faces
        !! Radiative heat flux contribution to the heat conduction face residual.
        allocate(flux(size(faces)))
        call this%vf_rad_prob(n)%heat_flux (t, qrad, Tface(faces), flux)
        do j = 1, size(faces)
          Fface(faces(j)) = Fface(faces(j)) + this%mesh%area(faces(j)) * flux(j)
        end do
        deallocate(flux)
        !! Residual of the algebraic radiosity system.
        call FHT_model_get_radiosity_view (this, n, f, rptr)
        call this%vf_rad_prob(n)%residual (t, qrad, Tface(faces), rptr)
        rptr = -rptr
      end do
    end if

    !! Return the on-process part of the heat conduction residuals (the rest is junk)
    call FHT_model_set_cell_temp (this, Fcell, f)
    call FHT_model_set_face_temp (this, Fface, f)

    call stop_timer ('FHT function')

  end subroutine FHT_model_compute_f

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! EVAL_UDOT
 !!
 !! Given initial cell temperatures, this routine "computes" the corresponding
 !! face temperatures and radiosities to complete the initial state vector,
 !! and also computes an approximation to the initial state time derivative.
 !!
 !! The results computed by this procedure are needed to set the initial state
 !! of the ODE integrator.  However the computed values are not of essential
 !! importance.  The result of the first time step depends only the heat
 !! densities corresponding to the initial cell temperatures, which we take
 !! as given, and not the other components of state computed here.  They are
 !! used only to generate the initial guess for the first time step.  It is
 !! somewhat important that the initial state be decently accurate so that the
 !! nonlinear solver has a chance of converging, but the influence of its
 !! time derivative can be made insignificant if one is willing to use a
 !! small enough initial time step (not always an acceptable approach).
 !!
 !! Currently we do something really crude.
 !!
 !! TODO: Regardless of above, we really want to start with a fully consistent
 !! initial state.  The enables correct interpolation for output, but more
 !! importantly will allow us to use the last step size when restarting
 !! calculations.  Computing the fully consistent initial state requires us
 !! to solve the algebraic components of the DAE system.  If view factor
 !! radiation is not involved, this is generally a linear problem (though
 !! simple radiation BC makes it nonlinear), but nonlinear otherwise.  Such
 !! a solver should be implemented in a separate class -- the FHT model
 !! should not involve solver parameters.
 !!
 !! TODO: Once we have a fully consistent initial state we can worry about
 !! its time derivative.  This is somewhat more difficult.  Instead of
 !! differentiating the the DAE system to get equations for the derivative
 !! of the algebraic variables, a simpler approach might be to take a small
 !! forward Euler step: advance the heat density, cell-by-cell solve for
 !! temperatures, then solve for the consistent face temperatures and
 !! radiosities exactly as was done to get the initial state.  This gives
 !! two states, whose difference gives us an approximation for the derivative.
 !!

  subroutine FHT_model_compute_udot (this, t, temp, u, udot)

    use index_partitioning

    type(FHT_model), intent(inout) :: this
    real(r8), intent(in)  :: t, temp(:)
    real(r8), intent(out) :: u(:), udot(:)
    target :: u, udot

    integer :: j, index, stat, numitr
    integer, pointer :: faces(:)
    real(r8), target :: f(size(u))
    real(r8) :: fdinc, H0, H1, error
    real(r8), pointer :: var(:), Fcell(:)
    real(r8), allocatable :: Tcell(:), Tface(:), Hdot(:)

    ASSERT(size(temp) == this%mesh%ncell_onP)
    ASSERT(size(u) == layout_size(this%layout))
    ASSERT(size(u) == size(udot))

    !! Set cell temperatures.
    call FHT_model_set_cell_temp (this, temp, u)
    call FHT_model_get_cell_temp_view (this, u, var)
    where (this%void_cell(:size(var))) var = 0.0_r8

    !! Here we ought to be solving for the remaining algebraically-coupled
    !! components of U, namely the face temperatures and radiosities. Instead
    !! we cheat and just approximate the face temperature by averaging the
    !! adjacent cell temperatures, and compute the radiosities accordingly.
    !! The end result is face temperatures and radiosities that do not give
    !! exact heat flux matching at the faces, and ignores heat flux due to
    !! view factor radiation.

    allocate(Tcell(this%mesh%ncell), Tface(this%mesh%nface))
    call FHT_model_get_cell_temp_copy (this, u, Tcell)
    call gather_boundary (this%mesh%cell_ip, Tcell)
    call eval_face_averages (Tcell, Tface)
    call FHT_model_set_face_temp(this, Tface, u)

    if (associated(this%vf_rad_prob)) then
      do index = 1, size(this%vf_rad_prob)
        call FHT_model_get_radiosity_view (this, index, u, var)
        faces => this%vf_rad_prob(index)%faces
        var = 0.0_r8
        call this%vf_rad_prob(index)%solve_radiosity (t, Tface(faces), var, stat, numitr, error)
        !! Go ahead and set radiosity derivatives to zero now.
        call FHT_model_get_radiosity_view (this, index, udot, var)
        var = 0.0_r8
      end do
    end if

    !! Compute the instantaneous time derivative of the heat density.
    !! This would be exact if U were consistent (it isn't).
    allocate(Hdot(this%mesh%ncell_onP))
    Hdot = 0.0_r8
    call FHT_model_compute_f (this, t, u, Hdot, f)
    call FHT_model_get_cell_temp_view (this, f, Fcell)
    do j = 1, this%mesh%ncell_onP
      if (this%void_cell(j)) then
        Hdot(j) = 0.0_r8
      else
        Hdot(j) = -Fcell(j) / this%mesh%volume(j)
      end if
    end do

    !! Time derivative of the cell temperatures; overwrite Tcell with the result.
    fdinc = sqrt(epsilon(1.0d0))  !TODO: FIXME
    do j = 1, this%mesh%ncell_onP
      if (this%void_cell(j)) then
        Tcell(j) = 0.0_r8
      else
        call pmf_eval (this%H_of_T, j, Tcell(j:j), H0)
        call pmf_eval (this%H_of_T, j, Tcell(j:j)+fdinc, H1)
        Tcell(j) = Hdot(j) / ((H1 - H0) / fdinc)
      end if
    end do
    call FHT_model_set_cell_temp (this, Tcell, udot)

    !! Approximate face derivative by average of adjacent cell derivatives.
    call gather_boundary (this%mesh%cell_ip, Tcell)
    call eval_face_averages (Tcell, Tface)
    call FHT_model_set_face_temp (this, Tface, udot)

  contains

    subroutine eval_face_averages (ucell, uface)

      real(r8), intent(in)  :: ucell(:)
      real(r8), intent(out) :: uface(:)

      integer :: j, scale(size(uface))

      ASSERT(size(ucell) == this%mesh%ncell)
      ASSERT(size(uface) == this%mesh%nface)

      uface = 0.0_r8
      scale = 0
      do j = 1, this%mesh%ncell
        if (this%void_cell(j)) cycle
        associate (cface => this%mesh%cface(this%mesh%xcface(j):this%mesh%xcface(j+1)-1))
          uface(cface) = uface(cface) + ucell(j)
          scale(cface) = scale(cface) + 1
        end associate
      end do
      call gather_boundary (this%mesh%face_ip, uface)
      call gather_boundary (this%mesh%face_ip, scale)
      do j = 1, this%mesh%nface
        if (scale(j) == 0) then
          uface(j) = 0.0_r8
        else
          uface(j) = uface(j) / scale(j)
        end if
      end do

    end subroutine eval_face_averages

  end subroutine FHT_model_compute_udot

  pure integer function FHT_model_size (this)
    type(FHT_model), intent(in) :: this
    FHT_model_size = layout_size(this%layout)
  end function FHT_model_size

  subroutine FHT_model_get_cell_temp_view (this, array, view)
    type(FHT_model), intent(in) :: this
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    call get_segment_view (this%layout, array, this%cell_temp_segid, view)
  end subroutine FHT_model_get_cell_temp_view

  subroutine FHT_model_get_face_temp_view (this, array, view)
    type(FHT_model), intent(in) :: this
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    call get_segment_view (this%layout, array, this%face_temp_segid, view)
  end subroutine FHT_model_get_face_temp_view

  subroutine FHT_model_get_radiosity_view (this, index, array, view)
    type(FHT_model), intent(in) :: this
    integer, intent(in) :: index
    real(r8), target, intent(in) :: array(:)
    real(r8), pointer :: view(:)
    ASSERT(associated(this%rad_segid))
    ASSERT(index > 0 .and. index <= size(this%rad_segid))
    call get_segment_view (this%layout, array, this%rad_segid(index), view)
  end subroutine FHT_model_get_radiosity_view

  subroutine FHT_model_get_cell_temp_copy (this, array, copy)
    type(FHT_model), intent(in) :: this
    real(r8), intent(in) :: array(:)
    real(r8), intent(inout) :: copy(:)
    call get_segment_copy (this%layout, array, this%cell_temp_segid, copy)
  end subroutine FHT_model_get_cell_temp_copy

  subroutine FHT_model_get_face_temp_copy (this, array, copy)
    type(FHT_model), intent(in) :: this
    real(r8), intent(in) :: array(:)
    real(r8), intent(inout) :: copy(:)
    call get_segment_copy (this%layout, array, this%face_temp_segid, copy)
  end subroutine FHT_model_get_face_temp_copy

  subroutine FHT_model_get_radiosity_copy (this, index, array, copy)
    type(FHT_model), intent(in) :: this
    integer, intent(in) :: index
    real(r8), intent(in) :: array(:)
    real(r8), intent(inout) :: copy(:)
    ASSERT(associated(this%rad_segid))
    ASSERT(index > 0 .and. index <= size(this%rad_segid))
    call get_segment_copy (this%layout, array, this%rad_segid(index), copy)
  end subroutine FHT_model_get_radiosity_copy

  subroutine FHT_model_set_cell_temp (this, source, array)
    type(FHT_model), intent(in) :: this
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout) :: array(:)
    call set_segment (this%layout, source, array, this%cell_temp_segid)
  end subroutine FHT_model_set_cell_temp

  subroutine FHT_model_set_face_temp (this, source, array)
    type(FHT_model), intent(in) :: this
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout) :: array(:)
    call set_segment (this%layout, source, array, this%face_temp_segid)
  end subroutine FHT_model_set_face_temp

  subroutine FHT_model_set_radiosity (this, index, source, array)
    type(FHT_model), intent(in) :: this
    integer, intent(in) :: index
    real(r8), intent(in) :: source(:)
    real(r8), intent(inout) :: array(:)
    ASSERT(associated(this%rad_segid))
    ASSERT(index > 0 .and. index <= size(this%rad_segid))
    call set_segment (this%layout, source, array, this%rad_segid(index))
  end subroutine FHT_model_set_radiosity

  integer function FHT_model_cell_temp_index (this, n) result (index)
    type(FHT_model), intent(in) :: this
    integer, intent(in) :: n
    index = layout_index(this%layout, this%cell_temp_segid, n)
  end function FHT_model_cell_temp_index

  integer function FHT_model_face_temp_index (this, n) result (index)
    type(FHT_model), intent(in) :: this
    integer, intent(in) :: n
    index = layout_index(this%layout, this%face_temp_segid, n)
  end function FHT_model_face_temp_index

end module FHT_model_type

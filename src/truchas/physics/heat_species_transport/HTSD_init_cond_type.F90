!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module HTSD_init_cond_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use unstr_mesh_type
  use data_layout_type
  use HTSD_model_type
  use truchas_logging_services
  use parallel_communication
  use parameter_list_type
  implicit none
  private
  
  type, public :: HTSD_init_cond
    private
    type(unstr_mesh), pointer :: mesh  => null()  ! reference only -- do not own
    type(HTSD_model), pointer :: model => null()  ! reference only -- do not own
    type(parameter_list), pointer :: params => null() ! reference only -- do not own
  contains
    procedure :: init
    procedure :: compute
    procedure :: compute_udot
  end type HTSD_init_cond

contains

  subroutine init (this, model, params)

    class(HTSD_init_cond), intent(out) :: this
    type(HTSD_model), intent(in), target :: model
    type(parameter_list), pointer :: params
    
    this%model => model
    this%mesh  => model%mesh
    this%params => params
    
  end subroutine init


  subroutine compute (this, t, temp, conc, u, udot)

    class(HTSD_init_cond), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), optional :: temp(:), conc(:,:)
    real(r8), intent(out), target :: u(:), udot(:)

    integer :: n
    real(r8), allocatable, target :: f(:)
    real(r8), allocatable :: Tface(:), Ccell(:), Cface(:), Tdot(:), dHdT(:)
    real(r8), pointer :: state(:,:), var(:), Hcell(:), Hdot(:), Cdot(:), Fcell(:)
    real(r8) :: dt
    
    ASSERT(size(u) == layout_size(this%model%layout))
    ASSERT(size(u) == size(udot))
    
    call TLS_info ('')
    call TLS_info ('Computing consistent initial state for HT/SD solver ...')
    
    u = 0.0_r8
    
    !! Set cell temperatures.
    if (associated(this%model%ht)) then
      ASSERT(present(temp))
      ASSERT(size(temp) == this%mesh%ncell_onP)
      call HTSD_model_set_cell_temp (this%model, temp, u)
      call HTSD_model_get_cell_temp_view (this%model, u, var)
      if (associated(this%model%void_cell)) &
          where (this%model%void_cell(:size(var))) var = this%model%void_temp
    end if
    
    !! Set cell concentrations.
    if (associated(this%model%sd)) then
      ASSERT(present(conc))
      ASSERT(size(conc,dim=1) == this%mesh%ncell_onP)
      ASSERT(size(conc,dim=2) == this%model%num_comp)
      do n = 1, this%model%num_comp
        call HTSD_model_set_cell_conc (this%model, n, conc(:,n), u)
        call HTSD_model_get_cell_conc_view (this%model, n, u, var)
        if (associated(this%model%void_cell)) &
            where (this%model%void_cell(:size(var))) var = 0.0_r8
      end do
    end if

    !! State variables (cell temps and concs) on which model parameters depend.
    state => HTSD_model_new_state_array(this%model, u)
    
    !! Set consistent face temperatures and radiosities.
    if (associated(this%model%ht)) then
      !! Approximate face temps by averaging adjacent cell temps; provides
      !! a cheap initial guess for the solution procedure that follows.
      allocate(Tface(this%mesh%nface))
      call average_to_faces (this%mesh, state(:,0), Tface, this%model%void_cell)
      !! Overwrite face temperatures with Dirichlet boundary data.
      if (allocated(this%model%ht%bc_dir)) then
        call this%model%ht%bc_dir%compute(t)
        Tface(this%model%ht%bc_dir%index) = this%model%ht%bc_dir%value
      end if
      !! Overwrite void face components with dummy value.
      if (associated(this%model%void_face)) where (this%model%void_face) Tface = this%model%void_temp
      call HTSD_model_set_face_temp (this%model, Tface, u)
      deallocate(Tface)
      !! Now solve the algebraic equations for face temps and radiosities.
      call compute_face_temp (this%model, t, state, u, this%params)
    end if
    
    !! Set face concentrations.  Here we ought to be solving the algebraic
    !! equations of species flux continuity for the consistent face concs
    !! as we do above for face temps, but instead we just approximate them
    !! by averaging adjacent cell concs.  (FIXME)
    if (associated(this%model%sd)) then
      allocate(Cface(this%mesh%nface))
      do n = 1, this%model%num_comp
        !! Approximate face concs by averaging adjacent cell concs.
        call average_to_faces (this%mesh, state(:,n), Cface, this%model%void_cell)
        !! Overwrite face concentrations with Dirichlet boundary data.
        if (allocated(this%model%sd(n)%bc_dir)) then
          call this%model%sd(n)%bc_dir%compute(t)
          Cface(this%model%sd(n)%bc_dir%index) = this%model%sd(n)%bc_dir%value
        end if
        !! Overwrite void face components with dummy value.
        if (associated(this%model%void_face)) where (this%model%void_face) Cface = 0.0_r8
        call HTSD_model_set_face_conc (this%model, n, Cface, u)
      end do
      deallocate(Cface)
    end if
    
    !! Set the cell heat density.
    if (associated(this%model%ht)) then
      call HTSD_model_get_cell_heat_view (this%model, u, Hcell)
      Hcell = 0.0_r8
    end if
    udot = 0.0_r8
    allocate(f(size(u)))
    call HTSD_model_compute_f (this%model, t, u, udot, f)
    if (associated(this%model%ht)) then
      !! Extract the heat densities from F.
      call HTSD_model_get_cell_heat_view (this%model, u, Hcell)
      call HTSD_model_get_cell_heat_view (this%model, f, Fcell)
      Hcell = -Fcell
      if (associated(this%model%void_cell)) &
          where (this%model%void_cell(:size(Hcell))) Hcell = 0.0_r8
    end if
    
    call TLS_info ('')
    call TLS_info ('Computing consistent initial state derivative for HT/SD solver ...')

    !! Set the cell heat density time derivative.
    if (associated(this%model%ht)) then
      !! Extract the heat density derivatives from F.
      call HTSD_model_get_cell_heat_view (this%model, udot, Hdot)
      call HTSD_model_get_cell_temp_view (this%model, f, Fcell)
      Hdot = -Fcell / this%mesh%volume(:size(Hdot))
      if (associated(this%model%void_cell)) &
          where (this%model%void_cell(:size(Hdot))) Hdot = 0.0_r8
    end if
    
    !! Set the cell concentration time derivative.  This is approximate
    !! because it depends on the face concentrations which are approximate.
    if (associated(this%model%sd)) then
      !! Extract the cell concentration derivatives from F.
      do n = 1, this%model%num_comp
        call HTSD_model_get_cell_conc_view (this%model, n, udot, Cdot)
        call HTSD_model_get_cell_conc_view (this%model, n, f, Fcell)
        Cdot = -Fcell / this%mesh%volume(:size(Cdot))
        if (associated(this%model%void_cell)) &
            where (this%model%void_cell(:size(Cdot))) Cdot = 0.0_r8
      end do
    end if

    !! Compute the cell temperature derivative from the temperature/enthalpy
    !! relation given the heat density derivative
    if (associated(this%model%ht)) then
      allocate(Tdot(this%mesh%ncell), dHdT(this%mesh%ncell))
      call this%model%ht%H_of_T%compute_deriv(state, 1, dHdT)
      if (associated(this%model%void_cell)) &
          where (this%model%void_cell) dHdT = 1.0_r8
      Tdot(:size(Hdot)) = Hdot / dHdT(:size(Hdot))
      call this%mesh%cell_imap%gather_offp(Tdot)
      call HTSD_model_set_cell_temp (this%model, Tdot, udot)
      deallocate(Tdot,dHdT)
    end if

    !! Approximate the initial time derivatives of the remaining algebraic
    !! variables (face temps, concs, and radiosities) by a finite difference.
    !! The temperatures and concentrations are advanced by a small time step
    !! using their time derivatives (forward Euler), and then the associated
    !! cell temps and face temps, concs, and radiosities solved for.
    
    call this%params%get ('dt', dt)
    f = u + dt * udot
    
    !! New state variables (cell temps and concs) on which model parameters depend.
    deallocate(state)
    state => HTSD_model_new_state_array(this%model, u)
    
    !! Set consistent advanced face temps and radiosities.  Use their initial
    !! conditions as the initial guess for the solution procedure.
    if (associated(this%model%ht)) then
      call compute_face_temp (this%model, t+dt, state, f, this%params)
    end if

    !! We should deal with face concs in the same manner as face temps, but
    !! for now we continue to use the simple approx of the initial face concs
    !! done later. 
    !if (associated(this%model%sd)) then
    !  do n = 1, this%model%num_comp
    !    call compute_face_conc (this%model, index, t+dt, state, f, this%params)
    !  end do
    !end if
    
    f = (f - u) / dt
    
    !! Finite difference approx of face temp and radiosity derivatives.
    if (associated(this%model%ht)) then
      call HTSD_model_get_face_temp_view (this%model, f, var)
      call HTSD_model_set_face_temp (this%model, var, udot)
      if (associated(this%model%ht%vf_rad_prob)) then
        do n = 1, size(this%model%ht%vf_rad_prob)
          call HTSD_model_get_radiosity_view (this%model, n, f, var)
          call HTSD_model_set_radiosity (this%model, n, var, udot)
        end do
      end if
    end if

    !! We should compute a similar finite difference approximation of the face
    !! concs, but for now we continue to approximate the face conc derivatives
    !! by averaging the adjacent cell conc derivatives.
    if (associated(this%model%sd)) then
      !! Approximate face concentration derivative by average of adjacent cell derivatives.
      allocate(Ccell(this%mesh%ncell), Cface(this%mesh%nface))
      do n = 1, this%model%num_comp
        call HTSD_model_get_cell_conc_copy (this%model, n, udot, Ccell)
        call this%mesh%cell_imap%gather_offp(Ccell)
        call average_to_faces (this%mesh, Ccell, Cface, this%model%void_cell)
        call HTSD_model_set_face_conc (this%model, n, Cface, udot)
      end do
      deallocate(Ccell, Cface)
    end if
    
    deallocate(state)

  end subroutine compute


  subroutine compute_udot (this, t, u, udot)

    class(HTSD_init_cond), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in),  target :: u(:)
    real(r8), intent(out), target :: udot(:)

    integer :: n
    real(r8), allocatable, target :: f(:)
    real(r8), allocatable :: Ccell(:), Cface(:), Tdot(:), dHdT(:)
    real(r8), pointer :: state(:,:), var(:), Hdot(:), Cdot(:), Fcell(:)
    real(r8) :: dt

    ASSERT(size(u) == layout_size(this%model%layout))
    ASSERT(size(u) == size(udot))

    call TLS_info ('')
    call TLS_info ('Computing consistent initial state derivative for HT/SD solver ...')

    !! State variables (cell temps and concs) on which model parameters depend.
    state => HTSD_model_new_state_array(this%model, u)

    udot = 0.0_r8
    allocate(f(size(u)))
    call HTSD_model_compute_f (this%model, t, u, udot, f)

    !! Set the cell heat density time derivative.
    if (associated(this%model%ht)) then
      !! Extract the heat density derivatives from F.
      call HTSD_model_get_cell_heat_view (this%model, udot, Hdot)
      call HTSD_model_get_cell_temp_view (this%model, f, Fcell)
      Hdot = -Fcell / this%mesh%volume(:size(Hdot))
      if (associated(this%model%void_cell)) &
          where (this%model%void_cell(:size(Hdot))) Hdot = 0.0_r8
    end if

    !! Set the cell concentration time derivative.  This is approximate
    !! because it depends on the face concentrations which are approximate.
    if (associated(this%model%sd)) then
      !! Extract the cell concentration derivatives from F.
      do n = 1, this%model%num_comp
        call HTSD_model_get_cell_conc_view (this%model, n, udot, Cdot)
        call HTSD_model_get_cell_conc_view (this%model, n, f, Fcell)
        Cdot = -Fcell / this%mesh%volume(:size(Cdot))
        if (associated(this%model%void_cell)) &
            where (this%model%void_cell(:size(Cdot))) Cdot = 0.0_r8
      end do
    end if

    !! Compute the cell temperature derivative from the temperature/enthalpy
    !! relation given the heat density derivative
    if (associated(this%model%ht)) then
      allocate(Tdot(this%mesh%ncell), dHdT(this%mesh%ncell))
      call this%model%ht%H_of_T%compute_deriv(state, 1, dHdT)
      if (associated(this%model%void_cell)) &
          where (this%model%void_cell) dHdT = 1.0_r8
      Tdot(:size(Hdot)) = Hdot / dHdT(:size(Hdot))
      call this%mesh%cell_imap%gather_offp(Tdot)
      call HTSD_model_set_cell_temp (this%model, Tdot, udot)
      deallocate(Tdot,dHdT)
    end if

    !! Approximate the initial time derivatives of the remaining algebraic
    !! variables (face temps, concs, and radiosities) by a finite difference.
    !! The temperatures and concentrations are advanced by a small time step
    !! using their time derivatives (forward Euler), and then the associated
    !! cell temps and face temps, concs, and radiosities solved for.

    call this%params%get ('dt', dt)
    f = u + dt * udot

    !! New state variables (cell temps and concs) on which model parameters depend.
    deallocate(state)
    state => HTSD_model_new_state_array(this%model, u)

    !! Set consistent advanced face temps and radiosities.  Use their initial
    !! conditions as the initial guess for the solution procedure.
    if (associated(this%model%ht)) then
      call compute_face_temp (this%model, t+dt, state, f, this%params)
    end if

    !! We should deal with face concs in the same manner as face temps, but
    !! for now we continue to use the simple approx of the initial face concs
    !! done later.
    !if (associated(this%model%sd)) then
    !  do n = 1, this%model%num_comp
    !    call compute_face_conc (this%model, index, t+dt, state, f, this%params)
    !  end do
    !end if

    f = (f - u) / dt

    !! Finite difference approx of face temp and radiosity derivatives.
    if (associated(this%model%ht)) then
      call HTSD_model_get_face_temp_view (this%model, f, var)
      call HTSD_model_set_face_temp (this%model, var, udot)
      if (associated(this%model%ht%vf_rad_prob)) then
        do n = 1, size(this%model%ht%vf_rad_prob)
          call HTSD_model_get_radiosity_view (this%model, n, f, var)
          call HTSD_model_set_radiosity (this%model, n, var, udot)
        end do
      end if
    end if

    !! We should compute a similar finite difference approximation of the face
    !! concs, but for now we continue to approximate the face conc derivatives
    !! by averaging the adjacent cell conc derivatives.
    if (associated(this%model%sd)) then
      !! Approximate face concentration derivative by average of adjacent cell derivatives.
      allocate(Ccell(this%mesh%ncell), Cface(this%mesh%nface))
      do n = 1, this%model%num_comp
        call HTSD_model_get_cell_conc_copy (this%model, n, udot, Ccell)
        call this%mesh%cell_imap%gather_offp(Ccell)
        call average_to_faces (this%mesh, Ccell, Cface, this%model%void_cell)
        call HTSD_model_set_face_conc (this%model, n, Cface, udot)
      end do
      deallocate(Ccell, Cface)
    end if

    deallocate(state)

  end subroutine compute_udot

  !! Given a cell-based field, this auxiliary subroutine produces a face-based
  !! field whose value on each face is the average of the values on adjacent
  !! cells.  SKIP is an optional cell-based mask array identifying cells whose
  !! values are to be ignored.  Faces with no adjacent cells have the value 0.

  subroutine average_to_faces (mesh, ucell, uface, skip)

    use bitfield_type

    type(unstr_mesh), intent(in) :: mesh
    real(r8), intent(in)  :: ucell(:)
    real(r8), intent(out) :: uface(:)
    logical, pointer, intent(in) :: skip(:)

    integer :: j
    integer :: scale(size(uface))

    ASSERT(size(ucell) == mesh%ncell)
    ASSERT(size(uface) == mesh%nface)

    
    uface = 0.0_r8
    scale = 0.0_r8
    do j = 1, mesh%ncell
      if (associated(skip)) then
        if (skip(j)) cycle
      end if
      associate (cface => mesh%cface(mesh%xcface(j):mesh%xcface(j+1)-1))
        uface(cface) = uface(cface) + ucell(j)
        scale(cface) = scale(cface) + 1
      end associate
    end do
    call mesh%face_imap%gather_offp(uface)
    call mesh%face_imap%gather_offp(scale)

    where (scale == 0)
      uface = 0.0_r8
    elsewhere
      uface = uface / scale
    end where

  end subroutine average_to_faces
  
  !! This auxillary procedure solves for the consistent face temperatures given
  !! cell temperatures (and concentrations -- heat equation parameters may
  !! depend on concentration).  It also computes the consistent radiosities if
  !! enclosure radiation problem exist.
  
  subroutine compute_face_temp (this, t, state, u, params)
  
    use mfd_diff_matrix_type
    use pcsr_precon_class
    use pcsr_precon_factory
    use parameter_list_type
    use nka_type
    use parallel_communication, only: global_maxval, global_all

    type(HTSD_model), intent(inout) :: this
    real(r8), intent(in) :: t, state(:,:) ! lower bound is now 1
    real(r8), intent(inout), target :: u(:)
    type(parameter_list) :: params
    
    integer :: n, stat, num_itr, max_iter, iter
    type(mfd_diff_matrix), target :: matrix
    real(r8) :: atol, rtol, error, r0_err, r_err, dT_max
    real(r8), pointer :: Tface(:), var(:), Fface(:)
    real(r8), allocatable :: z(:), udot(:), f(:)
    type(nka) :: nka_obj
    procedure(pardp), pointer :: dp
    target :: f
    character(80) :: string
    character(:), allocatable :: errmsg
    class(pcsr_precon), allocatable :: precon
    
    call TLS_info ('  computing consistent face temperatures and radiosities ...')
    
    !! Solve for the radiosity components given the (approx) face temperatures.
    call HTSD_model_get_face_temp_view (this, u, Tface)
    if (associated(this%ht%vf_rad_prob)) then
      do n = 1, size(this%ht%vf_rad_prob)
        call HTSD_model_get_radiosity_view (this, n, u, var)
        var = 0.0_r8
        associate (faces => this%ht%vf_rad_prob(n)%faces)
          call this%ht%vf_rad_prob(n)%solve_radiosity (t, Tface(faces), var, stat, num_itr, error)
        end associate
        if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
          write(string,'(2x,a,i0,a,es9.2," (",i0,")")') 'radiosity[', n, ']: |r|/|b|=', error, num_itr
          call TLS_info (string)
        end if
        if (stat /= 0) call TLS_info ('  WARNING: radiosities not converged')
      end do
    end if
    
    !! Initial residual of the face temperature equations.
    allocate(udot(size(u)), f(size(u)))
    udot = 0.0_r8
    call HTSD_model_compute_f (this, t, u, udot, f)
    call HTSD_model_get_face_temp_view (this, f, Fface)
    r0_err = sqrt(global_dot_product(Fface,Fface))
    if (r0_err == 0.0_r8) return
    if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
      write(string,'(2x,a,es9.2)') '||Rface(0)||_2 = ', r0_err
      call TLS_info (string)
    end if
    
    !! Form the Jacobian (approx) of the face temperature system.
    call matrix%init (this%disc)
    call make_matrix (matrix)
    
    call alloc_pcsr_precon (precon, matrix%a22, params, stat, errmsg)
    if (stat /= 0) call TLS_fatal('COMPUTE_FACE_TEMP: ' // errmsg)
    call precon%compute
    
    call params%get ('atol-temp', atol, default=0.0_r8)
    call params%get ('rtol-temp', rtol, default=0.0_r8)
    call params%get ('max-iter', max_iter)
    
    !! Simple Picard iteration.
    allocate(z(size(Fface)))
    call nka_obj%init (size(z), 10)
    dp => pardp
    call nka_obj%set_dot_prod (dp)
    
    iter = 0
    do ! until converged
    
      iter = iter + 1
      if (iter > max_iter) exit

      z = Fface
      call precon%apply (z)
      call nka_obj%accel_update (z) ! nonlinear krylov acceleration (NKA)
      
      Tface = Tface - z

      dT_max = global_maxval(abs(z))
      
      !! Solve the radiosity components given the new face temperatures.
      if (associated(this%ht%vf_rad_prob)) then
        do n = 1, size(this%ht%vf_rad_prob)
          call HTSD_model_get_radiosity_view (this, n, u, var)
          associate (faces => this%ht%vf_rad_prob(n)%faces)
            call this%ht%vf_rad_prob(n)%solve_radiosity (t, Tface(faces), var, stat, num_itr, error)
          end associate
          if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
            write(string,'(2x,a,i0,a,es9.2," (",i0,")")') 'radiosity[', n, ']: |r|/|b|=', error, num_itr
            call TLS_info (string)
          end if
          if (stat /= 0) call TLS_info ('  WARNING: radiosities not converged')
        end do
      end if
      
      !! Recompute the face temp residual.
      call HTSD_model_compute_f (this, t, u, udot, f)
      r_err = sqrt(global_dot_product(Fface,Fface))
      if (TLS_VERBOSITY >= TLS_VERB_NOISY) then
        write(string,'(2x,a,i0,2(a,es9.2))') &
            '||Rface(', iter, ')||=', r_err, ', ||ΔTface||_max=', dT_max
        call TLS_info (string)
      end if
      
      !! Check for convergence.
      if (global_all(abs(z) <= atol + rtol * abs(Tface))) exit
    end do
    
    write(string,'(2x,a,i0,3(a,es9.2))') &
        '||Rface(', iter, ')||=', r_err, ', ||Rface(0)||=', r0_err, &
        ', ||ΔTface||_max=', dT_max
    call TLS_info (string)
    
    if (iter > max_iter) call TLS_info ('  WARNING: face temperatures not converged')
    
  contains
  
    !! This auxillary routine computes an approximate Jacobian matrix for the
    !! time-independent heat equation (diffusion operator). This is a cell and
    !! face system matrix, but we are only interested in the (2,2) block (face-
    !! face).  That part is approximate in that it ignores the dependence of
    !! the radiosity on face temperatures (assuming the radiosities have been
    !! solved for). Moreover the Jacobian would be constant except for the T^4
    !! terms present in radiation BC.
    !!
    !! Note that this is nearly identical to code in HTSD_precon_type:compute
    !! that computes the preconditioner matrix.  CAN WE AVOID DUPLICATION?

    subroutine make_matrix (matrix)
    
      type(mfd_diff_matrix), intent(inout) :: matrix
      
      integer :: j, n, index
      integer, allocatable :: more_dir_faces(:)
      real(r8) :: D(this%mesh%ncell), Tface(this%mesh%nface)
      real(r8), allocatable :: values(:)
      
      !! Generate list of void faces; these are treated like Dirichlet BC.
      if (associated(this%void_face)) then
        n = count(this%void_face)
        allocate(more_dir_faces(n))
        n = 0
        do j = 1, this%mesh%nface
          if (this%void_face(j)) then
            n = n + 1
            more_dir_faces(n) = j
          end if
        end do
      else
        allocate(more_dir_faces(0))
      end if
      
      !! Compute the diffusion coefficient.
      call this%ht%conductivity%compute_value(state, D)
      if (associated(this%void_cell)) then
        where (this%void_cell) D = 0.0_r8
      end if
     
      !! Compute the basic diffusion matrix
      call matrix%compute (D)
      
      !! Dirichlet boundary condition fixups.
      if (allocated(this%ht%bc_dir)) then
        call this%ht%bc_dir%compute(t)
        call matrix%set_dir_faces(this%ht%bc_dir%index)
      end if

      call HTSD_model_get_face_temp_copy (this, u, Tface)
      call this%mesh%face_imap%gather_offp(Tface)

      !! External HTC boundary condition contribution.
      if (allocated(this%ht%bc_htc)) then
        call this%ht%bc_htc%compute(t, Tface)
        associate (index => this%ht%bc_htc%index, &
                   deriv => this%ht%bc_htc%deriv)
          call matrix%incr_face_diag(index, deriv)
        end associate
      end if
     
      !! Simple radiation boundary condition contribution.
      if (allocated(this%ht%bc_rad)) then
        call this%ht%bc_rad%compute(t, Tface)
        associate (index => this%ht%bc_rad%index, &
                   deriv => this%ht%bc_rad%deriv)
          call matrix%incr_face_diag(index, deriv)
        end associate
      end if

      !! Experimental evaporation heat flux contribution.
      if (allocated(this%ht%evap_flux)) then
        call this%ht%evap_flux%compute_deriv(t, Tface)
        associate (index => this%ht%evap_flux%index, &
                   deriv => this%ht%evap_flux%deriv)
          call matrix%incr_face_diag(index, this%mesh%area(index)*deriv)
        end associate
      endif

      !! Internal HTC interface condition contribution.
      if (allocated(this%ht%ic_htc)) then
        call this%ht%ic_htc%compute(t, Tface)
        associate (index => this%ht%ic_htc%index, &
                   deriv => this%ht%ic_htc%deriv)
          if (associated(this%void_face)) then
            do j = 1, size(index,2) !FIXME? Bad form to modify deriv?
              if (any(this%void_face(index(:,j)))) deriv(:,j) = 0.0_r8
            end do
          end if
          call matrix%incr_interface_flux3(index, deriv) !TODO: rename these methods
        end associate
      end if

      !! Internal gap radiation condition contribution.
      if (allocated(this%ht%ic_rad)) then
        call this%ht%ic_rad%compute(t, Tface)
        associate (index => this%ht%ic_rad%index, &
                   deriv => this%ht%ic_rad%deriv)
          if (associated(this%void_face)) then
            do j = 1, size(index,2) !FIXME? Bad form to modify deriv?
              if (any(this%void_face(index(:,j)))) deriv(:,j) = 0.0_r8
            end do
          end if
          call matrix%incr_interface_flux3(index, deriv) !TODO: rename these methods
        end associate
      end if

      !! Dirichlet fix-ups for void faces.
      call matrix%set_dir_faces (more_dir_faces)

      !! Enclosure radiation contributions to the preconditioner.
      !! This captures T^4 emitted heat (simple) but ignores absorbed heat
      !! from other faces (complex and non-local).
      if (associated(this%ht%vf_rad_prob)) then
        do index = 1, size(this%ht%vf_rad_prob)
          associate (faces => this%ht%vf_rad_prob(index)%faces)
            allocate(values(size(faces)))
            call this%ht%vf_rad_prob(index)%rhs_deriv (t, Tface(faces), values)
            where (.not.this%ht%vf_rad_prob(index)%fmask) values = 0
            call matrix%incr_face_diag (faces, this%mesh%area(faces) * values)
            deallocate(values)
          end associate
        end do
      end if

    end subroutine make_matrix

  end subroutine compute_face_temp
  
  function pardp (a, b) result (dp)
    use parallel_communication, only: global_dot_product
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: dp
    dp = global_dot_product(a, b)
  end function pardp

end module HTSD_init_cond_type

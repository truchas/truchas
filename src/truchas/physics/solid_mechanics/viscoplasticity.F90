!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define VERBOSE_BDF2 .false.

#include "f90_assert.fpp"

Module VISCOPLASTICITY
!-----------------------------------------------------------------------------
!Purpose:
!   This module has the routines for calculating elastic stresses, and elastic 
!   and plastic strains. 
!
!  Public interfaces:
!     Call PLASTIC_STRAIN_INCREMENT (Pl_Strain)
!     Call VISCOPLASTIC_STRAIN_RATE(Strain_Rate, Stress)
!     Call MTS_STRAIN_RATE(stress, temp, mm, strain_Rate)
!     Call MATERIAL_STRAINS (u)
!     Call DISPLACEMENT_GRADIENT (Q, dQ_dx, dQ_dy, dQ_dz)
!     Call MATERIAL_STRESSES (Strain, Stress)
!
!Authors: Dave Korzekwa (dak@lanl.gov)
!----------------------------------------------------------------------------- 
  use kinds, only: r8
  use time_step_module, only: dt
  use VP_model_class
  use solid_mechanics_mesh, only: ncomps
  use material_model_driver, only: matl_model
  use truchas_logging_services
  implicit none
  private

  ! Public procedures
  Public :: VISCOPLASTICITY_INIT,     &
            MATERIAL_STRESSES,        &
            MATERIAL_STRAINS,         &
            VISCOPLASTIC_STRAIN_RATE, &
            PLASTIC_STRAIN_INCREMENT, &
            DEVIATORIC_STRESS

  INTERFACE MATERIAL_STRESSES
     MODULE PROCEDURE MATERIAL_STRESSES_ONE_CELL
     MODULE PROCEDURE MATERIAL_STRESSES_ALL_CELLS
  END INTERFACE

  interface viscoplastic_strain_rate
    procedure viscoplastic_strain_rate_one
    procedure viscoplastic_strain_rate_all
  end interface

  ! Data needed by the strain rate routine for the bdf2 integrator,
  ! that cannot be passed to the bdf2 routine directly
  ! Cell number - needed for temperatue and material data
  integer :: cell_no, bdf2_unit
  integer, save :: bdf2_seq = 0
  ! Effective stress at the beginning of the time step
  real(r8) :: E_Stress_0
  real(r8), dimension(ncomps) :: D_Stress_0
  ! Total strain tensors for a single cell at the beginning and end of the time step
  real(r8), dimension(ncomps) :: T_Strain_0
  real(r8), dimension(ncomps) :: T_Strain_1
  ! Plastic strain tensor for a single cell at the beginning of the time step
  real(r8), dimension(ncomps) :: P_Strain_0

  ! The VP array holds the viscoplastic model for each of the materials.
  type :: VP_model_box
    class(VP_model), pointer :: model => null()
  end type
  type(VP_model_box), allocatable, save :: vp(:)

Contains

  subroutine viscoplasticity_init (plastic)

    use viscoplastic_model_namelist

    logical, intent(out) :: plastic

    integer :: m!, phase_id

    allocate(vp(matl_model%nphase))
    plastic = .false.
    do m = 1, matl_model%nphase_real
      vp(m)%model => get_VP_model(matl_model%phase_name(m))
      plastic = plastic .or. associated(vp(m)%model)
    end do

  end subroutine viscoplasticity_init

  !----------------------------------------------------------------------------
  !
  Subroutine PLASTIC_STRAIN_INCREMENT (Pl_Strain_Inc, Strain_Rate, Stress, LHS_Stress, &
                                  Dev_Stress_old, Stress_old, Strain_Rate_old, &
                                Eff_Stress_old, Tot_Strain_new, Tot_Strain_old, &
                                Plastic_Strain_old)
  !
  !----------------------------------------------------------------------------

    use legacy_mesh_api,       only: ncells, Mesh, GAP_ELEMENT_1
    use solid_mechanics_data,  only: Thermal_Strain, PC_Strain, plasticity
    use solid_mechanics_input, only: strain_limit
    Use zone_module,           only: zone
    Use time_step_module,      Only: cycle_number
    use bdf2_kinds
    use bdf2_integrator

    !Arguments
    real(r8), Dimension(ncomps,ncells), Intent(OUT)   :: Pl_Strain_Inc
    real(r8), Dimension(ncells),        Intent(OUT)   :: Strain_Rate
    real(r8), Dimension(ncomps,ncells), Intent(OUT)   :: Stress
    real(r8), Dimension(ncomps,ncells), Intent(OUT)   :: LHS_Stress

    real(r8), Dimension(ncomps,ncells), Intent(IN)    :: Dev_Stress_old
    real(r8), Dimension(ncomps,ncells), Intent(IN)    :: Stress_old
    real(r8), Dimension(ncells),        Intent(IN)    :: Strain_Rate_old
    real(r8), Dimension(ncells),        Intent(IN)    :: Eff_Stress_old
    real(r8), Dimension(ncomps,ncells), Intent(IN)    :: Tot_Strain_new
    real(r8), Dimension(ncomps,ncells), Intent(IN)    :: Tot_Strain_old
    real(r8), Dimension(ncomps,ncells), Intent(IN)    :: Plastic_Strain_old

    ! Local Variables
    real(r8), Dimension(ncells)         :: Eff_Stress
    real(r8), Dimension(ncomps,ncells)  :: Strain
    real(r8), Dimension(ncomps,ncells)  :: Dev_Stress
    integer  :: icomp, icell
    real(r8) :: max_strain, rate_change
    integer  :: nplas, nbdf

    ! bdf2 integrator
    type(bdf2_control) :: control
    type(bdf2_state)   :: state
    integer :: stat
    real(kind=rk) :: tout, t, atol(ncomps), depsilon(ncomps), hstep
    logical :: bdf2_failed, verbose

    ! Initialize output
    Pl_Strain_Inc = 0.0

    nplas = 0
    nbdf = 0
       
    ! If it is the initial elastic calculation or all materials are elastic, 
    ! do not calculate the plastic strain - this needs to be handled differently
    ! when we add more plasticity models
    if (plasticity .and. (cycle_number /= 0)) then

       bdf2_seq = bdf2_seq + 1
       bdf2_failed = .false.
       call bdf2_create_state (state, size(depsilon))
       
       do icell = 1, ncells

          if ((Strain_Rate_old(icell) > 1.0e-12) .and. &
               (Mesh(icell)%Cell_Shape < GAP_ELEMENT_1)) then

             nplas = nplas + 1

             ! Compute strain rate at midpoint of the time step using initial strain rate
             do icomp = 1,ncomps
                if (Strain_Rate_old(icell) > 1.0e-12) then
                   Pl_Strain_Inc(icomp,icell) = 0.5 * dt * 1.5 * Dev_Stress_old(icomp,icell) * &
                        Strain_Rate_old(icell)/Eff_Stress_old(icell)
                end if
             end do
                
             ! Elastic strain and elastic stress at midpoint - using the (questionable?)
             ! assumption that elastic strain components at the beginning and end of the time step 
             ! can be averaged.  We need to account for the thermal strain here.
             Strain(:,icell) = 0.5 * (Tot_Strain_new(:,icell) + Tot_Strain_old(:,icell)) - &
                               Plastic_Strain_old(:,icell) - Pl_Strain_Inc(:,icell) - &
                               Thermal_Strain(:,icell) - PC_Strain(:,icell)

             Call MATERIAL_STRESSES (Strain(:,icell), Stress(:,icell),icell)
             
             ! Effective stress, deviatoric stress and strain rate at the midpoint
             Eff_Stress(icell) = sqrt(((Stress(1,icell) - Stress(2,icell))**2 + &
                  (Stress(2,icell) - Stress(3,icell))**2 + &
                  (Stress(3,icell) - Stress(1,icell))**2 + &
                  6.0 * (Stress(4,icell)**2 + Stress(5,icell)**2 + Stress(6,icell)**2))/2.)

             Call DEVIATORIC_STRESS(Stress(:,icell), Dev_Stress(:,icell))

             Call VISCOPLASTIC_STRAIN_RATE(icell, Eff_Stress(icell), Zone(icell)%temp, Strain_Rate(icell))

             rate_change = Strain_Rate(icell)/ Strain_Rate_old(icell)
             if (rate_change == 0.0) then
                rate_change = 1.0e6
             else if (rate_change < 1.0) then
                rate_change = 1.0/rate_change
             end if
             max_strain = MAX(Strain_Rate(icell), Strain_Rate_old(icell)) * dt
             if ((rate_change < 1.1) .or. (max_strain < strain_limit)) then

                ! Strain increment for the whole time step using the strain rate and deviatoric stress
                ! from the midpoint projection
                do icomp = 1,ncomps
                   Pl_Strain_Inc(icomp,icell) = dt * 1.5 * Dev_Stress(icomp,icell) * &
                                                Strain_Rate(icell)/Eff_Stress(icell)
                end do

             else

                nbdf = nbdf + 1

                Pl_Strain_Inc(:,icell) = 0.0
                Eff_Stress(icell) = 0.0
                Strain(:,icell) = 0.0
                Stress(:,icell) = 0.0
                Dev_Stress(:,icell) = 0.0

                ! This should be set to a number on the order of the smallest strain increment of interest
                atol = 1.0e-12

                call bdf2_set_param (control, atol=atol, rtol=1.0d-3, mvec=1, ntol=0.01d0)

                t = 0.0
                depsilon = 0.0

                ! Pick hstart based on time step, strain limit input parameter and initial strain rate
                hstep = strain_limit / Strain_Rate_old(icell) / 10.0
                if (hstep > dt/10.0) hstep = dt/10.0 ! This shouldn't happen often

                tout = dt
                D_stress_0 = Dev_Stress_old(:,icell)
                E_Stress_0 = Eff_Stress_old(icell)
                T_Strain_0 = Tot_Strain_old(:,icell)
                T_Strain_1 = Tot_Strain_new(:,icell)
                P_Strain_0 = Plastic_Strain_old(:,icell)
                cell_no    = icell

                verbose = VERBOSE_BDF2
                if (.not.verbose) then
                   call bdf2_init_state (state, depsilon, t, hstart=hstep)
                   call bdf2_integrate (state, control, tout=tout, stat=stat, rhs=dstrain_dt)
                   select case (stat)
                   case (SOLN_AT_TOUT)
                      depsilon = bdf2_interpolate_solution (state, tout)
                   case default
                      bdf2_failed = .true.
                      verbose = .true. ! force a reintegration to collect more diagnostics
                   end select
                end if

                if (verbose) then
                   call create_bdf2_diagnostics_file (bdf2_unit)
                   call bdf2_init_state (state, depsilon, t, hstart=hstep, profile=.true.)
                   call bdf2_integrate (state, control, tout=tout, stat=stat, rhs=dstrain_dt, user=userout)
                   select case (stat)
                   case (SOLN_AT_TOUT)
                      depsilon = bdf2_interpolate_solution (state, tout)
                      call userout (tout, depsilon)
                   case default
                      bdf2_failed = .true.
                   end select
                   call bdf2_write_profile (state, bdf2_unit)
                   close(bdf2_unit)
                end if

                if (bdf2_failed) exit

                do icomp = 1,ncomps
                   Pl_Strain_Inc(icomp,icell) = depsilon(icomp)
                end do

             end if
          end if
       end do
       
       call destroy (state)
       call destroy (control)

       call TLS_fatal_if_any (bdf2_failed, 'PLASTIC_STRAIN_INCREMENT: BDF2_INTEGRATE failure: see bdf2-*.err files for details')

    end if

    ! Get elastic stress and plastic strain rate at the end of the step
    Strain = Tot_Strain_new - Plastic_Strain_old - Pl_Strain_Inc - Thermal_Strain - PC_Strain
    Call MATERIAL_STRESSES (Strain, Stress)

    Eff_Stress(:) = sqrt(((Stress(1,:) - Stress(2,:))**2 + &
         (Stress(2,:) - Stress(3,:))**2 + &
         (Stress(3,:) - Stress(1,:))**2 + &
         6.0 * (Stress(4,:)**2 + Stress(5,:)**2 + Stress(6,:)**2))/2.)

    call viscoplastic_strain_rate (eff_stress, zone%temp, strain_rate)

    ! Finally, we need to remove the thermal stress before returning the LHS_stress values to 
    ! the residual calculation (thermal stresses are already included in the RHS).
    Strain = Strain + Thermal_Strain + PC_Strain
    Call MATERIAL_STRESSES (Strain, LHS_Stress)

!    write(*,*) nplas, nbdf

  contains

    subroutine create_bdf2_diagnostics_file (unit)

      use string_utilities, only: i_to_c
      use parallel_info_module, only: p_info
      use legacy_mesh_api, only: unpermute_mesh_vector
      use truchas_env, only: output_dir

      integer, intent(out) :: unit

      open(newunit=unit,file=trim(output_dir)//'bdf2-'//i_to_c(p_info%thisPE)//'.'// & 
           i_to_c(unpermute_mesh_vector(cell_no))//'.'//i_to_c(bdf2_seq)//'.err', &
           status='replace', position='rewind', action='write')

      write(unit,fmt='(a,i4)') 'PE=', p_info%thisPE
      write(unit,fmt='(a,i7)') 'LOC_CELL=', cell_no
      write(unit,fmt='(a,i7)') 'USR_CELL=', unpermute_mesh_vector(cell_no)
      write(unit,fmt='(a,i3)') 'BDF2_STAT=', stat

      write(unit,fmt='(/,a,es12.4)') 'E_Stress_0=', E_Stress_0
      write(unit,fmt='(a,6es12.4)') 'D_stress_0=', D_stress_0
      write(unit,fmt='(a,6es12.4)') 'T_Strain_0=', T_Strain_0
      write(unit,fmt='(a,6es12.4)') 'T_Strain_1=', T_Strain_1
      write(unit,fmt='(a,6es12.4)') 'P_Strain_0=', P_Strain_0

      write(unit,fmt='(/,a,es12.4)') 'HSTART=', hstep

      write(unit,fmt='(/,a)') 'TIME DEPSILON DEPSILON_DOT'

    end subroutine create_bdf2_diagnostics_file

  end Subroutine PLASTIC_STRAIN_INCREMENT

  Subroutine USEROUT (t, epsilon)
    use bdf2_kinds
    real(kind=rk), intent(in) :: t, epsilon(:)
    real(kind=rk) :: epsilon_dot(size(epsilon))
    call dstrain_dt (t, epsilon, epsilon_dot)
    write(bdf2_unit, fmt='(/,1es12.4)') t
    write(bdf2_unit, fmt='(6es12.4)') epsilon
    write(bdf2_unit, fmt='(6es12.4)') epsilon_dot
  end Subroutine USEROUT

  Subroutine DSTRAIN_DT(t, depsilon, depsdt)
    !
    ! Purpose:  Calculate the time derivative of the plastic strain for the bdf2
    ! ODE integrator.
    !
    !Authors: Dave Korzekwa (dak@lanl.gov)
    !-----------------------------------------------------------------------------  

    use bdf2_kinds
    use solid_mechanics_data, only: Thermal_Strain, PC_Strain, Thermal_Strain_Inc
    use zone_module, only: Zone
    
    real(kind=rk), intent(in)  :: t
    real(kind=rk), intent(in), dimension(:)  :: depsilon
    real(kind=rk), intent(out), dimension(:) :: depsdt

    ! Local variables
    integer :: icomp
    real(r8), dimension(ncomps) :: T_Strain_t
    real(r8), dimension(ncomps) :: Stress_t
    real(r8), dimension(ncomps) :: D_Stress_t
    real(r8) :: E_Stress_t
    real(r8) :: epsdot_eff

    ! Elastic strain and elastic stress at the current time - using the assumption 
    ! that elastic strain components at the beginning and end of the time step 
    ! can be averaged.  We need to account for the thermal and phase change strains here.
    do icomp = 1,ncomps
       T_Strain_t(icomp) = (dt - t)/dt * T_Strain_0(icomp) + t/dt * T_Strain_1(icomp) - &
                           P_Strain_0(icomp) - depsilon(icomp) - Thermal_Strain(icomp,cell_no) + &
                           (dt - t) * Thermal_Strain_Inc(icomp,cell_no) - PC_Strain(icomp,cell_no)
    end do

    Call MATERIAL_STRESSES (T_Strain_t, Stress_t, cell_no)

    E_Stress_t = sqrt(((Stress_t(1) - Stress_t(2))**2 + &
                 (Stress_t(2) - Stress_t(3))**2 + &
                 (Stress_t(3) - Stress_t(1))**2 + &
                 6.0 * (Stress_t(4)**2 + Stress_t(5)**2 + Stress_t(6)**2))/2.)

    Call VISCOPLASTIC_STRAIN_RATE(cell_no, E_Stress_t, Zone(cell_no)%temp, epsdot_eff)

    Call DEVIATORIC_STRESS(Stress_t, D_Stress_t)

    do icomp = 1,ncomps
       depsdt(icomp) =  epsdot_eff * 1.5 * D_Stress_t(icomp) / E_Stress_t
    end do
   
  END Subroutine DSTRAIN_DT

  !!
  !! Calculate the average plastic strain rate on a single cell.
  !! Volume fraction weighted average over cell materials.
  !!

  subroutine viscoplastic_strain_rate_one (icell, stress, temp, strain_rate)

    use parameter_module, only: mat_slot
    use matl_module, only: matl
    use time_step_module, Only: dt

    integer,  intent(in)  :: icell
    real(r8), intent(in)  :: stress
    real(r8), intent(in)  :: temp
    real(r8), intent(out) :: strain_rate

    integer :: s, m
    real(r8) :: rate_mat

    strain_rate = 0.0_r8
    do s = 1, mat_slot
      m = matl(s)%cell(icell)%id
      if (m == 0) cycle ! TODO: can we exit?  Are the unused slots always at the end?
      if (.not.associated(vp(m)%model)) cycle
      rate_mat = vp(m)%model%strain_rate(stress, temp, dt)
      strain_rate = strain_rate + matl(s)%cell(icell)%vof * rate_mat
    end do

  end subroutine viscoplastic_strain_rate_one

  !!
  !! Calculate the average plastic strain rate for all cells.
  !! Volume fraction weighted average over cell materials.
  !!

  subroutine viscoplastic_strain_rate_all (stress, temp, strain_rate)

    use legacy_mesh_api, only: ncells
    use matl_module, only: gather_vof
    use time_step_module, only: dt

    real(r8), intent(in)  :: stress(:)
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: strain_rate(:)

    integer  :: m, j
    real(r8) :: mrate, vofm(ncells)

    ASSERT(size(stress) == ncells)
    ASSERT(size(temp) == ncells)
    ASSERT(size(strain_rate) == ncells)

    strain_rate = 0.0_r8
    do m = 1, matl_model%nphase_real
      if (matl_model%is_fluid(m)) cycle
      if (.not.associated(vp(m)%model)) cycle
      call gather_vof (m, vofm)
      do j = 1, ncells
        if (vofm(j) > 0.0_r8) then
          mrate = vp(m)%model%strain_rate(stress(j), temp(j), dt)
          strain_rate(j) = strain_rate(j) + mrate * vofm(j)
        end if
      end do
    end do

  end subroutine viscoplastic_strain_rate_all

  !
  !-----------------------------------------------------------------------------
  !
  Subroutine MATERIAL_STRAINS (u, Strain, Rotation)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    Calculate the cell-centered solid material total strain field
    !
    !---------------------------------------------------------------------------

    use legacy_mesh_api, only: nnodes, ncells
    use solid_mechanics_mesh, only: ndim

    ! Argument list
    real(r8), Dimension(ndim*nnodes),   Intent(IN)  :: u
    real(r8), Dimension(ncomps,ncells), Intent(OUT)  :: Strain
    real(r8), Dimension(ncells), Intent(OUT), Optional  :: Rotation

    ! Local variables
    integer                              :: inodes
    real(r8), Dimension(nnodes)           :: Q
    real(r8), Dimension(ncells)           :: dQ_dx, dQ_dy, dQ_dz
    real(r8), Dimension(ndim,ndim,ncells) :: Grad_u
    real(r8), Allocatable, Dimension(:)   :: R1, R2, R3

    !---------------------------------------------------------------------------

    ! Calculate the cell-centered displacement gradients
    Do inodes = 1, nnodes
     Q(inodes) = u((inodes - 1) * ndim + 1)
    End Do

    Call DISPLACEMENT_GRADIENT (Q, dQ_dx, dQ_dy, dQ_dz)

    Grad_u(1,1,:) = dQ_dx(:)
    Grad_u(2,1,:) = dQ_dy(:)
    Grad_u(3,1,:) = dQ_dz(:)

    Do inodes = 1, nnodes
     Q(inodes) = u((inodes - 1) * ndim + 2)
    End Do

    Call DISPLACEMENT_GRADIENT (Q, dQ_dx, dQ_dy, dQ_dz)

    Grad_u(1,2,:) = dQ_dx(:)
    Grad_u(2,2,:) = dQ_dy(:)
    Grad_u(3,2,:) = dQ_dz(:)

    Do inodes = 1, nnodes
     Q(inodes) = u((inodes - 1) * ndim + 3)
    End Do

    Call DISPLACEMENT_GRADIENT (Q, dQ_dx, dQ_dy, dQ_dz)

    Grad_u(1,3,:) = dQ_dx(:)
    Grad_u(2,3,:) = dQ_dy(:)
    Grad_u(3,3,:) = dQ_dz(:)

    ! Calculate the cell-centered total strain field
    Strain(1,:) = Grad_u(1,1,:)
    Strain(2,:) = Grad_u(2,2,:)
    Strain(3,:) = Grad_u(3,3,:)
    Strain(4,:) = 0.5*(Grad_u(1,2,:) + Grad_u(2,1,:))
    Strain(5,:) = 0.5*(Grad_u(1,3,:) + Grad_u(3,1,:))
    Strain(6,:) = 0.5*(Grad_u(2,3,:) + Grad_u(3,2,:))

    ! Optionally calculate the magnitude of the rotation vector
    if (PRESENT(Rotation)) then
       Allocate (R1(ncells))
       Allocate (R2(ncells))
       Allocate (R3(ncells))
       R1(:) = 0.5*(Grad_u(1,2,:) - Grad_u(2,1,:))
       R2(:) = 0.5*(Grad_u(1,3,:) - Grad_u(3,1,:))
       R3(:) = 0.5*(Grad_u(2,3,:) - Grad_u(3,2,:))
       Rotation(:) = sqrt(R1(:)*R1(:) + R2(:)*R2(:) + R3(:)*R3(:))
       Deallocate(R1, R2, R3)
    end if

  End Subroutine MATERIAL_STRAINS
  !
  !-----------------------------------------------------------------------------
  !
  Subroutine DISPLACEMENT_GRADIENT (Q, dQ_dx, dQ_dy, dQ_dz)
    !---------------------------------------------------------------------------
    ! Purpose:
    ! 
    !   Calculate the cell-centered gradient {dQ_dx, dQ_dy, dQ_dz} of the
    !   vertex-centered displacement component Q by averaging the gradient, 
    !   defined by the solution of a 3x3 system of equations, over the volume.  
    !   This averaging procedure is accomplished through an integral over the 
    !   cell volume.  The gradient at any logical point within the cell is a 
    !   ratio of Jacobian-like determinants.
    !
    !---------------------------------------------------------------------------
    use legacy_mesh_api, only: nnodes, ncells, Cell, EN_GATHER, gather_vertex_coord
    use solid_mechanics_mesh, only: nvc

    ! Argument list
    real(r8), Dimension(nnodes), Intent(IN)  :: Q
    real(r8), Dimension(ncells), Intent(OUT) :: dQ_dx
    real(r8), Dimension(ncells), Intent(OUT) :: dQ_dy
    real(r8), Dimension(ncells), Intent(OUT) :: dQ_dz

    ! Local Variables
    real(r8), Dimension(ncells)     :: Q_f
    real(r8), Dimension(nvc,ncells) :: Q_v
    real(r8), Dimension(nvc,ncells) :: X_v
    real(r8), Dimension(nvc,ncells) :: Y_v
    real(r8), Dimension(nvc,ncells) :: Z_v
    !
    !---------------------------------------------------------------------------
    !
    ! Gather the vertex-centered scalar Q into the cell-vector quantity Q_v
    Call EN_GATHER (Q_v, Q)

    ! Gather the vertex coordinates
    call gather_vertex_coord (X_v, dim=1)
    call gather_vertex_coord (Y_v, dim=2)
    call gather_vertex_coord (Z_v, dim=3)

    ! Calculate the volume-averaged gradient
    Call DETERMINANT_VOL_AVG (Q_v, Y_v, Z_v, dQ_dx)
    Call DETERMINANT_VOL_AVG (X_v, Q_v, Z_v, dQ_dy)
    Call DETERMINANT_VOL_AVG (X_v, Y_v, Q_v, dQ_dz)

    ! Normalize the volume-averaged gradient by the cell volume
    Q_f(:)   = 1.0_r8/Cell(:)%Volume
    dQ_dx(:) = dQ_dx(:)*Q_f(:)
    dQ_dy(:) = dQ_dy(:)*Q_f(:)
    dQ_dz(:) = dQ_dz(:)*Q_f(:)

  End Subroutine DISPLACEMENT_GRADIENT



  SUBROUTINE DETERMINANT_VOL_AVG (Xv, Yv, Zv, Avg)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute the volume-averaged value of a Jacobian determinant |J|:
    !   Avg = Integral(|J|dV), where J is given by column vectors X, Y,
    !   and Z:             -                       -
    !                      |  X_xi   Y_xi   Z_xi   |
    !                J =   |  X_eta  Y_eta  Z_eta  |
    !                      |  X_zeta Y_zeta Z_zeta |
    !                      -                       -
    !   where X_xi == dX/dxi, X_eta == dX/deta, X_zeta == dX/dzeta, etc.
    !   The volume integral is converted to a surface integral and
    !   evaluated following J. Dukowicz, JCP 74: 493-496 (1988).
    !
    !=======================================================================
    use legacy_mesh_api, only: ncells, nfc, nvc

    ! Arguments
    real(r8), dimension(nvc,ncells), intent(IN)  :: Xv, Yv, Zv
    real(r8), dimension(ncells),     intent(OUT) :: Avg

    ! Local Variables
    integer :: f
    integer :: v1, v2, v3, v4, v5, v6

    real(r8), dimension(ncells) :: X1, Y1, Z1
    real(r8), dimension(ncells) :: X2, Y2, Z2
    real(r8), dimension(ncells) :: X3, Y3, Z3

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize relevant quantities
    Avg = 0.0_r8

    ! Loop over faces, accumulating the average
    do f = 1,nfc
       select case (f)
          case (1)
          ! Side 1 (vertices 4-8-7-3)
          v1 = 8; v2 = 4; v3 = 7; v4 = 8; v5 = 3; v6 = 4
          case (2)
          ! Side 2 (vertices 1-2-6-5)
          v1 = 6; v2 = 2; v3 = 5; v4 = 6; v5 = 1; v6 = 2
          case (3)
          ! Side 3 (vertices 4-1-5-8)
          v1 = 5; v2 = 1; v3 = 8; v4 = 5; v5 = 4; v6 = 1
          case (4)
          ! Side 4 (vertices 3-7-6-2)
          v1 = 7; v2 = 3; v3 = 6; v4 = 7; v5 = 2; v6 = 3
          case (5)
          ! Side 5 (vertices 4-3-2-1)
          v1 = 3; v2 = 4; v3 = 2; v4 = 3; v5 = 1; v6 = 4
          case (6)
          ! Side 6 (vertices 8-5-6-7)
          v1 = 6; v2 = 5; v3 = 7; v4 = 6; v5 = 8; v6 = 5
       end select

       X1 = Xv(v1,:) + Xv(v2,:)
       Y1 = Yv(v1,:) + Yv(v2,:)
       Z1 = Zv(v1,:) + Zv(v2,:)

       X2 = Xv(v3,:) + Xv(v4,:)
       Y2 = Yv(v3,:) + Yv(v4,:)
       Z2 = Zv(v3,:) + Zv(v4,:)

       X3 = Xv(v5,:) + Xv(v6,:)
       Y3 = Yv(v5,:) + Yv(v6,:)
       Z3 = Zv(v5,:) + Zv(v6,:)

       Avg = Avg + X1*(Y2*Z3 - Y3*Z2) + Y1*(X3*Z2 - X2*Z3) + Z1*(X2*Y3 - X3*Y2)
    end do

    Avg = Avg / 12.0_r8

  END SUBROUTINE DETERMINANT_VOL_AVG

  !
  !-----------------------------------------------------------------------------
  !
  Subroutine MATERIAL_STRESSES_ALL_CELLS (Strain, Stress)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    Calculate the cell-centered solid material stress field - all cells
    !
    !---------------------------------------------------------------------------
    use legacy_mesh_api, only: ncells, ndim, Mesh, GAP_ELEMENT_1
    use solid_mechanics_data, only: Lame1, Lame2
    use solid_mechanics_mesh, only: ncomps

    ! Argument list
    real(r8), Dimension(ncomps,ncells), Intent(IN) :: Strain
    real(r8), Dimension(ncomps,ncells), Intent(out):: Stress
    !
    ! Local variables
    integer                    :: icomps
    real(r8), Dimension(ncells) :: Dilatation
    !
    !---------------------------------------------------------------------------
    !
    ! Compute the dilatation
    Dilatation(:) = 0.0_r8

    Do icomps = 1, ndim
     Dilatation(:) = Dilatation(:) + Strain(icomps,:)
    End Do

    ! Calculate the cell-centered stresses
    Do icomps = 1, ndim
     Stress(icomps,:) = Lame1(:)*Dilatation(:) + 2.0*Lame2(:)*Strain(icomps,:)
    End Do

    Do icomps = ndim + 1, ncomps
     Stress(icomps,:) =                          2.0*Lame2(:)*Strain(icomps,:)
    End Do

    ! Return zero  stresses for gap elements
    Do icomps = 1, ncomps
       where (Mesh(:)%Cell_Shape >= GAP_ELEMENT_1) Stress(icomps,:) = 0.0
    End Do

  End Subroutine MATERIAL_STRESSES_ALL_CELLS
  !-----------------------------------------------------------------------------
  !
  Subroutine MATERIAL_STRESSES_ONE_CELL (Strain, Stress, icell)
    !---------------------------------------------------------------------------
    ! Purpose:
    !
    !    Calculate the cell-centered solid material stress field - all cells
    !
    !---------------------------------------------------------------------------
    use solid_mechanics_data, only: Lame1, Lame2
    use solid_mechanics_mesh, only: ndim, ncomps
    use legacy_mesh_api,      only: Mesh, GAP_ELEMENT_1

    ! Argument list
    real(r8), Dimension(ncomps), Intent(IN) :: Strain
    real(r8), Dimension(ncomps), Intent(out):: Stress
    integer,Intent(IN)                  :: icell
    !
    ! Local variables
    integer                    :: icomps
    real(r8)                    :: Dilatation
    !
    !---------------------------------------------------------------------------
    !
    ! Return zero  stresses for gap elements

    if (Mesh(icell)%Cell_Shape >= GAP_ELEMENT_1) then 
       Stress = 0.0
    else
       ! Compute the dilatation
       Dilatation = 0.0_r8

       Do icomps = 1, ndim
          Dilatation = Dilatation + Strain(icomps)
       End Do

       ! Calculate the cell-centered stresses
       Do icomps = 1, ndim
          Stress(icomps) = Lame1(icell)*Dilatation + 2.0*Lame2(icell)*Strain(icomps)
       End Do

       Do icomps = ndim + 1, ncomps
          Stress(icomps) =                           2.0*Lame2(icell)*Strain(icomps)
       End Do
    end if

  End Subroutine MATERIAL_STRESSES_ONE_CELL

  Subroutine DEVIATORIC_STRESS(Stress, Dev_Stress)
    ! Calculate the deviatoric stress tensor from the elastic stress tensor
    use solid_mechanics_mesh, only: ndim, ncomps

    ! Argument list
    real(r8), Dimension(ncomps), Intent(IN) :: Stress
    real(r8), Dimension(ncomps), Intent(out):: Dev_Stress
    !
    ! Local variables
    integer                    :: idim
    real(r8)                    :: Mean_Stress

    Mean_Stress = (Stress(1) + Stress(2) + Stress(3))/3.0
    do idim = 1,ndim
       Dev_Stress(idim) = Stress(idim) - Mean_Stress
    end do
    do idim = ndim+1,ncomps
       Dev_Stress(idim) = Stress(idim)
    end do

  end Subroutine DEVIATORIC_STRESS

end Module VISCOPLASTICITY

#define VERBOSE_BDF2 .false.

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
  Use parameter_module, only: ncomps
  use time_step_module, only: dt
  Use kind_module,      only: int_kind, real_kind
  use truchas_logging_services

  Private

  ! Public procedures
  Public :: MATERIAL_STRESSES,        &
            MATERIAL_STRAINS,         &
            VISCOPLASTIC_STRAIN_RATE, &
            PLASTIC_STRAIN_INCREMENT, &
            DEVIATORIC_STRESS

  INTERFACE MATERIAL_STRESSES
     MODULE PROCEDURE MATERIAL_STRESSES_ONE_CELL
     MODULE PROCEDURE MATERIAL_STRESSES_ALL_CELLS
  END INTERFACE


  ! Data needed by the strain rate routine for the bdf2 integrator,
  ! that cannot be passed to the bdf2 routine directly
  ! Cell number - needed for temperatue and material data
  integer(KIND = int_kind)                   :: cell_no, bdf2_unit
  integer(KIND = int_kind), save             :: bdf2_seq = 0
  ! Effective stress at the beginning of the time step
  Real(KIND = real_kind)                     :: E_Stress_0
  Real(KIND = real_kind),dimension(ncomps)   :: D_Stress_0
  ! Total strain tensors for a single cell at the beginning and end of the time step
  Real(KIND = real_kind),dimension(ncomps)   :: T_Strain_0
  Real(KIND = real_kind),dimension(ncomps)   :: T_Strain_1
  ! Plastic strain tensor for a single cell at the beginning of the time step
  Real(KIND = real_kind),dimension(ncomps)   :: P_Strain_0

Contains

  !----------------------------------------------------------------------------
  !
  Subroutine PLASTIC_STRAIN_INCREMENT (Pl_Strain_Inc, Strain_Rate, Stress, LHS_Stress, &
                                  Dev_Stress_old, Stress_old, Strain_Rate_old, &
                                Eff_Stress_old, Tot_Strain_new, Tot_Strain_old, &
                                Plastic_Strain_old)
  !
  !----------------------------------------------------------------------------

    use parameter_module,     only: ncells
    use solid_mechanics_data, only: Thermal_Strain, PC_Strain, plasticity, strain_limit
    Use time_step_module,     Only: cycle_number
    use bdf2_kinds
    use bdf2_integrator
    use mesh_module,          only: Mesh, GAP_ELEMENT_1

    !Arguments
    Real(KIND = real_kind), Dimension(ncomps,ncells), Intent(OUT)   :: Pl_Strain_Inc
    Real(KIND = real_kind), Dimension(ncells),        Intent(OUT)   :: Strain_Rate
    Real(KIND = real_kind), Dimension(ncomps,ncells), Intent(OUT)   :: Stress
    Real(KIND = real_kind), Dimension(ncomps,ncells), Intent(OUT)   :: LHS_Stress

    Real(KIND = real_kind), Dimension(ncomps,ncells), Intent(IN)    :: Dev_Stress_old
    Real(KIND = real_kind), Dimension(ncomps,ncells), Intent(IN)    :: Stress_old
    Real(KIND = real_kind), Dimension(ncells),        Intent(IN)    :: Strain_Rate_old
    Real(KIND = real_kind), Dimension(ncells),        Intent(IN)    :: Eff_Stress_old
    Real(KIND = real_kind), Dimension(ncomps,ncells), Intent(IN)    :: Tot_Strain_new
    Real(KIND = real_kind), Dimension(ncomps,ncells), Intent(IN)    :: Tot_Strain_old
    Real(KIND = real_kind), Dimension(ncomps,ncells), Intent(IN)    :: Plastic_Strain_old

    ! Local Variables
    Real(KIND = real_kind), Dimension(ncells)         :: Eff_Stress
    Real(KIND = real_kind), Dimension(ncomps,ncells)  :: Strain
    Real(KIND = real_kind), Dimension(ncomps,ncells)  :: Dev_Stress
    integer(KIND = int_kind)                          :: icomp, icell
    Real(KIND = real_kind)                            :: max_strain, rate_change
    integer(KIND = int_kind)                          :: nplas, nbdf

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

             Call VISCOPLASTIC_STRAIN_RATE(Strain_Rate(icell), Eff_Stress(icell), icell)

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

    do icell = 1, ncells
       Call VISCOPLASTIC_STRAIN_RATE(Strain_Rate(icell), Eff_Stress(icell), icell)
    end do

    ! Finally, we need to remove the thermal stress before returning the LHS_stress values to 
    ! the residual calculation (thermal stresses are already included in the RHS).
    Strain = Strain + Thermal_Strain + PC_Strain
    Call MATERIAL_STRESSES (Strain, LHS_Stress)

!    write(*,*) nplas, nbdf

  contains

    subroutine create_bdf2_diagnostics_file (unit)

      use string_utilities, only: i_to_c
      use parallel_info_module, only: p_info
      use mesh_module, only: unpermute_mesh_vector
#ifdef SUPPORTS_NEWUNIT
      use truchas_env, only: output_dir
#else
      use truchas_env, only: output_dir, new_unit
#endif

      integer, intent(out) :: unit

#ifdef SUPPORTS_NEWUNIT
      open(newunit=unit,file=trim(output_dir)//'bdf2-'//i_to_c(p_info%thisPE)//'.'// & 
           i_to_c(unpermute_mesh_vector(cell_no))//'.'//i_to_c(bdf2_seq)//'.err', &
           status='replace', position='rewind', action='write')
#else
      call new_unit (unit)
      open(unit, status='replace', position='rewind', action='write', &
           file=trim(output_dir)//'bdf2-'//i_to_c(p_info%thisPE)//'.'// & 
           i_to_c(unpermute_mesh_vector(cell_no))//'.'//i_to_c(bdf2_seq)//'.err')
#endif

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
    use solid_mechanics_data, only: Thermal_Strain, PC_Strain, &
                                    Thermal_Strain_Inc
    
    real(kind=rk), intent(in)  :: t
    real(kind=rk), intent(in), dimension(:)  :: depsilon
    real(kind=rk), intent(out), dimension(:) :: depsdt

    ! Local variables
    integer(kind=int_kind)                  :: icomp
    real(kind=real_kind), dimension(ncomps) :: T_Strain_t
    real(kind=real_kind), dimension(ncomps) :: Stress_t
    real(kind=real_kind), dimension(ncomps) :: D_Stress_t
    real(kind=real_kind)                    :: E_Stress_t
    real(kind=real_kind)                    :: epsdot_eff

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

    Call VISCOPLASTIC_STRAIN_RATE(epsdot_eff, E_Stress_t, cell_no)

    Call DEVIATORIC_STRESS(Stress_t, D_Stress_t)

    do icomp = 1,ncomps
       depsdt(icomp) =  epsdot_eff * 1.5 * D_Stress_t(icomp) / E_Stress_t
    end do
   
  END Subroutine DSTRAIN_DT

  !
  !--------------------------------------------------------------------------
  !
  Subroutine VISCOPLASTIC_STRAIN_RATE(Strain_Rate, Stress, icell)
  !
  !--------------------------------------------------------------------------
  !
  ! Purpose:
  !   Calculate the plastic strain rate in a cell averaged by volume fraction.
  !
  !--------------------------------------------------------------------------
    use fluid_data_module,    only: isImmobile
    use parameter_module,     only: nmat, mat_slot
    Use matl_module,          Only: Matl
    Use zone_module,          Only: Zone
    use property_data_module, only: Viscoplastic_Model
    Use time_step_module,     Only: dt

    Real(KIND = real_kind),Intent(OUT)   :: Strain_Rate
    integer(KIND = int_kind),Intent(IN)  :: icell             
    Real(KIND = real_kind),Intent(IN)    :: Stress

    Real(KIND = real_kind)               :: Mid_Temp
    Real(KIND = real_kind)               :: rate_mat
    integer(KIND = int_kind)             :: imat, islot

    ! Initialize Strain_Rate
    Strain_Rate = 0.0
    Mid_Temp = 0.5 * (Zone(icell)%Temp + Zone(icell)%Temp_old)
    do imat = 1,nmat
       if (isImmobile(imat)) then
          MAT_SLOT_LOOP: do islot = 1, mat_slot
             if (Matl(islot)%Cell(icell)%Id == imat) then
                select case (Viscoplastic_Model(imat))
                case('elastic_only')
                   exit MAT_SLOT_LOOP
                case('mts')
                   Call MTS_STRAIN_RATE(Stress, Zone(icell)%Temp, &
                        imat, rate_mat, dt)
                   Strain_Rate = Strain_Rate + &
                        Matl(islot)%Cell(icell)%Vof * rate_mat
                case('power_law')
                   Call POWER_LAW_STRAIN_RATE(Stress, Zone(icell)%Temp, &
                        imat, rate_mat, dt)
                   Strain_Rate = Strain_Rate + &
                        Matl(islot)%Cell(icell)%Vof * rate_mat
                end select
             end if
          end do MAT_SLOT_LOOP
       end if
    end do

  end Subroutine VISCOPLASTIC_STRAIN_RATE
  !
  !--------------------------------------------------------------------------
  !
  Subroutine MTS_STRAIN_RATE(stress, temp, mm, strain_rate, dt)
  !
  !--------------------------------------------------------------------------
  !
  ! Purpose:
  !
  !    Calculate the strain rate for a material using the MTS model
  !
  !--------------------------------------------------------------------------
    Use property_data_module, only : MTS_k, MTS_mu_0, MTS_sig_a,     &
                                     MTS_d, MTS_temp_0, MTS_b,       &
                                     MTS_edot_0i, MTS_g_0i, MTS_q_i, &
                                     MTS_p_i, MTS_sig_i

    ! Argument list
    Real(KIND = real_kind),intent(IN)        :: stress, temp
    Integer(KIND = int_kind),intent(IN)      :: mm
    Real(KIND = real_kind),intent(IN)        :: dt
    Real(KIND = real_kind),intent(OUT)       :: strain_rate

    ! Local Variables
    Real(KIND = real_kind)        :: a, c, mu, epsdot_min, eps_sig_a, k
    real (real_kind)              :: tmp

    ! Temperature dependent shear modulus
    mu = MTS_mu_0(mm) - MTS_d(mm)/(exp(MTS_temp_0(mm)/temp) - 1.0)

    if (stress > (MTS_sig_a(mm) + MTS_sig_i(mm)*mu/MTS_mu_0(mm))) then
       strain_rate = MTS_edot_0i(mm)
       return
    end if

    ! We need these terms in any case
    c = MTS_k(mm) * temp / (mu * MTS_b(mm)**3 * MTS_g_0i(mm))
    epsdot_min = 1.0e-20 / dt

    if (stress < MTS_sig_a(mm)) then
       ! Get strain rate for the athermal stress MTS_sig_a
       eps_sig_a = MTS_edot_0i(mm) / exp(1.0 /c)
       if (eps_sig_a < epsdot_min) then
          strain_rate = 0.0
          return
       end if
       ! Assume a sigma^5 relationship for strain rate
       k = log(eps_sig_a) - 5 * log(MTS_sig_a(mm))
       k = exp(k)
       strain_rate = k * stress**5
    else

    ! Rate dependence of the flow (yield) strength
       a = MTS_mu_0(mm)/mu/MTS_sig_i(mm)
       if (a*(stress-MTS_sig_a(mm)) < 0.0_real_kind &
       .or. (a*(stress-MTS_sig_a(mm)) == 0.0_real_kind .and. MTS_p_i(mm) == 0.0_real_kind)) then
          call TLS_panic ('MTS_STRAIN_RATE: invalid exponentiation')
       else
          tmp = 1.0-(a*(stress - MTS_sig_a(mm)))**MTS_p_i(mm)
       end if
       if (tmp < 0.0_real_kind .or. (tmp == 0.0_real_kind .and. MTS_q_i(mm) == 0.0_real_kind)) then
          call TLS_panic ('MTS_STRAIN_RATE: invalid exponentiation')
       else
          strain_rate = MTS_edot_0i(mm) / (exp(tmp**MTS_q_i(mm) /c))
       end if
    end if

    if (strain_rate < epsdot_min) strain_rate = 0.0
    return

  end Subroutine MTS_STRAIN_RATE
  !
  !
  !--------------------------------------------------------------------------
  !
  Subroutine POWER_LAW_STRAIN_RATE(stress, temp, mm, strain_rate, dt)
  !
  !--------------------------------------------------------------------------
  !
  ! Purpose:
  !
  !    Calculate the strain rate for a material using power law model
  !
  !--------------------------------------------------------------------------
    Use property_data_module, only : Pwr_Law_A, Pwr_Law_N, Pwr_Law_Q, Pwr_Law_R

    ! Argument list
    Real(KIND = real_kind),intent(IN)        :: stress, temp
    Integer(KIND = int_kind),intent(IN)      :: mm
    Real(KIND = real_kind),intent(IN)      :: dt
    Real(KIND = real_kind),intent(OUT)       :: strain_rate

    ! Local Variables
    Real(KIND = real_kind)        :: epsdot_min
    
    strain_rate = Pwr_Law_A(mm) * stress**Pwr_Law_N(mm) * exp(- Pwr_Law_Q(mm) / Pwr_Law_R(mm) / temp)

    if (strain_rate > 1.0e7) strain_rate = 1.0e7

    epsdot_min = 1.0e-20 / dt

    if (strain_rate < epsdot_min) strain_rate = 0.0

  end Subroutine POWER_LAW_STRAIN_RATE
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

    Use parameter_module, Only: ndim, nnodes, ncells

    Implicit None

    ! Argument list
    Real(KIND = real_kind), Dimension(ndim*nnodes),   Intent(IN)  :: u
    Real(KIND = real_kind), Dimension(ncomps,ncells), Intent(OUT)  :: Strain
    Real(KIND = real_kind), Dimension(ncells), Intent(OUT), Optional  :: Rotation

    ! Local variables
    Integer(KIND =  int_kind)                              :: inodes
    Real   (KIND = real_kind), Dimension(nnodes)           :: Q
    Real   (KIND = real_kind), Dimension(ncells)           :: dQ_dx, dQ_dy, dQ_dz
    Real   (KIND = real_kind), Dimension(ndim,ndim,ncells) :: Grad_u
    Real   (KIND = real_kind), Allocatable, Dimension(:)   :: R1, R2, R3

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
    Return

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

    Use constants_module,   Only: one
    Use discrete_op_module, Only: DETERMINANT_VOL_AVG
    Use gs_module,          Only: EN_GATHER
    Use mesh_module,        Only: Cell, Vertex, Vrtx_Bdy
    Use parameter_module,   Only: nnodes, ncells, nvc

    Implicit None

    ! Argument list
    Real(KIND = real_kind), Dimension(nnodes), Intent(IN)  :: Q
    Real(KIND = real_kind), Dimension(ncells), Intent(OUT) :: dQ_dx
    Real(KIND = real_kind), Dimension(ncells), Intent(OUT) :: dQ_dy
    Real(KIND = real_kind), Dimension(ncells), Intent(OUT) :: dQ_dz

    ! Local Variables
    Real(KIND = real_kind), Dimension(ncells)     :: Q_f
    Real(KIND = real_kind), Dimension(nvc,ncells) :: Q_v
    Real(KIND = real_kind), Dimension(nvc,ncells) :: X_v
    Real(KIND = real_kind), Dimension(nvc,ncells) :: Y_v
    Real(KIND = real_kind), Dimension(nvc,ncells) :: Z_v
    !
    !---------------------------------------------------------------------------
    !
    ! Gather the vertex-centered scalar Q into the cell-vector quantity Q_v
    Call EN_GATHER (Q_v, Q)

    ! Gather the vertex coordinates
    Call EN_GATHER (X_v, Vertex(:)%Coord(1), BOUNDARY=Vrtx_Bdy(1)%Data)
    Call EN_GATHER (Y_v, Vertex(:)%Coord(2), BOUNDARY=Vrtx_Bdy(2)%Data)
    Call EN_GATHER (Z_v, Vertex(:)%Coord(3), BOUNDARY=Vrtx_Bdy(3)%Data)

    ! Calculate the volume-averaged gradient
    Call DETERMINANT_VOL_AVG (Q_v, Y_v, Z_v, dQ_dx)
    Call DETERMINANT_VOL_AVG (X_v, Q_v, Z_v, dQ_dy)
    Call DETERMINANT_VOL_AVG (X_v, Y_v, Q_v, dQ_dz)

    ! Normalize the volume-averaged gradient by the cell volume
    Q_f(:)   = one/Cell(:)%Volume
    dQ_dx(:) = dQ_dx(:)*Q_f(:)
    dQ_dy(:) = dQ_dy(:)*Q_f(:)
    dQ_dz(:) = dQ_dz(:)*Q_f(:)

    Return

  End Subroutine DISPLACEMENT_GRADIENT
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

    Use constants_module,     Only: zero
    Use parameter_module,     Only: ncells, ndim, ncomps
    use solid_mechanics_data, only: Lame1, Lame2
    use mesh_module,          only: Mesh, GAP_ELEMENT_1
    Implicit None

    ! Argument list
    Real(KIND = real_kind), Dimension(ncomps,ncells), Intent(IN) :: Strain
    Real(KIND = real_kind), Dimension(ncomps,ncells), Intent(out):: Stress
    !
    ! Local variables
    Integer(KIND =  int_kind)                    :: icomps
    Real   (KIND = real_kind), Dimension(ncells) :: Dilatation
    !
    !---------------------------------------------------------------------------
    !
    ! Compute the dilatation
    Dilatation(:) = zero

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

    Return

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

    Use constants_module,     Only: zero
    Use parameter_module,     Only: ndim, ncomps
    use solid_mechanics_data, only: Lame1, Lame2
    use mesh_module,          only: Mesh, GAP_ELEMENT_1
    Implicit None

    ! Argument list
    Real(KIND = real_kind), Dimension(ncomps), Intent(IN) :: Strain
    Real(KIND = real_kind), Dimension(ncomps), Intent(out):: Stress
    Integer(KIND =  int_kind),Intent(IN)                  :: icell
    !
    ! Local variables
    Integer(KIND =  int_kind)                    :: icomps
    Real   (KIND = real_kind)                    :: Dilatation
    !
    !---------------------------------------------------------------------------
    !
    ! Return zero  stresses for gap elements

    if (Mesh(icell)%Cell_Shape >= GAP_ELEMENT_1) then 
       Stress = 0.0
    else
       ! Compute the dilatation
       Dilatation = zero

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

    Return

  End Subroutine MATERIAL_STRESSES_ONE_CELL

  Subroutine DEVIATORIC_STRESS(Stress, Dev_Stress)
    ! Calculate the deviatoric stress tensor from the elastic stress tensor
    Use parameter_module,     Only: ndim, ncomps
    Implicit None

    ! Argument list
    Real(KIND = real_kind), Dimension(ncomps), Intent(IN) :: Stress
    Real(KIND = real_kind), Dimension(ncomps), Intent(out):: Dev_Stress
    !
    ! Local variables
    Integer(KIND =  int_kind)                    :: idim
    Real   (KIND = real_kind)                    :: Mean_Stress

    Mean_Stress = (Stress(1) + Stress(2) + Stress(3))/3.0
    do idim = 1,ndim
       Dev_Stress(idim) = Stress(idim) - Mean_Stress
    end do
    do idim = ndim+1,ncomps
       Dev_Stress(idim) = Stress(idim)
    end do

  end Subroutine DEVIATORIC_STRESS

end Module VISCOPLASTICITY

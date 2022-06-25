!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

MODULE EDIT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define quantities and procedures for short edits
  !   to the "prefix.log" output file.
  !
  !   Public Interface:
  !
  !     * call EDIT_SHORT ()
  !
  !         Perform a "short edit" by computing and printing global
  !         diagnostics.
  !
  ! Contains: EDIT_SHORT
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use parameter_module, only: mops
  use truchas_logging_services
  use scalar_func_class
  implicit none
  private

  public :: EDIT_SHORT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! OUTPUTS namelist edit variables.
  integer, dimension(mops), save, public :: Short_Output_Dt_Multiplier

  ! Edit flags.
  logical, save, public :: short_edit

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  !! Acquire the mesh to pass as an argument to the real edit_short
  subroutine edit_short
    use mesh_manager, only: unstr_mesh_ptr
    use unstr_mesh_type
    type(unstr_mesh), pointer :: mesh
    mesh => unstr_mesh_ptr('MAIN')
    INSIST(associated(mesh))
    call edit_short_aux(mesh)
  end subroutine

  SUBROUTINE EDIT_SHORT_AUX (mesh)
    !=======================================================================
    ! Purpose(s):
    !
    !   Perform a "short edit" by computing and printing global diagnostics.
    !
    !=======================================================================
    use cutoffs_module,         only: alittle
    use matl_module,            only: GATHER_VOF
    use truchas_env,            only: output_file_name
    use parameter_module,       only: nmat
    use parallel_communication
    use material_model_driver,  only: matl_model
    use time_step_module,       only: cycle_number, t
    use zone_module,            only: Zone
    use output_utilities,       only: ANNOUNCE_FILE_WRITE
    use unstr_mesh_type

    type(unstr_mesh), intent(in) :: mesh

    ! Local Variables
    integer, parameter :: ndim = 3
    character(LEN = 128) :: string, string2
    integer :: i, m, n, variables = 2*ndim + 6, nmechvar = 4
    integer, dimension(1) :: MaxLoc_L, MinLoc_L
    real(r8), dimension(mesh%ncell_onP) :: Enthalpy, KE, Mass, Matl_Mass, Tmp, Matl_Vol
    real(r8) :: Temperature
    !type(CELL_MECH_INVARIANT), pointer, dimension(:) :: mech_info => NULL()

    real(r8), dimension(nmat) :: Material_Enthalpy, Material_KE, Material_Volume, Material_Mass
    real(r8) :: Material_Momentum(ndim,nmat)
    real(r8) :: Total_Enthalpy, total_KE, total_mass, Total_Momentum(ndim), total_volume
    class(scalar_func), allocatable :: f

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Write out the header, time, and cycle
    write (string, 10) t, cycle_number
10  format (1x,14('* '),'SHORT EDIT: t = ',1pe13.5,', step = ',i6,14(' *'))
    call TLS_info ('')
    call TLS_info (string)
    call TLS_info ('')

    ! Compute mass, momentum, kinetic energy, enthalpy, and
    ! total energy for each material. Also compute the sum
    ! over all materials. Print out the results to out_tbrook.
    call TLS_info ('')
    call TLS_info (repeat(' ',32) // 'Mass-Momentum-Energy Summary')
    call TLS_info ('')
    call TLS_info (' Matl    Name     Volume      Mass               (X,Y,Z) Momentum               KE        Enthalpy')
    call TLS_info (' ----    ----     ------      ----               ----------------               --        --------')

    ! Preinitialize
    Mass = 0.0_r8; Enthalpy = 0.0_r8; KE = 0.0_r8
    total_mass = 0.0_r8; total_momentum = 0.0_r8
    total_KE = 0.0_r8; total_enthalpy = 0.0_r8
 
    ! Zone kinetic energy density
    do n = 1,ndim
       KE = KE + Zone%Vc(n)**2
    end do
    KE = 0.5_r8*KE

    ! Loop over each material
    MATERIAL_SUMS: do m = 1, matl_model%nphase_real  ! nmat or nmat-1 if have_void

       ! Gather vof for this material
       call GATHER_VOF (m, Tmp)

       ! Sum material volume.
       Matl_Vol = mesh%volume(:mesh%ncell_onP) * Tmp
       Material_Volume(m) = global_sum(Matl_Vol)
       if (ABS(Material_Volume(m)) <= alittle) Material_Volume(m) = 0.0_r8

       ! Sum material mass.
       Matl_Mass = Matl_Vol*matl_model%const_phase_prop(m, 'density')
       Material_Mass(m) = global_sum(Matl_Mass)
       if (ABS(Material_Mass(m)) <= alittle) Material_Mass(m) = 0.0_r8

       ! Accumulate the total mass.
       Mass = Mass + Matl_Mass
       total_mass = total_mass + Material_Mass(m)

       ! Compute material momentum.
       MOMENTUM: do n = 1,ndim
          Material_Momentum(n,m) = global_sum(Matl_Mass*Zone%Vc(n))
          if (ABS(Material_Momentum(n,m)) <= alittle .or. .not.matl_model%is_fluid(m)) Material_Momentum(n,m) = 0.0_r8
          Total_Momentum(n) = Total_Momentum(n) + Material_Momentum(n,m)
       end do MOMENTUM

       ! Compute material kinetic energy.
       Material_KE(m) = global_sum(Matl_Mass*KE)
       if (ABS(Material_KE(m)) <= alittle .or. .not.matl_model%is_fluid(m)) Material_KE(m) = 0.0_r8
       total_KE = total_KE + Material_KE(m)

       ! Get the material enthalpy.
       call matl_model%alloc_phase_prop(m, 'enthalpy', f)
       ENTHALPY_LOOP: do n=1,mesh%ncell_onP
         Tmp(n) = Matl_Vol(n) * f%eval([Zone(n)%Temp])
       end do ENTHALPY_LOOP

       ! Accumulate the material enthalpy.
       Material_Enthalpy(m) = global_sum(Tmp)
       if (ABS(Material_Enthalpy(m)) <= alittle) Material_Enthalpy(m) = 0.0_r8

       ! Accumulate the total enthalpy.
       Enthalpy = Enthalpy + Tmp
       total_enthalpy = total_enthalpy + Material_Enthalpy(m)

       ! Write out the summaries for this material.
       write (string, 25) TRIM(matl_model%phase_name(m)), Material_Volume(m),      &
                                 Material_Mass(m), (Material_Momentum(n,m),n=1,ndim), &
                                 Material_KE(m), Material_Enthalpy(m)
25     format (6x,a8,1x,1pe11.4,1x,1pe11.4,1x,'(',1pe11.4,',',1pe11.4,',', &
               1pe11.4,')',1pe11.4,1x,1pe11.4)
       call TLS_info (string)

    end do MATERIAL_SUMS

    total_volume = global_sum(mesh%volume(:mesh%ncell_onP))

    ! Write out the totals.
    write (string,30) total_volume, total_mass, (Total_Momentum(n),n=1,ndim), total_KE, total_enthalpy
30  format (1x,'Totals',8x,1pe11.4,1x,1pe11.4,1x,'(',1pe11.4,',',1pe11.4,',',1pe11.4,')',1pe11.4,1x,1pe11.4)
    call TLS_info ('                ----------  ---------- ------------------------------------- ----------  ----------')
    call TLS_info (string)       
    
    ! Compute and write out global extrema.
    call TLS_info ('')
    call TLS_info ('                                    Global Extrema Summary')
    call TLS_info ('')
    call TLS_info ('                        Quantity      Minimum      Cell    Maximum      Cell')
    call TLS_info ('                        --------      -------      ----    -------      ----')

    ! Loop over all variables whose global extrema will be printed out.
    VARIABLE_EXTREMA: do i = 1,variables

       ! Select the appropriate variable
       select case (i)
          case (1:ndim)
             ! Velocity
             Tmp = Zone%Vc(i)
             if (i == 1) then
                string = 'X-Velocity'
             else if (i == 2) then
                string = 'Y-Velocity'
             else if (i == 3) then
                string = 'Z-Velocity'
             end if
          case (ndim+1:2*ndim)
             ! Momentum
             n = i - ndim
             Tmp = Mass*Zone%Vc(n)
             if (n == 1) then
                string = 'X-Momentum'
             else if (n == 2) then
                string = 'Y-Momentum'
             else if (n == 3) then
                string = 'Z-Momentum'
             end if
          case (2*ndim+1)
             ! Kinetic Energy
             Tmp = Mass*KE
             string = 'KE'
          case (2*ndim+2)
             ! Mass
             Tmp = Mass
             string = 'Mass'
          case (2*ndim+3)
             ! Density
             Tmp = Zone%Rho
             string = 'Density'
          case (2*ndim+4)
             ! Pressure
             Tmp = Zone%P
             string = 'Pressure'
          case (2*ndim+5)
             ! Temperature
             Tmp = Zone%Temp
             string = 'Temperature'
          case (2*ndim+6)
             ! Enthalpy
             Tmp = Enthalpy
             string = 'Enthalpy'
       end select

       ! Find the extrema location
       MinLoc_L = global_minloc(Tmp)
       MaxLoc_L = global_maxloc(Tmp)

       ! Write out the extrema location and values
       write (string2, 40) string, &
                                 global_minval(Tmp), MinLoc_L, &
                                 global_maxval(Tmp), MaxLoc_L
40     format (23x,a11,2x,1pe11.4,1x,i8,2x,1pe11.4,1x,i8)
       call TLS_info (string2)

    end do VARIABLE_EXTREMA

    ! Write out the trailing header.
    write (string, 10) t, cycle_number
    call TLS_info ('')
    call TLS_info (string)
    call TLS_info ('')

    call ANNOUNCE_FILE_WRITE ('Short edit', output_file_name('log'))

  END SUBROUTINE EDIT_SHORT_AUX

END MODULE EDIT_MODULE

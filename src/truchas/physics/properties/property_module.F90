!!
!! PROPERTY_MODULE
!!
!! NNC, Dec 1012. The module have been extensively restructured. Most obvious
!! is the removal of the PROPERTY subroutine and switch to using property
!! definitions from the phase property table (PHASE namelist).  A number of
!! methods for evaluating specific cell-averaged properties are provided.
!! Still here are several flow specific initialization routines.  Further
!! refactoring is planned.
!!

#include "f90_assert.fpp"

module property_module

  use kinds, only: r8
  use parameter_module, only: maxmat
  use truchas_logging_services
  implicit none
  private
 
  public :: set_user_material_id, get_user_material_id, get_truchas_material_id
  public :: get_viscosity, FLUID_PROPERTIES
  public :: get_density, get_density_delta
  public :: EM_permittivity, EM_permeability, EM_conductivity
  public :: request_fluid_property
  public :: density_material, enthalpy_density_material

  ! Private source of data for user material number
  integer, dimension(0:maxmat) :: User_Material_Number = -1
 
contains

  subroutine set_user_material_id (id_truchas, id_user)
    use parameter_module, only: maxmat
    integer, intent(in) :: id_truchas, id_user
    INSIST(id_truchas >= 1 .and. id_truchas <= maxmat)
    user_material_number(id_truchas) = id_user
  end subroutine set_user_material_id

  integer function get_user_material_id (id_truchas) result (id_user)
    use parameter_module, only: nmat
    integer, intent(in) :: id_truchas
    if (id_truchas >= 1 .and. id_truchas <= nmat) then
      id_user = user_material_number(id_truchas)
    else if (id_truchas == 0) then
      id_user = 0 
    else
      INSIST(.false.)
    end if
  end function get_user_material_id

  integer function get_truchas_material_id (id_user) result (id_truchas)
    use parameter_module, only: nmat
    use constants_module, only: ipreset
    integer, intent(in) :: id_user
    if (id_user == 0) then
      id_truchas = 0
    else 
      do id_truchas = 1, nmat
        if (user_material_number(id_truchas) == id_user) return
      end do
      id_truchas = ipreset  ! This is actually used!
    end if
  end function get_truchas_material_id

  subroutine get_density (temp, value)
    use kinds, only: r8
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: value(:)
    call compute_cell_property ('density', temp, value)
  end subroutine get_density


  subroutine get_density_delta (temp, drho)

    use kinds, only: r8
    use parameter_module, only: ncells, nmat
    use phase_property_table
    use material_interop, only: void_material_index, material_to_phase
    use matl_module, only: gather_vof
    use fluid_data_module, only: isImmobile
    use scalar_functions

    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: drho(:)

    integer :: m, j, phase_id, prop_id
    real(r8) :: vofm(ncells), state(1), mval, rhom
    type(scafun), pointer :: prop_fun

    ASSERT(size(temp) == ncells)
    ASSERT(size(drho) == ncells)

    ASSERT(ppt_has_property('density deviation'))
    prop_id = ppt_property_id('density deviation')

    drho = 0.0_r8
    do m = 1, nmat
      if (m == void_material_index .or. isImmobile(m)) cycle
      phase_id = material_to_phase(m)
      ASSERT(phase_id > 0)
      call ppt_get_phase_property (phase_id, prop_id, prop_fun)
      ASSERT(associated(prop_fun))
      call gather_vof (m, vofm)
      rhom = density_material(m)
      if (is_const_scafun(prop_fun)) then
        mval = eval(prop_fun, state)  ! state is ignored, but required
        drho = drho + (rhom*mval)*vofm
      else
        do j = 1, ncells
          if (vofm(j) > 0.0_r8) then
            state(1) = temp(j)
            mval = eval(prop_fun, state)
            drho(j) = drho(j) + rhom*mval*vofm(j)
          end if
        end do
      end if
    end do

  end subroutine get_density_delta


  subroutine get_viscosity (mu)
  
    use parameter_module, only: ncells, nmat
    use fluid_data_module, only: IsImmobile, FluidRho
    use matl_module, only: gather_vof
    use turbulence_module, only: turbulence_model, Nu_Turb
    use zone_module, only: Zone

    real(r8), intent(out) :: mu(:)
    
    integer :: m
    real(r8) :: vofm(ncells), fluid_vof(ncells)
    
    call compute_cell_property ('viscosity', Zone%Temp, mu, fluid=.true.)
    
    !! Get the fluid volume fraction, including void.
    fluid_vof = 0.0_r8
    do m = 1, nmat
      if (IsImmobile(m)) cycle
      call gather_vof (m, vofm)
      fluid_vof = fluid_vof + vofm
    end do
    
    where (fluid_vof > 0.0_r8) mu = mu / fluid_vof
    
    if (turbulence_model == 'alg') mu = mu + FluidRho*Nu_Turb
 
  end subroutine get_viscosity


  function EM_permittivity () result (value)
    use kinds, only: r8
    use parameter_module, only: ncells
    use zone_module, only: zone
    real(r8) :: value(ncells)
    call compute_cell_property ('electric susceptibility', zone%temp, value)
    value = 1.0_r8 + value
  end function EM_permittivity

  function EM_permeability () result (value)
    use kinds, only: r8
    use parameter_module, only: ncells
    use zone_module, only: zone
    real(r8) :: value(ncells)
    call compute_cell_property ('magnetic susceptibility', zone%temp, value)
    value = 1.0_r8 + value
  end function EM_permeability

  function EM_conductivity () result (value)
    use kinds, only: r8
    use parameter_module, only: ncells
    use zone_module, only: zone
    real(r8) :: value(ncells)
    call compute_cell_property ('electrical conductivity', zone%temp, value)
  end function EM_conductivity


  function enthalpy_density_material (matl_id, temp) result (value)
    use material_interop, only: ds_enthalpy_density
    integer, intent(in) :: matl_id
    real(r8), intent(in) :: temp
    real(r8) :: value
    if (matl_id == 0) then ! somebody calls this with a zero ID, argh!
      value = 0.0_r8
    else
      value = ds_enthalpy_density(matl_id, [temp])
    end if
  end function enthalpy_density_material

  function density_material (matl_id) result (value)
    use property_data_module, only: density
    use material_interop, only: ds_density
    integer, intent(in) :: matl_id
    real(r8) :: value, state(1)
    if (matl_id == 0) then ! somebody calls this with a zero ID, argh!
      value = density(matl_id)
    else
      value = ds_density(matl_id, state)
    end if
  end function density_material


  SUBROUTINE FLUID_PROPERTIES (SkipFlow)
    !=======================================================================
    ! Purpose(s):
    !
    !   Evaluate cell properties EXCLUDING immobile materials.  Check that
    !   there are some cells in which there are flow equations to solve
    !
    !======================================================================
    use bc_module,              only: BC_Mat, BC_Vel
    use constants_module,       only: zero, ipreset
    use cutoffs_module,         only: cutvof
    use fluid_data_module,      only: fluid_flow,                             &
                                      fluidRho, fluidVof, fluidDeltaRho,      &
                                      RealFluidVof, cutRho, fluid_cutoff,     &
                                      Solid_Face, isPureImmobile, isImmobile, &
                                      boussinesq_approximation,               &
                                      Cell_isnt_Void, Ngbr_isnt_Void,         &
                                      Fluxing_Velocity,MinFluidRho
    use gs_module,              only: EE_GATHER
    use matl_module,            only: GATHER_VOF
    use mesh_module,            only: Cell
    use parameter_module,       only: ncells, nfc, ndim, nmat
    use projection_data_module, only: Boundary_Flag
    use pgslib_module,          only: PGSLIB_GLOBAL_ALL, PGSLIB_GLOBAL_ANY, PGSLIB_GLOBAL_MINVAL
    use zone_module,            only: Zone
    use material_interop,       only: void_material_index

    ! Argument List
    logical :: SkipFlow  ! Flag indicating that the flow solution should abort

    ! Local Variables
    integer :: status, m, n, f, nc, c 
    real(r8) :: BC_V_Dot_N, rhom
    real(r8), allocatable :: Prop(:), DeltaProp(:), DeltaFluidRho(:), mVof(:)
    logical,  allocatable :: FlowInCell(:), isPureImmobileNgbr(:,:), Mask(:)

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (Prop(ncells),                   &
              DeltaProp(ncells),              &
              FlowInCell(ncells),             &
              isPureImmobileNgbr(nfc,ncells), &
              Mask(ncells),                   &
              DeltaFluidRho(ncells),          &
              mVof(ncells),                   STAT = status)
    call TLS_fatal_if_any (status/=0, 'FLUID_PROPERTIES: failed allocation')

    RealFluidVof  = zero  ! Vof of fluids that are not void
    cutRho        = zero  ! used for mass cutoff
    DeltaFluidRho = zero  ! Need to be sure to zero the density change from reference
    
    fluidRho = zero      ! density of (remaining) fluid in cell
    fluidVof = zero      ! volume fraction of fluid in cell
    Zone%Rho = zero      ! cell density
    !fluidDeltaRho = zero ! density variation (due to temperature)

    do m = 1, nmat
       call gather_vof(m, mVof)
       !mMask(:) = (mVof > zero)
     
       !call PROPERTY (Zone%Temp, m, 'density', &
       !               Value = Prop, material_mask = mMask)
       rhom = density_material(m)
      
       ! calculate the contribution of material m to each cell density (Zone%Rho)
       !Zone%Rho = Zone%Rho + mVof * Prop
       Zone%Rho = Zone%Rho + mVof * rhom

       ! if material m flows, calculate fluid density & fraction (fluidRho, fluidVof);
       ! also calculate the deviation from the reference density, and add it into fluidRho
       ! the Boussinesq Option will remove the deviation immediately after the body force calculation
       if (isImmobile(m)) CYCLE
       !call PROPERTY (Zone%Temp, m, 'density', &
       !               PROPERTY_TYPE_EVALUATION = 'delta', &
       !               Value = DeltaProp, material_mask = mMask)

       !fluidRho      = fluidRho      + mVof * Prop
       fluidRho      = fluidRho      + mVof * rhom
       !fluidDeltaRho = fluidDeltaRho + mVof * DeltaProp
       fluidVof      = fluidVof      + mVof
       
       ! realFluidVof is kept as a fluid-flow variable since it's useful there
       ! cutRho is used as a cutoff metric for the velocity update in the predictor
       !where (Prop > zero) 
       !   RealFluidVof = RealFluidVof + mVof
       !   cutRho       = cutRho + cutvof*Prop
       !end where
       if (m /= void_material_index) then
          RealFluidVof = RealFluidVof + mVof
          cutRho = cutRho + merge(cutvof*rhom, 0.0_r8, mVof > 0)
       end if

    end do
    
    !call get_density (Zone%Temp, Zone%Rho)
    !call get_density (Zone%Temp, fluidRho, fluid=.true.)
    
    call get_density_delta (Zone%Temp, fluidDeltaRho)

    ! Set void cell indicator arrays.
    Cell_isnt_Void = Zone%Rho > zero
    call EE_GATHER(Ngbr_isnt_Void, Cell_isnt_Void)

    if (boussinesq_approximation) then
       where (fluidVof > zero)
          fluidRho      = fluidRho / fluidVof
          fluidDeltaRho = fluidDeltaRho / fluidVof
       end where
    else
       where (fluidVof > zero)
          fluidRho       = (fluidRho+fluidDeltaRho) / fluidVof
          fluidDeltaRho  = zero
       end where
    endif

    ! Evaluate and store the minimum 'real' fluid density for evaluation of the
    ! face density limit
    Mask = RealFluidVof > zero
    Prop = zero
    where(Mask) Prop = fluidRho*fluidVof/RealFluidVof
    MinFluidRho = PGSLIB_GLOBAL_MINVAL(prop, MASK=Mask)

    ! Where cells are void (rho = zero), zero out velocity & pressure
    Mask = .false.
    where (FluidRho == zero ) Mask = .true.
    do n = 1,ndim
       where (Mask) Zone%Vc(n) = zero
    end do
!!!!   remove this -- we must let the projection solve impose this condition
!!!!      where (Mask) Zone%P = void_pressure

    ! Since the last flow call, the phase change routines may have solidified material
    ! that was fluid during the previous timestep; so we need to recalculate arrays
    ! that tell us where the fluid material is.
    isPureImmobile = .false.
    call TURN_OFF_FLOW(IsPureImmobile)
    where (fluidVof < fluid_cutoff) isPureImmobile = .true.
    call EE_Gather (isPureImmobileNgbr,isPureImmobile)
    Solid_Face = .false.
    do f = 1, nfc
       where (isPureImmobile) Solid_Face(f,:) = .true.
       where (isPureImmobileNgbr(f,:)) Solid_Face(f,:) = .true.
    end do

    ! Zero out velocities for any cells that have completely solidified
    ! between the last flow call and this one.
    do nc = 1, ncells
       if (isPureImmobile(nc)) then
          Zone(nc)%Vc(:)  = zero
          Zone(nc)%Vc_old(:) = zero
       endif
       do f= 1, nfc
          if(Solid_Face(f,nc)) then
             Fluxing_Velocity(f,nc)   = zero
          endif
       end do
    end do
    ! Check that there is some fluid to flow, that there are cells with
    ! fluidVof > fluid_cutoff that also have one or more Solid_Face's == .false.

    SkipFlow = .true.       ! true if there is nowhere to solve for flow

    if(fluid_flow) then
       ! Check bcs for possible inflow of real fluid
       SKIPFLOWLOOP: do c = 1,ncells
          if(isPureImmobile(c)) cycle
          do f = 1, nfc
              if(Boundary_Flag(f,c) == 1) then
              ! Dirichlet Pressure BC
                  if(BC_Mat(f,c) /= ipreset) then
                     if(DENSITY_MATERIAL(BC_Mat(f,c)) > zero ) then
                        SkipFlow = .false.
                        EXIT SKIPFLOWLOOP
                    endif
                  endif
              endif
              if(Boundary_Flag(f,c) == 2) then
              ! Dirichlet Velocity BC
                  if(BC_Mat(f,c) /= ipreset) then
                     if(DENSITY_MATERIAL(BC_Mat(f,c)) > zero ) then
                        !Inward Flow?
                        BC_V_Dot_N = 0.0
                        do n = 1, ndim
                           BC_V_Dot_N = BC_V_Dot_N + BC_Vel(n,f,c)*Cell(c)%Face_Normal(n,f)
                        enddo
                        if(BC_V_Dot_N < zero) then
                           SkipFlow = .false.
                           EXIT SKIPFLOWLOOP
                        endif
                     endif
                  endif
              endif
          enddo 
       enddo SKIPFLOWLOOP
       SkipFlow = PGSLIB_GLOBAL_ALL(SkipFlow)

       ! There are no inflow bcs with real fluid, check for flow within the mesh
       if(SkipFlow) then
           FlowInCell = .false.
           do f = 1,nfc
              where (.not. Solid_Face(f,:) .and. Boundary_Flag(f,:) /= 0) FlowInCell = .true.
           end do
           where (RealFluidVof < fluid_cutoff) FlowInCell = .false.
           if (PGSLIB_GLOBAL_ANY(FlowInCell)) SkipFlow = .false.
       endif
    endif

    DEALLOCATE (Prop)
    DEALLOCATE (DeltaProp) 
    DEALLOCATE (FlowInCell)
    DEALLOCATE (isPureImmobileNgbr)
    DEALLOCATE (Mask)
    DEALLOCATE (DeltaFluidRho)
    DEALLOCATE (mVof)

  END SUBROUTINE FLUID_PROPERTIES

  SUBROUTINE TURN_OFF_FLOW(IsPureImmobile)
    !=======================================================================
    ! Purpose:
    !
    !  reset saved flow solution quantities at the start of restart calculations
    !
    !    Jim Sicilian, LANL CCS-2 (sicilian@lanl.gov)   Dec 2002
    !
    !=======================================================================

    use mesh_module,           only: Cell
    use parameter_module,      only: ncells, ndim
    use region_data,           only: Regions, nregion

    ! calling arguments
    logical, intent(inout) :: IsPureImmobile(:)

    ! Local Variables
    integer :: n, status, i, d
    logical :: inside

    logical, allocatable :: flow_off_cell(:)

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ALLOCATE (flow_off_cell(ncells), STAT = status)
    call TLS_fatal_if_any (status/=0, 'TURN_OFF_FLOW: failed allocation')

    ! initialize flow_off_cell to false
    flow_off_cell(:) = .false.
    REGION_LOOP: do n=1,nregion
      if(Regions(n)%flow_off) then
         ! do the geometry check with the cell center
         do i=1,ncells
            inside = .true.
            do d=1,ndim
               if (Cell(i)%Centroid(d) < Regions(n)%x1(d) .or. Cell(i)%Centroid(d) > Regions(n)%x2(d)) then
                  inside = .false.
               end if
            end do
            if (inside) then
               flow_off_cell(i) = .true.
            end if
         end do
      end if
    end do REGION_LOOP

    do i=1,ncells
      if (flow_off_cell(i)) then
        IsPureImmobile(i) = .true.
      end if
    end do
  
    DEALLOCATE(flow_off_cell)

  END SUBROUTINE TURN_OFF_FLOW

  !!
  !! COMPUTE_CELL_PROPERTY
  !!
  !! This computes the specified property for the cells on the original Truchas
  !! mesh. The property is one given in a PHASE namelist PROPERTY_NAME variable,
  !! or one created internally by Truchas, and the associated property value is
  !! either constant or a function of temperature only.  The routine essentially
  !! handles the material averaging of the property over a cell using the volume
  !! fraction data from MATL.  The void material (one with a MATERIAL namelist
  !! density of zero) contributes nothing to the property value; e.g., zero is
  !! returned for an entirely void cell.
  !! 

  subroutine compute_cell_property (prop, temp, value, fluid)

    use kinds, only: r8
    use parameter_module, only: ncells, nmat
    use fluid_data_module, only: isImmobile
    use phase_property_table
    use material_interop, only: void_material_index, material_to_phase
    use scalar_functions
    use matl_module, only: gather_vof

    character(*), intent(in) :: prop
    real(r8), intent(in) :: temp(:)
    real(r8), intent(out) :: value(:)
    logical, intent(in), optional :: fluid

    integer :: m, j, phase_id, prop_id
    real(r8) :: vofm (ncells), state(1), mval
    type(scafun), pointer :: prop_fun
    logical :: fluids_only

    ASSERT(size(temp) == ncells)
    ASSERT(size(value) == ncells)

    ASSERT(ppt_has_property(prop))
    prop_id = ppt_property_id(prop)

    fluids_only = .false.
    if (present(fluid)) fluids_only = fluid

    value = 0.0_r8
    do m = 1, nmat
      if (m == void_material_index) cycle
      if (fluids_only .and. isImmobile(m)) cycle
      phase_id = material_to_phase(m)
      ASSERT(phase_id > 0)
      call ppt_get_phase_property (phase_id, prop_id, prop_fun)
      ASSERT(associated(prop_fun))
      call gather_vof (m, vofm)
      if (is_const_scafun(prop_fun)) then
        mval = eval(prop_fun, state)  ! state is ignored, but required
        value = value + mval*vofm
      else
        do j = 1, ncells
          if (vofm(j) > 0.0_r8) then
            state(1) = temp(j)
            mval = eval(prop_fun, state)
            value(j) = value(j) + mval*vofm(j)
          end if
        end do
      end if
    end do

  end subroutine compute_cell_property

  !!
  !! REQUEST_FLUID_PROPERTY
  !!
  !! This routine ensures that the specified property is defined for all fluid
  !! phases.  If the property is not defined for such a phase, execution halts
  !! with an error message, unless the optional argument DEFAULT is specified.
  !! In this case the property is defined for the phase using this default value.
  !!

  subroutine request_fluid_property (prop, default)

    use phase_property_table
    use parameter_module, only: nmat
    use fluid_data_module, only: isImmobile
    use material_interop, only: void_material_index, material_to_phase
    use scalar_functions

    character(*), intent(in) :: prop
    real(r8), intent(in), optional :: default

    integer :: prop_id, phase_id, m
    type(scafun) :: f_default
    character(128) :: message

    call ppt_add_property (prop, prop_id)
    if (present(default)) call create_scafun_const (f_default, default)

    do m = 1, nmat
      if (m == void_material_index .or. isImmobile(m)) cycle
      phase_id = material_to_phase(m)
      ASSERT(phase_id > 0)
      if (ppt_has_phase_property(phase_id, prop_id)) cycle
      if (present(default)) then
        call ppt_assign_phase_property (phase_id, prop_id, f_default)
        write(message, '(2x,3a,es10.3,3a)') 'Using default value "', trim(prop), &
            '" =', default, ' for phase "', trim(ppt_phase_name(phase_id)), '"'
        call TLS_info (message)
      else
        write(message,'(5a)') 'missing property "', trim(prop), &
            '" for phase "', trim(ppt_phase_name(phase_id)), '"'
        call TLS_fatal (message)
      end if
    end do

    if (present(default)) call destroy (f_default)

  end subroutine request_fluid_property

end module property_module

!!
!! FLOW_PROPERTY_MODULE
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! NNC, Dec 2019. This is the remnant of the original property_module.F90. It
!! adapts the new material model used elsewhere in Truchas to the original API
!! used by the legacy flow model. It is not, and should not be, used elsewhere.
!!

#include "f90_assert.fpp"

module flow_property_module

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use material_model_driver, only: matl_model
  use scalar_func_containers, only: scalar_func_box
  use truchas_logging_services
  implicit none
  private
 
  public :: get_material_id
  public :: get_viscosity, FLUID_PROPERTIES
  public :: get_density, get_density_delta
  public :: density_material, enthalpy_material
  public :: material_name
  public :: have_void
  public :: void_material_index
  public :: isImmobile

  real(r8), allocatable, public :: matl_density(:)
  type(scalar_func_box), allocatable :: matl_enthalpy(:)

contains

  function material_name(n) result(name)
    integer, intent(in) :: n
    character(:), allocatable :: name
    name = matl_model%phase_name(n)
  end function material_name

  integer function get_material_id(name) result(id)
    character(*), intent(in) :: name
    id = matl_model%phase_index(name)
  end function get_material_id

  logical function have_void()
    have_void = matl_model%have_void
  end function have_void

  integer function void_material_index()
    void_material_index = matl_model%void_index
  end function

  logical function isImmobile(n)
    integer, intent(in) :: n
    isImmobile = .not.matl_model%is_fluid(n)
  end function isImmobile

  subroutine get_density (temp, value)
    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: value(:)
    call compute_cell_property('density', temp, value)
  end subroutine get_density


  subroutine get_density_delta (temp, drho)

    use legacy_mesh_api, only: ncells
    use matl_module, only: gather_vof
    use scalar_func_class

    real(r8), intent(in)  :: temp(:)
    real(r8), intent(out) :: drho(:)

    integer :: m, j
    real(r8) :: vofm(ncells)
    class(scalar_func), allocatable :: prop_fun

    ASSERT(size(temp) == ncells)
    ASSERT(size(drho) == ncells)

    drho = 0.0_r8
    do m = 1, matl_model%nphase_real
      if (.not.matl_model%is_fluid(m)) cycle
      call matl_model%alloc_phase_prop(m, 'density-delta', prop_fun)
      ASSERT(allocated(prop_fun))
      call gather_vof (m, vofm)
      do j = 1, ncells
        if (vofm(j) > 0.0_r8) drho(j) = drho(j) + vofm(j)*prop_fun%eval([temp(j)])
      end do
    end do

  end subroutine get_density_delta


  subroutine get_viscosity (mu)
  
    use legacy_mesh_api, only: ncells
    use fluid_data_module, only: FluidRho
    use matl_module, only: gather_vof
    use turbulence_module, only: turbulence_model, Nu_Turb
    use zone_module, only: Zone

    real(r8), intent(out) :: mu(:)
    
    integer :: m
    real(r8) :: vofm(ncells), fluid_vof(ncells)
    
    call compute_cell_property('viscosity', Zone%Temp, mu, fluid=.true.)
    
    !! Get the fluid volume fraction, including void.
    fluid_vof = 0.0_r8
    do m = 1, matl_model%nphase
      if (matl_model%is_fluid(m)) then
        call gather_vof(m, vofm)
        fluid_vof = fluid_vof + vofm
      end if
    end do
    
    where (fluid_vof > 0.0_r8) mu = mu / fluid_vof
    
    if (turbulence_model == 'alg') mu = mu + FluidRho*Nu_Turb
 
  end subroutine get_viscosity


  function enthalpy_material (matl_id, temp) result (value)
    integer, intent(in) :: matl_id
    real(r8), intent(in) :: temp
    real(r8) :: value
    integer :: i
    if (.not.allocated(matl_enthalpy)) then
      allocate(matl_enthalpy(matl_model%nphase_real))
      do i = 1, matl_model%nphase_real
        call matl_model%alloc_phase_prop(i, 'enthalpy', matl_enthalpy(i)%f)
      end do
    end if
    if (matl_id == 0) then ! somebody calls this with a zero ID, argh!
      value = 0.0_r8
    else if (matl_id == matl_model%void_index) then
      value = 0.0_r8
    else
      value = matl_enthalpy(matl_id)%f%eval([temp])
    end if
  end function enthalpy_material

  function density_material (matl_id) result (value)
    integer, intent(in) :: matl_id
    real(r8) :: value
    integer :: j
    if (.not.allocated(matl_density)) then
      allocate(matl_density(0:matl_model%nphase))
      do j = 1, matl_model%nphase_real
        matl_density(j) = matl_model%const_phase_prop(j, 'density')
      end do
      if (matl_model%have_void) matl_density(matl_model%nphase) = 0.0_r8
      matl_density(0) = 0.0_r8 ! somebody calls this with a zero ID, argh!
    end if
    value = matl_density(matl_id)
  end function density_material


  SUBROUTINE FLUID_PROPERTIES (SkipFlow, t)
    !=======================================================================
    ! Purpose(s):
    !
    !   Evaluate cell properties EXCLUDING immobile materials.  Check that
    !   there are some cells in which there are flow equations to solve
    !
    !======================================================================
    use bc_module,              only: BC_Mat, bndry_vel !BC_Vel
    use input_utilities,        only: NULL_I
    use cutoffs_module,         only: cutvof
    use fluid_data_module,      only: fluid_flow,                             &
                                      fluidRho, fluidVof, fluidDeltaRho,      &
                                      RealFluidVof, cutRho, fluid_cutoff,     &
                                      Solid_Face, isPureImmobile,             &
                                      boussinesq_approximation,               &
                                      Cell_isnt_Void, Ngbr_isnt_Void,         &
                                      Fluxing_Velocity,MinFluidRho
    use matl_module,            only: GATHER_VOF
    use legacy_mesh_api,        only: ncells, nfc, ndim, Cell, EE_GATHER
    use projection_data_module, only: Boundary_Flag
    use pgslib_module,          only: PGSLIB_GLOBAL_ALL, PGSLIB_GLOBAL_ANY, PGSLIB_GLOBAL_MINVAL
    use zone_module,            only: Zone

    ! Argument List
    logical :: SkipFlow  ! Flag indicating that the flow solution should abort
    real(r8), intent(in) :: t ! time

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

    RealFluidVof  = 0.0_r8  ! Vof of fluids that are not void
    cutRho        = 0.0_r8  ! used for mass cutoff
    DeltaFluidRho = 0.0_r8  ! Need to be sure to zero the density change from reference
    
    fluidRho = 0.0_r8      ! density of (remaining) fluid in cell
    fluidVof = 0.0_r8      ! volume fraction of fluid in cell
    Zone%Rho = 0.0_r8      ! cell density
    !fluidDeltaRho = 0.0_r8 ! density variation (due to temperature)

    do m = 1, matl_model%nphase
       call gather_vof(m, mVof)
       rhom = density_material(m)
      
       ! calculate the contribution of material m to each cell density (Zone%Rho)
       Zone%Rho = Zone%Rho + mVof * rhom

       ! if material m flows, calculate fluid density & fraction (fluidRho, fluidVof);
       ! also calculate the deviation from the reference density, and add it into fluidRho
       ! the Boussinesq Option will remove the deviation immediately after the body force calculation
       if (.not.matl_model%is_fluid(m)) CYCLE

       fluidRho      = fluidRho      + mVof * rhom
       fluidVof      = fluidVof      + mVof
       
       ! realFluidVof is kept as a fluid-flow variable since it's useful there
       ! cutRho is used as a cutoff metric for the velocity update in the predictor
       if (m /= matl_model%void_index) then
          RealFluidVof = RealFluidVof + mVof
          cutRho = cutRho + merge(cutvof*rhom, 0.0_r8, mVof > 0)
       end if

    end do
    
    !call get_density (Zone%Temp, Zone%Rho)
    !call get_density (Zone%Temp, fluidRho, fluid=.true.)
    
    call get_density_delta (Zone%Temp, fluidDeltaRho)

    ! Set void cell indicator arrays.
    Cell_isnt_Void = Zone%Rho > 0.0_r8
    call EE_GATHER(Ngbr_isnt_Void, Cell_isnt_Void)

    if (boussinesq_approximation) then
       where (fluidVof > 0.0_r8)
          fluidRho      = fluidRho / fluidVof
          fluidDeltaRho = fluidDeltaRho / fluidVof
       end where
    else
       where (fluidVof > 0.0_r8)
          fluidRho       = (fluidRho+fluidDeltaRho) / fluidVof
          fluidDeltaRho  = 0.0_r8
       end where
    endif

    ! Evaluate and store the minimum 'real' fluid density for evaluation of the
    ! face density limit
    Mask = RealFluidVof > 0.0_r8
    Prop = 0.0_r8
    where(Mask) Prop = fluidRho*fluidVof/RealFluidVof
    MinFluidRho = PGSLIB_GLOBAL_MINVAL(prop, MASK=Mask)

    ! Where cells are void (rho = 0.0_r8), zero out velocity & pressure
    ! NNC, Mar 2013.  Oh my! Why not just "mask = FluidRho == 0?
    Mask = .false.
    where (FluidRho == 0.0_r8 ) Mask = .true.
    do n = 1,ndim
       where (Mask) Zone%Vc(n) = 0.0_r8
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
          Zone(nc)%Vc(:)  = 0.0_r8
          Zone(nc)%Vc_old(:) = 0.0_r8
       endif
       do f= 1, nfc
          if(Solid_Face(f,nc)) then
             Fluxing_Velocity(f,nc)   = 0.0_r8
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
                  if(BC_Mat(f,c) /= NULL_I) then
                     if(DENSITY_MATERIAL(BC_Mat(f,c)) > 0.0_r8 ) then
                        SkipFlow = .false.
                        EXIT SKIPFLOWLOOP
                    endif
                  endif
              endif
              if(Boundary_Flag(f,c) == 2) then
              ! Dirichlet Velocity BC
                  if(BC_Mat(f,c) /= NULL_I) then
                     if(DENSITY_MATERIAL(BC_Mat(f,c)) > 0.0_r8 ) then
                        !Inward Flow?
                        !! NNC, Jan 2014. Time-dependent dirichlet velocity
                        !ORIG: BC_V_Dot_N = 0.0
                        !ORIG: do n = 1, ndim
                        !ORIG:    BC_V_Dot_N = BC_V_Dot_N + BC_Vel(n,f,c)*Cell(c)%Face_Normal(n,f)
                        !ORIG: enddo
                        BC_V_Dot_N = dot_product(bndry_vel%get(f,c,t), Cell(c)%Face_Normal(:,f))
                        if(BC_V_Dot_N < 0.0_r8) then
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

    use legacy_mesh_api, only: ncells, ndim, Cell
    use region_data, only: Regions, nregion

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

  subroutine compute_cell_property(prop, temp, value, fluid)

    use legacy_mesh_api, only: ncells
    use scalar_func_class
    use scalar_func_tools, only: is_const
    use matl_module, only: gather_vof

    character(*), intent(in) :: prop
    real(r8), intent(in) :: temp(:)
    real(r8), intent(out) :: value(:)
    logical, intent(in), optional :: fluid

    integer :: m, j
    real(r8) :: vofm (ncells), state(1), mval
    class(scalar_func), allocatable :: prop_fun
    logical :: fluids_only

    ASSERT(size(temp) == ncells)
    ASSERT(size(value) == ncells)

    fluids_only = .false.
    if (present(fluid)) fluids_only = fluid

    value = 0.0_r8
    do m = 1, matl_model%nphase_real
      if (fluids_only .and. .not.matl_model%is_fluid(m)) cycle
      call matl_model%alloc_phase_prop(m, prop, prop_fun)
      ASSERT(allocated(prop_fun))
      call gather_vof (m, vofm)
      if (is_const(prop_fun)) then
        mval = prop_fun%eval(state)  ! state is ignored, but required
        value = value + mval*vofm
      else
        do j = 1, ncells
          if (vofm(j) > 0.0_r8) then
            state(1) = temp(j)
            mval = prop_fun%eval(state)
            value(j) = value(j) + mval*vofm(j)
          end if
        end do
      end if
    end do

  end subroutine compute_cell_property

end module flow_property_module

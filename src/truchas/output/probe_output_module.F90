#include "f90_assert.fpp"

MODULE PROBE_OUTPUT_MODULE
  !-----------------------------------------------------------------------------
  ! Purpose:
  !
  !   Write output for all existing probes 
  !
  !   Public Interface:
  !
  !      call PROBES_OUTPUT ()
  !
  ! Contains: PROBES_OUTPUT
  !
  ! Author(s): Sharen Cummins (scummins@lanl.gov)
  !-----------------------------------------------------------------------------
  use kinds, only: r8
  implicit none
  private

  ! public procedures
  public :: PROBES_OUTPUT, PROBE_INIT
#ifdef USE_DANU
  public :: PROBE_INIT_DANU
#endif

  integer, Save, Public :: Probe_Output_Cycle_Multiplier = 1

CONTAINS

  SUBROUTINE PROBES_OUTPUT ()
    !---------------------------------------------------------------------------
    ! Purpose:
    
    ! For all probes we need to append probe field data to each respective 
    ! binary probe field data file.
    !
    ! Append in the following order:
    !       1. all scalar variables
    !       2. all vector variables
    !       3. all tensor variables
    !
    !---------------------------------------------------------------------------
    use probe_module,           only: probes
    use parameter_module,       only: nprobes
    use parallel_info_module,   only: p_info
    use time_step_module,       only: t2, cycle_number
#ifdef USE_TBROOK
    use output_data_module,     only: enable_tbrook_output
    use tbrook_module,          only: TBROOK_WRITE, TB_SCOPE_LOCAL
#endif
    use diagnostics_module,     only: PROBES_FIELDS
#ifdef USE_DANU
    use truchas_danu_output_data, only: sid
    use danu_module
    use parallel_communication, only: is_IOP, broadcast
#endif

    ! Local Variables
    integer :: iStatus, i, j, count, thesize
    real(r8), pointer, dimension(:) :: probe_cycle
#ifdef USE_DANU
    integer :: stat
    real(r8), allocatable :: probe_data(:,:)
#endif

    call PROBES_FIELDS ()

    iStatus = 0

    do i=1,nprobes

       probe_cycle  => NULL() !array containing scalar,vector or tensor probe field value for a given cycle

       count = 1
       do j=1,SIZE(probes(i)%ScalarVarLU)
          
          thesize = 1
          ALLOCATE(probe_cycle(thesize+2))
          
          probe_cycle(1)           = cycle_number
          probe_cycle(2)           = t2
          probe_cycle(3:thesize+2) = probes(i)%ScalarVarLU(j)%field

#ifdef USE_TBROOK
          if (enable_tbrook_output) then
          if ( iStatus == 0 .and. p_info%iop) then
             call TBrook_Write(probes(i)%BrookLU(count),probe_cycle,iStatus=iStatus,SCOPE=TB_SCOPE_LOCAL)
          end if
          end if
#endif
          
#ifdef USE_DANU
          ! Create the probe data set in the HDF file
          allocate(probe_data(size(probe_cycle),1))
          probe_data(:,1) = probe_cycle
          if (is_IOP) call probe_data_write (probes(i)%pid(count), probe_data, stat)
          call broadcast (stat)
          INSIST(stat == DANU_SUCCESS)
          deallocate(probe_data)
#endif

          if (ASSOCIATED(probe_cycle)) DEALLOCATE(probe_cycle)

          count = count + 1

       end do

       do j=1,SIZE(probes(i)%VectorVarLU)
          
          thesize = SIZE(probes(i)%VectorVarLU(1)%field)
          ALLOCATE(probe_cycle(thesize+2))
          
          probe_cycle(1)           = cycle_number
          probe_cycle(2)           = t2
          probe_cycle(3:thesize+2) = probes(i)%VectorVarLU(j)%field(:)

#ifdef USE_TBROOK
          if (enable_tbrook_output) then
          if ( iStatus == 0 .and. p_info%iop) then
             call TBrook_Write(probes(i)%BrookLU(count),probe_cycle,iStatus=iStatus,SCOPE=TB_SCOPE_LOCAL)
          end if
          end if
#endif
          
#ifdef USE_DANU
          ! Create the probe data set in the HDF file
          allocate(probe_data(size(probe_cycle),1))
          probe_data(:,1) = probe_cycle
          if (is_IOP) call probe_data_write (probes(i)%pid(count), probe_data, stat)
          call broadcast (stat)
          INSIST(stat == DANU_SUCCESS)
          deallocate(probe_data)
#endif

          if (ASSOCIATED(probe_cycle)) DEALLOCATE(probe_cycle)

          count = count + 1

       end do

       do j=1,SIZE(probes(i)%TensorVarLU)
          
          thesize = SIZE(probes(i)%TensorVarLU(1)%field)
          ALLOCATE(probe_cycle(thesize+2))
          
          probe_cycle(1)           = cycle_number
          probe_cycle(2)           = t2
          probe_cycle(3:thesize+2) = probes(i)%TensorVarLU(j)%field(:)

#ifdef USE_TBROOK
          if (enable_tbrook_output) then
          if ( iStatus == 0 .and. p_info%iop) then
             call TBrook_Write(probes(i)%BrookLU(count),probe_cycle,iStatus=iStatus,SCOPE=TB_SCOPE_LOCAL)
          end if
          end if
#endif

#ifdef USE_DANU
          ! Create the probe data set in the HDF file
          allocate(probe_data(size(probe_cycle),1))
          probe_data(:,1) = probe_cycle
          if (is_IOP) call probe_data_write (probes(i)%pid(count), probe_data, stat)
          call broadcast (stat)
          INSIST(stat == DANU_SUCCESS)
          deallocate(probe_data)
#endif
          
          if (ASSOCIATED(probe_cycle)) DEALLOCATE(probe_cycle)

          count = count + 1

       end do

    end do

  END SUBROUTINE PROBES_OUTPUT

  !-----------------------------------------------------------------------------


  SUBROUTINE PROBE_INIT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Initialize all quantities in the probe_variable structure
    !   (refer to setup/base_types/probe_module.F90 for this structure)
    !
    ! Author(s): Sharen Cummins (scummins@lanl.gov)
    !=======================================================================
    use parameter_module,           only: nprobes, nmat
    use solid_mechanics_data,       only: solid_mechanics
    use EM_data_proxy,              only: EM_is_on
    use input_utilities,            only: NULL_I
    use property_module,            only: GET_USER_MATERIAL_ID
    use probe_data_module
    use probe_module
#ifdef USE_TBROOK
    use output_data_module,         only: enable_tbrook_output
    use tbrook_module,              only: Brook
#endif

    ! Local Variables
    integer :: i, j, count, index,               &
               cellf_count, nodef_count,cellv_count, nodev_count, &
               cellt_count, nodet_count
    character(LEN=256) :: str
#ifdef USE_TBROOK
    type(Brook) :: RHOB, DELTA_RHOB, TEMPB, dT_dtB, GRADTB,           &
                   ENTHALPYB, PB, JOULE_PB, EPSDOTB, GAPDISPB,        &
                   GAPFORCEB, VCB, DISPLACEMENTB, SIGMAB, EPSILONB,   &
                   EPSTHERMB, EPSPCB
    type(Brook), pointer, dimension(:) :: VOF_Bs
#endif

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    ! A given probe contains :
    !
    ! 1. name                      [user supplied name of the probe used to designate the probe in the parser]
    ! 2. description               [user supplied optional description]
    ! 3. coordinates               [user supplied, defining the location of the probe]
    ! 4. coordinate scaling factor [user supplied, required if a mesh scaling factor is employed]
    !
    ! 5. nearest cell index and coordinates
    !               [calculated in subroutine PROBES_POSITIONS in diagnostics_module.F90]
    ! 6. nearest node index and coordinates
    !               [calculated in subroutine PROBES_POSITIONS in diagnostics_module.F9]
    ! 7. names of all field variables whose values will be calculated at the probe location ('NameLU')
    ! 8. brook streams ('BrookLU')
    !               [pointers to lookaside binary files, a file contains probe data for one field variable]
    ! 9. list of scalar field values at the probe location ('ScalarVarLU')
    !               [calculated in subroutine PROBES_POSITIONS in diagnostics_module.F90]
    ! 10.list of vector field values at the probe location ('VectorVarLU')
    !               [calculated in subroutine PROBES_POSITIONS in diagnostics_module.F90]
    ! 11.list of tensor field values at the probe location ('TensorVarLU')
    !               [calculated in subroutine PROBES_POSITIONS in diagnostics_module.F90]
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    !In assigning names to 'NameLU' and brooks to 'BrookLU' the following order is used:

    !Scalar probe fields are stored first in the following order:
    !    a) cell fields  : [Rho,Delta_Rho,Temp,dT_dt,|Grad(T)|,Enthalpy,P,PC,Vof,Joule,Epsdot]
    !    b) nodal fields : [GapDisp,GapForce]

    !Vector probe fields are stored next in the following order:
    !    a) cell fields  : [Vc]
    !    b) nodal fields : [Disp]

    !Tensor probe fields are stored next in the following order:
    !    a) cell fields  : [Elastic Stress, Total Strain, Thermal Strain, PC Strain]
    !    b) nodal fields : []

    !***********************************************************************************
    !
    !NOTE:This ordering must remain consistent with the ordering provided in the
    !diagnostics subroutine 'PROBES_FIELDS' in src/output/diagnostics_module.F90
    !
    !***********************************************************************************


    do i = 1, nprobes

      probes(i)%name           = probe_name(i)
      probes(i)%description    = probe_description(i)
      probes(i)%coords_scale   = probe_coords_scale(i)
      probes(i)%coords(:)      = probes(i)%coords_scale*probe_coords(:, i)
      probes(i)%node%index     = NULL_I
      probes(i)%node%coords(:) = probe_coords(:,i)
      probes(i)%cell%index     = NULL_I
      probes(i)%cell%coords(:) = probe_coords(:,i)

      !Store all names of scalar fields, vector fields and tensor fields for each probe in a NameLU table

      !first all scalar fields......
      !            cell centered fields.....
      probes(i)%NameLU(1)      = 'RHO'
      probes(i)%NameLU(2)      = 'DELTA_RHO'
      probes(i)%NameLU(3)      = 'TEMP'
      probes(i)%NameLU(4)      = 'dT_dt'
      probes(i)%NameLU(5)      = 'Grad(T)'
      probes(i)%NameLU(6)      = 'ENTHALPY'
      probes(i)%NameLU(7)      = 'P'

      count = 8

61    FORMAT('VOF',i4.4)
      do j=1,nmat
         write(str,61) GET_USER_MATERIAL_ID(j)
         probes(i)%NameLU(count) = str
         count = count + 1
      end do

      if (EM_is_on()) then
         probes(i)%NameLU(count)  = 'JOULE_P'
         count = count + 1
      end if
      if (solid_mechanics) then
         probes(i)%NameLU(count)  = 'EPSDOT'
         count = count + 1
      end if

      cellf_count                 = count-1               !number of cell centered scalar fields

      !node centered fields....
      if (solid_mechanics) then
         probes(i)%NameLU(count)  = 'GAPDISP'
         count = count + 1
         probes(i)%NameLU(count)  = 'GAPFORCE'
         count = count + 1
      end if
      nodef_count                 = count-1-cellf_count    !number of node centered scalar fields

      !now all vector fields......
      !cell centered vector fields...
      probes(i)%NameLU(count)    = 'VC'
      count = count + 1
      cellv_count                = count-1-(nodef_count+cellf_count) !number of cell centered vector fields

      !node centered vector fields...
      if (solid_mechanics) then
         probes(i)%NameLU(count) = 'DISPLACEMENT'
         count = count + 1
      end if
      nodev_count                = count-1-(cellv_count+ nodef_count+cellf_count) !number of node centered vector fields

      !now all tensor fields......
      !cell centered tensor fields....
      if (solid_mechanics) then
         probes(i)%NameLU(count)  = 'SIGMA'
         count = count + 1
         probes(i)%NameLU(count)  = 'EPSILON'
         count = count + 1
         probes(i)%NameLU(count)  = 'EPSTHERM'
         count = count + 1
         probes(i)%NameLU(count)  = 'EPSPC'
         count = count + 1
      end if
      cellt_count              = count-1-(nodef_count+cellf_count) - &
                                 (nodev_count+cellv_count)           !number of cell centered tensor fields

      !node centered tensor fields...
      nodet_count              = 0                                   !number of node centered tensor fields

#ifdef USE_TBROOK
      if (enable_tbrook_output) then
      !store all scalar field brooks, all vector field brooks and all tensor field brooks for each probe in a BrookLU table

      VOF_Bs       => NULL()

      !first all scalar field brooks...
      probes(i)%BrookLU(1)     = RHOB
      probes(i)%BrookLU(2)     = DELTA_RHOB
      probes(i)%BrookLU(3)     = TEMPB
      probes(i)%BrookLU(4)     = dT_dtB
      probes(i)%BrookLU(5)     = GRADTB
      probes(i)%BrookLU(6)     = ENTHALPYB
      probes(i)%BrookLU(7)     = PB

      count = 8

      index = 1
      ALLOCATE(VOF_Bs(nmat))
       do j=1,nmat
         probes(i)%BrookLU(count) = VOF_Bs(index)
         count = count + 1
         index = index + 1
      end do
      if (ASSOCIATED(VOF_Bs))       DEALLOCATE(VOF_Bs)

      if (EM_is_on()) then
         probes(i)%BrookLU(count)  = JOULE_PB
         count = count + 1
      end if
      if (solid_mechanics) then
         probes(i)%BrookLU(count)  = EPSDOTB
         count = count + 1
      end if

      if (solid_mechanics) then
         probes(i)%BrookLU(count)  = GAPDISPB
         count = count + 1
         probes(i)%BrookLU(count)  = GAPFORCEB
         count = count + 1
      end if
      probes(i)%BrookLU(count)  = VCB
      count = count + 1
      if (solid_mechanics) then
         probes(i)%BrookLU(count)  = DISPLACEMENTB
         count = count + 1
         probes(i)%BrookLU(count)  = SIGMAB
         count = count + 1
         probes(i)%BrookLU(count)  = EPSILONB
         count = count + 1
         probes(i)%BrookLU(count)  = EPSTHERMB
         count = count + 1
         probes(i)%BrookLU(count)  = EPSPCB
         count = count + 1
      end if
      end if
#endif

      !initialise all scalar probe variables for each probe
      !cell centered variables...
      do j=1,cellf_count
         probes(i)%ScalarVarLU(j)%field     = 0.0
         probes(i)%ScalarVarLU(j)%meshspace = 'cell'
      end do
      !node centered variables...
      do j=1,nodef_count
         probes(i)%ScalarVarLU(j+cellf_count)%field     = 0.0
         probes(i)%ScalarVarLU(j+cellf_count)%meshspace = 'node'
      end do

      !initialise all vector probe variables for each probe
      !cell centered variables...
      do j=1,cellv_count
         probes(i)%VectorVarLU(j)%field(:)  = 0.0
         probes(i)%VectorVarLU(j)%meshspace = 'cell'
      end do
      !node centered variables
      do j=1,nodev_count
         probes(i)%VectorVarLU(cellv_count+j)%field(:)  = 0.0
         probes(i)%VectorVarLU(cellv_count+j)%meshspace = 'node'
      end do

      !initialise all tensor probe variables for each probe
      !cell centered variables...
      do j=1,cellt_count
         probes(i)%TensorVarLU(j)%field(:)  = 0.0
         probes(i)%TensorVarLU(j)%meshspace = 'cell'
      end do
      do j=1,nodet_count
         probes(i)%TensorVarLU(j+cellt_count)%field(:)  = 0.0
         probes(i)%TensorVarLU(j+cellt_count)%meshspace = 'node'
      end do

      !Each brook in BrookLU is prepared for output in the tbrook_utility 'TBU_WRITEPROBES' subroutine
      !probes(i)%ScalarVarLU, VectorVarLU, TensorVarLU recalculated
      !at each cycle in the diagnostics_module subroutine 'PROBES_FIELDS'

    enddo

  END SUBROUTINE PROBE_INIT

#ifdef USE_DANU
  ! A version of TBU_WRITEPROBES that writes to HDF5 output instead
  SUBROUTINE PROBE_INIT_DANU

    use probe_module,           only: probes
    use parameter_module,       only: nprobes
    use time_step_module,       only: t, cycle_number
    use diagnostics_module,     only: PROBES_POSITIONS                   
    use truchas_danu_output_data, only: sid
    use danu_module
    use parallel_communication, only: is_IOP, broadcast
    use,intrinsic :: iso_c_binding, only: c_ptr 

    ! Local Variables 
    integer :: i, j, count, scalarsize, vectorsize, tensorsize, stat
    real(r8), pointer, dimension(:) :: probe_cycle, probe_cycleV, probe_cycleT
    type(c_ptr) :: pid
    character(64) :: dataset
    real(r8), allocatable :: probe_data(:,:)

    call PROBES_POSITIONS ()

    do i=1,nprobes

       probe_cycle  => NULL() !array containing cycle, time, scalar probe value for a given cycle
       probe_cycleV => NULL() !array containing cycle, time, vector probe values for a given cycle
       probe_cycleT => NULL() !array containing cycle, time, tensor probe values for a given cycle

       scalarsize        = 3
       ALLOCATE(probe_cycle(scalarsize))
       probe_cycle       = 0.0
       probe_cycle(1)    = cycle_number
       probe_cycle(2)    = t

       vectorsize        = SIZE(probes(i)%VectorVarLU(1)%field) + 2
       ALLOCATE(probe_cycleV(vectorsize))
       probe_cycleV      = 0.0
       probe_cycleV(1:2) = probe_cycle(1:2)

       tensorsize        = SIZE(probes(i)%TensorVarLU(1)%field) + 2
       ALLOCATE(probe_cycleT(tensorsize))
       probe_cycleT      = 0.0
       probe_cycleT(1:2) = probe_cycle(1:2)

       count = 1
       do j=1,SIZE(probes(i)%ScalarVarLU)
          probe_cycle(3) = probes(i)%ScalarVarLU(j)%field
          ! Create the probe data set in the HDF file
          allocate(probe_data(size(probe_cycle),1))
          probe_data(:,1) = probe_cycle
          write(dataset,'(a,i3.3,2a)') 'P', i, ':', trim(adjustl(probes(i)%NameLU(count)))
          if (is_IOP) call probe_create_data (sid, trim(dataset), probe_data, probes(i)%pid(count), stat)
          call broadcast (stat)
          INSIST(stat == DANU_SUCCESS)
          deallocate(probe_data)
          if (is_IOP) call probe_data_open (sid, dataset, pid, stat)
          call broadcast (stat)
          INSIST(stat == DANU_SUCCESS)
          if (is_IOP) then
            call attribute_write (pid, 'NAME', trim(probes(i)%name), stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'X', probes(i)%coords(1), stat)
            call attribute_write (pid, 'Y', probes(i)%coords(2), stat)
            call attribute_write (pid, 'Z', probes(i)%coords(3), stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'DESCRIPTION', trim(probes(i)%description), stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'CELL_INDEX', probes(i)%cell%index, stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'CELL_X', probes(i)%cell%coords(1), stat)
            call attribute_write (pid, 'CELL_Y', probes(i)%cell%coords(2), stat)
            call attribute_write (pid, 'CELL_Z', probes(i)%cell%coords(3), stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'NODE_INDEX', probes(i)%node%index, stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'NODE_X', probes(i)%node%coords(1), stat)
            call attribute_write (pid, 'NODE_Y', probes(i)%node%coords(2), stat)
            call attribute_write (pid, 'NODE_Z', probes(i)%node%coords(3), stat)
            INSIST(stat == DANU_SUCCESS)
          end if
          count = count + 1
       end do

       do j=1,SIZE(probes(i)%VectorVarLU)
          probe_cycleV(3:vectorsize) = probes(i)%VectorVarLU(j)%field
          ! Create the probe data set in the HDF file
          allocate(probe_data(size(probe_cycleV),1))
          probe_data(:,1) = probe_cycleV
          write(dataset,'(a,i3.3,2a)') 'P', i, ':', trim(adjustl(probes(i)%NameLU(count)))
          if (is_IOP) call probe_create_data (sid, trim(dataset), probe_data, probes(i)%pid(count), stat)
          call broadcast (stat)
          INSIST(stat == DANU_SUCCESS)
          deallocate(probe_data)
          if (is_IOP) call probe_data_open (sid, dataset, pid, stat)
          call broadcast (stat)
          INSIST(stat == DANU_SUCCESS)
          if (is_IOP) then
            call attribute_write (pid, 'NAME', trim(probes(i)%name), stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'X', probes(i)%coords(1), stat)
            call attribute_write (pid, 'Y', probes(i)%coords(2), stat)
            call attribute_write (pid, 'Z', probes(i)%coords(3), stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'DESCRIPTION', trim(probes(i)%description), stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'CELL_INDEX', probes(i)%cell%index, stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'CELL_X', probes(i)%cell%coords(1), stat)
            call attribute_write (pid, 'CELL_Y', probes(i)%cell%coords(2), stat)
            call attribute_write (pid, 'CELL_Z', probes(i)%cell%coords(3), stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'NODE_INDEX', probes(i)%node%index, stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'NODE_X', probes(i)%node%coords(1), stat)
            call attribute_write (pid, 'NODE_Y', probes(i)%node%coords(2), stat)
            call attribute_write (pid, 'NODE_Z', probes(i)%node%coords(3), stat)
            INSIST(stat == DANU_SUCCESS)
          end if
          count = count + 1
       end do

       do j=1,SIZE(probes(i)%TensorVarLU)
          probe_cycleT(3:tensorsize) = probes(i)%TensorVarLU(j)%field
          ! Create the probe data set in the HDF file
          allocate(probe_data(size(probe_cycleT),1))
          probe_data(:,1) = probe_cycleT
          write(dataset,'(a,i3.3,2a)') 'P', i, ':', trim(adjustl(probes(i)%NameLU(count)))
          if (is_IOP) call probe_create_data (sid, trim(dataset), probe_data, probes(i)%pid(count), stat)
          call broadcast (stat)
          INSIST(stat == DANU_SUCCESS)
          deallocate(probe_data)
          if (is_IOP) call probe_data_open (sid, dataset, pid, stat)
          call broadcast (stat)
          INSIST(stat == DANU_SUCCESS)
          if (is_IOP) then
            call attribute_write (pid, 'NAME', trim(probes(i)%name), stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'X', probes(i)%coords(1), stat)
            call attribute_write (pid, 'Y', probes(i)%coords(2), stat)
            call attribute_write (pid, 'Z', probes(i)%coords(3), stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'DESCRIPTION', trim(probes(i)%description), stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'CELL_INDEX', probes(i)%cell%index, stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'CELL_X', probes(i)%cell%coords(1), stat)
            call attribute_write (pid, 'CELL_Y', probes(i)%cell%coords(2), stat)
            call attribute_write (pid, 'CELL_Z', probes(i)%cell%coords(3), stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'NODE_INDEX', probes(i)%node%index, stat)
            INSIST(stat == DANU_SUCCESS)
            call attribute_write (pid, 'NODE_X', probes(i)%node%coords(1), stat)
            call attribute_write (pid, 'NODE_Y', probes(i)%node%coords(2), stat)
            call attribute_write (pid, 'NODE_Z', probes(i)%node%coords(3), stat)
            INSIST(stat == DANU_SUCCESS)
          end if
          count = count + 1
       end do

       if (ASSOCIATED(probe_cycle))  DEALLOCATE(probe_cycle)
       if (ASSOCIATED(probe_cycleV)) DEALLOCATE(probe_cycleV)
       if (ASSOCIATED(probe_cycleT)) DEALLOCATE(probe_cycleT)

    end do

  END SUBROUTINE PROBE_INIT_DANU
#endif

END MODULE PROBE_OUTPUT_MODULE

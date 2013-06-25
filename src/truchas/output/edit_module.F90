MODULE EDIT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define quantities and procedures for long and short edits
  !   to the "prefix.out" output file.
  !
  !   Public Interface:
  !
  !     * call EDIT_SHORT ()
  !
  !         Perform a "short edit" by computing and printing global
  !         diagnostics.
  !
  !     * call EDIT_LONG ()
  !
  !         Perform a "long edit" by printing all components of all
  !         base derived types in user-specified cells.
  !
  ! Contains: EDIT_LONG
  !           EDIT_LONG_INIT
  !           EDIT_SHORT
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use long_edit_data_types
  use kind_module,      only: int_kind, log_kind, real_kind
  use parameter_module, only: max_long_edit_boxes, max_long_edit_cells, &
                              mops, ndim
  use truchas_logging_services
  implicit none

  ! Private Module
  private

  ! Public Procedures
  public :: EDIT_LONG, EDIT_SHORT

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! OUTPUTS namelist edit variables.
  integer(KIND = int_kind), dimension(mops), &
                            save, public :: Long_Output_Dt_Multiplier
  integer(KIND = int_kind), dimension(mops), &
                            save, public :: Short_Output_Dt_Multiplier

  real(KIND = real_kind), dimension(2,ndim,max_long_edit_boxes), &
                          save, public :: Long_Edit_Bounding_Coords

  ! Edit flags.
  logical(KIND = log_kind), save, public :: long_edit
  logical(KIND = log_kind), save, public :: short_edit

  ! Long edit variables.
  integer(KIND = int_kind), save, public :: long_edit_cells
  integer(KIND = int_kind), dimension(max_long_edit_cells), &
                            save, public :: Long_Edit_Cell_List

  integer(int_kind), save :: num_local_long_edits
  integer(int_kind), save :: tot_local_long_edits
  
  type (LONG_EDIT_LIST), target, save, public :: Local_Long_Edit_List
  type (LONG_EDIT_LIST), target, save, public :: Collated_Long_Edit_List
  type (MECH_EDIT_LIST), target, save, public :: Local_Mech_Edit_List
  type (MECH_EDIT_LIST), target, save, public :: Collated_Mech_Edit_List

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE EDIT_LONG ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Perform a "long edit" by printing all components of all
    !   base derived types in specified cells.
    !
    !=======================================================================
    use kind_module,      only: int_kind, log_kind
    use matl_module,      only: Matl
    use mesh_module,      only: UnPermute_Mesh_Vector
    use truchas_env,      only: output_file_name
    use parameter_module, only: mat_slot
    use time_step_module, only: cycle_number, t
    use zone_module,      only: Zone
    use output_utilities, only: ANNOUNCE_FILE_WRITE
    use solid_mechanics_module, only: STRESS_STRAIN_INVARIANTS
    use solid_mechanics_data, only:  solid_mechanics, Cell_Mech_Invariant
    use property_module,  only: Get_User_Material_ID
    use gap_output,   only: SET_GAP_ELEMENT_OUTPUT

    implicit none

    ! Argument List

    ! Local Variables
    logical(KIND = log_kind), save :: first_time = .true.
    logical (log_kind), save :: no_long_edit_cells
    integer(KIND = int_kind) :: i, j, s
    integer                  :: status

    type(CELL_AVG), POINTER :: zone_info => NULL()
    type(MATERIAL), POINTER :: matl_info => NULL()
    ! Solid mechanics data
    type(CELL_MECH_INVARIANT), pointer, dimension(:) :: mech_info => NULL()
    type(CELL_MECH_INVARIANT), pointer               :: mech_item => NULL()
    integer(int_kind) :: cell_number
    character(128) :: string1, string2

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Write out the leading header.
    write (string1, 10) t, cycle_number
10  format (1x,14('* '),'LONG EDIT: t = ',1pe13.5,', step = ',i6,14(' *'))
    call TLS_info ('')
    call TLS_info (string1)
    call TLS_info ('')

    ! If this is the first time here, construct the long edit cell list.
    if (first_time) then
       call EDIT_LONG_INIT (no_long_edit_cells)
       first_time = .false.
    end if

    if (no_long_edit_cells) return

    ! Write out the header.
    call TLS_info ('')
    call TLS_info (repeat(' ',37) // 'Fluid and Thermodynamic Information')
    call TLS_info ('')
    call TLS_info ('')
    call TLS_info (repeat(' ',82) // 'Matl Matl    Matl')
    call TLS_info ('   Cell   Density   Temperature   Pressure' // repeat(' ',17) // &
        'Velocity' // repeat(' ',15) // 'Slot  Id     VOF')
    call TLS_info ('   ----   -------   -----------   --------' // repeat(' ',17) // &
        '--------' // repeat(' ',15) // '----  --     ---')

    ! Loop over all long edit cells, printing out
    ! the Zone and Matl base types.
    LONG_EDIT_CELL_LOOP: do i = 1,long_edit_cells

       ! Get the cell index.
       j = Long_Edit_Cell_List(i)

       do s = 1,mat_slot
          if (s == 1) then
             ! For first material put in zone information and cell number.
             ! The original cell number is actually the integer
             ! in UnPermute_Mesh_Vector
             call SET(LOCAL_LONG_EDIT_LIST, ITEM = i,     &
                      CELLNUMBER=UnPermute_Mesh_Vector(j), &
                      ZONE = Zone(j))

          endif
          call SET(LOCAL_LONG_EDIT_LIST, ITEM = i,     &
                   MATL_SLOT = s,                       &
                   MATL = Matl(s)%Cell(j))

       end do

    end do LONG_EDIT_CELL_LOOP

    ! COLLATE the long edit lists
    call COLLATE(COLLATED_LONG_EDIT_LIST, LOCAL_LONG_EDIT_LIST)

    ! Now sort the long edit list
    call SORT(COLLATED_LONG_EDIT_LIST)

    ! Now we can output the list, and it will be all shown in the proper order
    LONG_EDIT_OUTPUT_LOOP: DO i = 1, SIZE(COLLATED_LONG_EDIT_LIST)
       do s = 1,mat_slot
          matl_info     => MATL_DATA (LIST = COLLATED_LONG_EDIT_LIST, ITEM = i, MATL_SLOT = s)
          if (s == 1) then ! First material print out lots of extra info
             cell_number   =  CELLNUMBER(LIST = COLLATED_LONG_EDIT_LIST, ITEM = i)
             zone_info     => ZONE_DATA (LIST = COLLATED_LONG_EDIT_LIST, ITEM = i)
             write(string1,'(i6,6es12.4)') cell_number, Zone_info%Rho, Zone_info%Temp, Zone_info%P, Zone_info%Vc(1:3)
          else
             string1 = ''
          endif
          write(string2,'(2i6,es12.4)') s, Get_User_Material_ID(Matl_info%ID), Matl_info%Vof
          call TLS_info (string1(:78) // string2)
       end do
    end do LONG_EDIT_OUTPUT_LOOP       

    ! If solid mechanics, write out data.

    SOLID_MECHANICS_LOOP: if (solid_mechanics) then
       ! Put reasonable values in gap_elements
       call SET_GAP_ELEMENT_OUTPUT()
       ! Calculate or update the scalar stress and strain invariants
       mech_info => STRESS_STRAIN_INVARIANTS()
       ! Write out the solid mechanics header.
       call TLS_info ('')
       call TLS_info ('               Solid Mechanics Information')
       call TLS_info ('')
       call TLS_info ('                     Effective                Total')
       call TLS_info ('         von Mises    Plastic      Mean     Volumetric')
       call TLS_info ('   Cell   Stress      Strain      Stress      Strain')
       call TLS_info ('   ----  --------    --------    --------    --------')
       
       MECH_SET_LOOP: do i = 1,long_edit_cells
          
          ! Get the cell index.
          j = Long_Edit_Cell_List(i)
          mech_item => mech_info(j)
          call SET(LOCAL_MECH_EDIT_LIST, i, UnPermute_Mesh_Vector(j), mech_item)
          
       end do MECH_SET_LOOP

       ! COLLATE the long edit lists
       call COLLATE(COLLATED_MECH_EDIT_LIST, LOCAL_MECH_EDIT_LIST)
       
       ! Now sort the long edit list
       call SORT(COLLATED_MECH_EDIT_LIST)

       ! Now we can output the list, and it will be all shown in the proper order       
       MECH_OUT_LOOP: DO i = 1, SIZE(COLLATED_MECH_EDIT_LIST)

          cell_number   =  CELLNUMBER(COLLATED_MECH_EDIT_LIST, i)
          mech_item     => LIST_ITEM_MECH (COLLATED_MECH_EDIT_LIST, i)
          
          write(string1,'(i6,4es12.4)') cell_number, mech_item%mises_stress, &
              mech_item%eff_plastic_strain, mech_item%mean_stress, mech_item%volumetric_strain
          call TLS_info (string1)

       end DO MECH_OUT_LOOP

       DEALLOCATE(mech_info)

    end if SOLID_MECHANICS_LOOP


    ! Write out the trailing header.
    write (string1, 10) t, cycle_number
    call TLS_info ('')
    call TLS_info (string1)
    call TLS_info ('')
    call ANNOUNCE_FILE_WRITE ('Long edit', output_file_name('log'))

999 continue

    return

  END SUBROUTINE EDIT_LONG

  SUBROUTINE EDIT_LONG_INIT (no_long_edit_cells)
    !=======================================================================
    ! Purpose:
    !
    !   Perform the first-time portion of a "long edit" by printing
    !   all components of the Mesh, Cell, and Vertex derived types.
    !=======================================================================
    use ArrayAllocate_Module, only: ARRAYCREATE, ARRAYDESTROY
    use gs_module,            only: EN_GATHER, EE_Gather
    use kind_module,          only: int_kind, log_kind, real_kind
    use mesh_module,          only: Cell, Mesh, Vertex, Vrtx_Bdy, &
                                    UnPermute_Mesh_Vector,        &
                                    UnPermute_Vertex_Vector
    use parameter_module,     only: max_long_edit_boxes, max_long_edit_cells, &
                                    ncells, ncells_tot, ndim, nfc, nvc
    use parallel_util_module, only: Is_IO_PE
    use pgslib_module,        only: PGSLIB_GLOBAL_ALL, PGSLIB_GLOBAL_SUM,     &
                                    PGSLIB_GLOBAL_MINVAL, PGSLIB_SUM_PREFIX,  &
                                    PGSLIB_PERMUTE, PGSLIB_GATHER, PGSLIB_GRADE_UP
    use solid_mechanics_data, only: solid_mechanics
    use zone_module,          only: CELL_AVG
    use matl_module,          only: MATERIAL

    implicit none

    ! Argument List
    logical(log_kind), intent(OUT) :: no_long_edit_cells

    ! Local Variables
    logical(log_kind) :: long_edit_cell
    integer(int_kind) :: b, f, i, j, n, v, cell_number, long_edit_cells_old, s
    real(real_kind)   :: distance
    type(CELL_GEOMETRY), pointer                      :: Cell_Info => NULL()
    type(VERTEX_DATA),   pointer                      :: Vrtx_Info => NULL()
    type(FACE_DATA),     pointer                      :: Face_Info => NULL()
    type(CELL_AVG)                                    :: zone_zero
    type(MATERIAL)                                    :: matl_zero
    logical(log_kind), dimension(max_long_edit_boxes) :: Long_Edit_Point
    logical(log_kind), dimension(ncells)              :: Long_Edit_Cell_Mask
    logical(log_kind), dimension(ncells)              :: Long_Edit_Cell_Mask_Ordered
    integer(int_kind), dimension(:,:), pointer        :: Mesh_Ngbr_Vrtx_Orig_Order => NULL()
    integer(int_kind), dimension(:,:), pointer        :: Mesh_Ngbr_Cell_Orig_Order => NULL()
    integer(int_kind), dimension(max_long_edit_boxes) :: Distance_Min_Loc = 0
    integer(int_kind), dimension(ncells)              :: Long_Edit_Cells_Orig_Number
    integer(int_kind), dimension(ncells)              :: Long_Edit_Cells_Ordered
    integer(int_kind), dimension(ncells)              :: Cell_Rank
    integer(int_kind), dimension(ncells)              :: Global_Cell_Count
    real(real_kind),   dimension(max_long_edit_boxes) :: Distance_Min     = HUGE(1_real_kind)
    real(real_kind), pointer, dimension(:,:,:)        :: Xv => NULL()
    character(LEN = 256) :: tmp_string
    character(LEN = 256) :: tmp_string2

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Check to see if we have any long edit points. Set flag appropriately.
    LONG_EDIT_POINTS: do b = 1,max_long_edit_boxes
       Long_Edit_Point(b) = .true.
       do n = 1,ndim
          Long_Edit_Point(b) = Long_Edit_Point(b) .and. &
                              (Long_Edit_Bounding_Coords(1,n,b) == Long_Edit_Bounding_Coords(2,n,b))
       end do
    end do LONG_EDIT_POINTS

    ! assume there are some long edit cells
    no_long_edit_cells = .false.
    long_edit_cell_list = ncells_tot + 1

    ! Construct the list of long edit cells
    ! Note: Long_Edit_Bounding_Coords(1:2,n) = (min,max) for coordinate n.
    LONG_EDIT_LIST: do i = 1,ncells
       LONG_EDIT_BOXES: do b = 1,max_long_edit_boxes
          if (Long_Edit_Point(b)) then
             long_edit_cell = .false.
             distance = 0_real_kind
             do n = 1,ndim
                distance = distance + (Cell(i)%Centroid(n) - Long_Edit_Bounding_Coords(1,n,b))**2
             end do
             distance = SQRT(distance)
             Distance_Min(b) = MIN(Distance_Min(b),distance)
             if (Distance_Min(b) == distance) Distance_Min_Loc(b) = i
          else
             long_edit_cell = .true.
             do n = 1,ndim
                long_edit_cell = long_edit_cell .and. &
                     Cell(i)%Centroid(n) >= Long_Edit_Bounding_Coords(1,n,b) .and. &
                     Cell(i)%Centroid(n) <= Long_Edit_Bounding_Coords(2,n,b)
             end do
          end if
          if (long_edit_cell) exit LONG_EDIT_BOXES
       end do LONG_EDIT_BOXES
       if (long_edit_cell .and. long_edit_cells < max_long_edit_cells) then
          long_edit_cells = long_edit_cells + 1
          Long_Edit_Cell_List(long_edit_cells) = i
       end if
    end do LONG_EDIT_LIST

    ! If we have any long edit points, save only those cells 
    ! whose distance is minimum to that point.
    LONG_EDIT_LOC: do b = 1,max_long_edit_boxes
       if (Long_Edit_Point(b)) then
          distance = PGSLIB_GLOBAL_MINVAL(Distance_Min(b))
          if (distance == Distance_Min(b) .and. long_edit_cells < max_long_edit_cells) then
             long_edit_cells = long_edit_cells + 1
             Long_Edit_Cell_List(long_edit_cells) = Distance_Min_Loc(b)
          end if
       end if
    end do LONG_EDIT_LOC

    ! Now make sure we don't have more long edit cells than max_long_edit_cells, globally.
    ! If we do, then take the first max_long_edit_cells based on global cell number.
    ! Sort cell list
    Long_Edit_Cells_Orig_Number = ncells_tot + 1
    do i = 1, long_edit_cells
       Long_Edit_Cells_Orig_Number(i) = UnPermute_Mesh_Vector(Long_Edit_Cell_List(i))
    end do
    Cell_Rank = PGSLIB_GRADE_UP (Long_Edit_Cells_Orig_Number)
    call PGSLIB_PERMUTE (DEST   = Long_Edit_Cells_Ordered,     &
                         SOURCE = Long_Edit_Cells_Orig_Number, &
                         INDEX  = Cell_Rank)

    ! Pick first max_long_edit_cells
    ! Get a global counter
    Global_Cell_Count = PGSLIB_SUM_PREFIX ((/ (1, i=1,ncells) /))
    ! Now pick only the first max_long_edit_cells, and also skip empty slots
    ! in Long_edit_Cell_List.
    Long_Edit_Cell_Mask_Ordered = (Global_Cell_Count <= Max_Long_Edit_Cells) .and. &
                                  (Long_Edit_Cells_Ordered <= ncells_tot)
    ! "Unpermute" the mask back into the local order (where the long edit cells were found)
    call PGSLIB_GATHER (DEST   = Long_Edit_Cell_Mask,         &
                        SOURCE = Long_Edit_Cell_Mask_Ordered, &
                        INDEX  = cell_rank)

    ! Pack the data in Long_Edit_Cell_List based on the mask we just found.
    ! Also, get a new count of how many long edit cells we are keeping.
    long_edit_cells_old = long_edit_cells
    long_edit_cells     = 0
    do i = 1, long_edit_cells_old
       if (Long_Edit_Cell_Mask(i)) then
          long_edit_cells = long_edit_cells + 1
          Long_Edit_Cell_List(long_edit_cells) = Long_Edit_Cell_List(i)
       else
          cycle
       end if
    end do

    ! If no long edit cells were found globally, print warning and
    ! disable the long edit flag.
    if (PGSLIB_GLOBAL_ALL(long_edit_cells == 0)) then
       call TLS_warn ('No long edit cells were found - disabling all long edits!')
       Long_Output_Dt_Multiplier = 0
       long_edit = .false.
       no_long_edit_cells = .true.
       return
    end if

    ! Write out a header.
    write (tmp_string, 15)
15  format (45x,'Geometry and Connectivity Information')
    call TLS_info ('')
    call TLS_info (tmp_string)
    call TLS_info ('')

    ! Allocate arrays to hold vertex and face centroid coordinates.
    call ARRAYCREATE (Xv, 1, ndim, 1, nvc, 1, ncells, 'Array Xv(ndim,nvc,ncells)')

    ! Gather vertex coordinates.
    do n = 1,ndim
       call EN_GATHER (Xv(n,:,:), Vertex%Coord(n), BOUNDARY=Vrtx_Bdy(n)%Data)
    end do

    ! To ensure same output no matter what reordering is done, output original vertex and cell numbers
    call ARRAYCREATE (mesh_ngbr_vrtx_orig_order, 1, nvc, 1, ncells, 'Array mesh_ngbr_vrtx_orig_order(nvc,ncells)')
    call ARRAYCREATE (mesh_ngbr_cell_orig_order, 1, nfc, 1, ncells, 'Array mesh_ngbr_cell_orig_order(nfc,ncells)')
    call en_gather(mesh_ngbr_vrtx_orig_order, UnPermute_Vertex_Vector)
    call ee_gather(mesh_ngbr_cell_orig_order, UnPermute_Mesh_Vector)

    ! Allocate the long edit lists 
    NUM_LOCAL_LONG_EDITS = long_edit_cells
    TOT_LOCAL_LONG_EDITS = PGSLib_GLOBAL_SUM(NUM_LOCAL_LONG_EDITS)
    if (.NOT. Is_IO_PE()) TOT_LOCAL_LONG_EDITS = 0

    call CREATE(LOCAL_LONG_EDIT_LIST, ITEMS=NUM_LOCAL_LONG_EDITS)
    call CREATE(COLLATED_LONG_EDIT_LIST, ITEMS= TOT_LOCAL_LONG_EDITS)
    if (solid_mechanics) then
       call CREATE(LOCAL_MECH_EDIT_LIST, ITEMS=NUM_LOCAL_LONG_EDITS)
       call CREATE(COLLATED_MECH_EDIT_LIST, ITEMS= TOT_LOCAL_LONG_EDITS)
    end if

    ! Set up a single zone entry to initialize the zone portion of LOCAL_LONG_EDIT_LIST
    zone_zero%Rho          = 0.0
    zone_zero%Rho_Old      = 0.0
    zone_zero%Temp         = 0.0
    zone_zero%Temp_Old     = 0.0
    zone_zero%Enthalpy     = 0.0
    zone_zero%Enthalpy_Old = 0.0
    zone_zero%P            = 0.0
    zone_zero%Vc           = 0.0
    zone_zero%Vc_Old       = 0.0

    ! Set up a single matl entry to initialize the matl portion of LOCAL_LONG_EDIT_LIST
    matl_zero%Id      = 0
    matl_zero%Vof     = 0.0
    matl_zero%Vof_Old = 0.0

    ! Loop over all long edit cells and print out the base types
    ! housing connectivity and geometric information.
    LONG_EDIT_CELL_LOOP: do i = 1,long_edit_cells

       ! Get the cell index.
       j = Long_Edit_Cell_List(i)

       ! Put the geometric information into the long edit structures
       call SET(LOCAL_LONG_EDIT_LIST, ITEM = i,      &
                CELLNUMBER= UnPermute_Mesh_Vector(j), &
                CELL_GEO = CELL(j))

       ! Initialize the zone portion to 0.0
       call SET(LOCAL_LONG_EDIT_LIST, ITEM = i,      &
                CELLNUMBER= UnPermute_Mesh_Vector(j), &
                ZONE = zone_zero)




       do s = 1,mat_slot
          call SET(LOCAL_LONG_EDIT_LIST, ITEM = i,     &
                   CELLNUMBER= UnPermute_Mesh_Vector(j), &
                   MATL_SLOT = s,                       &
                   MATL = matl_zero)
       end do


       ! Write out the cell geometry and vertex information.
       do v = 1,nvc

          ! Put vertex geometrix information into the long edit structures
          call SET(LOCAL_LONG_EDIT_LIST, ITEM = i,        &
                   VERTEX_NUMBER = v,                   &
                   VERTEX        = VERTEX_DATA(Mesh_ngbr_vrtx_orig_order(v,j), Xv(:,v,j) ) )
       end do

       ! Write out face-centered connectivity and geometry information.
       do f = 1,nfc
          ! Put face neighbor info into the long edit structures
          call SET(LOCAL_LONG_EDIT_LIST, ITEM = i,                                  &
                   FACE_NUMBER = f,                                                 &
                   FACE = FACE_DATA(mesh_ngbr_cell_orig_order(f,j),                 &
                                    Mesh(j)%Ngbr_Face(f),                           &
                                    (/ (Cell(j)%Face_Normal(n,f), n = 1,ndim) /),   &
                                    (/ (cell(j)%face_centroid(n,f), n = 1,ndim) /), &
                                    Cell(j)%Face_Area(f)) )
       end do

    end do LONG_EDIT_CELL_LOOP

    ! COLLATE the long edit lists
    call COLLATE(COLLATED_LONG_EDIT_LIST, LOCAL_LONG_EDIT_LIST)

    ! Now sort the long edit list
    call SORT(COLLATED_LONG_EDIT_LIST)

    ! Now we can output the list, and it will be all shown in the proper order

    DO i = 1, SIZE(COLLATED_LONG_EDIT_LIST)

       write(tmp_string,*) '  Cell     ' ,              &
                           '        Centroid          ', &
                           '           Volume' ,         &
                           '     Vertex Number ',        & 
                           '             Coordinates '
       write(tmp_string2,*)'   ----     ',               &
                           '        --------          ', &
                           '           ------',          & 
                           '     ------  ----- ',        &
                           '             ----------- ' 
       call TLS_info (tmp_string)
       call TLS_info (tmp_string2)
       do v = 1, nvc
          vrtx_info   => VERTEX_GEO_DATA(LIST = COLLATED_LONG_EDIT_LIST, ITEM = i, VERTEX_NUMBER = v)
          if (v == 1) then ! First vertex print out cell info also
             cell_number =  CELLNUMBER(LIST = COLLATED_LONG_EDIT_LIST, ITEM = i)
             cell_info   => CELL_GEO_DATA(LIST = COLLATED_LONG_EDIT_LIST, ITEM = i)
             write(tmp_string,'(i6,"   (",3es12.4,")",es12.4)') cell_number, Cell_info%centroid(1:3), Cell_info%Volume
          else
             tmp_string = ''
          end if
          write(tmp_string2,'(2i6,"   (",3es12.4,")")') v, vrtx_info%ngbr_vrtx, vrtx_info%Vertex_Coord(1:3)
          call TLS_info (tmp_string(1:59) // tmp_string2)
       end do
       call TLS_info ('')

       ! Output face neighbor information
       call TLS_info ('    Neighbor')
       tmp_string='   Face  Cell  Face           Face Unit Normal                       Face Centroid                Face Area'
       tmp_string2='   ----  ----  ----           ----------------                       -------------                ---------'
       call TLS_info (tmp_string)
       call TLS_info (tmp_string2)
       do f = 1, nfc
          face_info => FACE_GEO_DATA(LIST = COLLATED_LONG_EDIT_LIST, ITEM = i, FACE_NUMBER = f)
          write(tmp_string,'(3i6," (",3es12.4,") (",3es12.4,")",1es12.4)') &
              f, face_info%Ngbr_Cell, face_info%Ngbr_Face, &
              face_info%Face_Normal(1:3), face_info%Face_Centroid(1:3), face_info%Face_Area
          call TLS_info (tmp_string)
       end do
       call TLS_info ('')
    end DO

    ! Deallocate the coordinate temporaries.
    call ARRAYDESTROY (Xv, 'Array Xv(ndim,nvc,ncells)')
    call ARRAYDESTROY (mesh_ngbr_vrtx_orig_order, 'Array mesh_ngbr_vrtx_orig_order(nvc,ncells)')
    call ARRAYDESTROY (mesh_ngbr_cell_orig_order, 'Array mesh_ngbr_cell_orig_order(nfc,ncells)')

    return

  END SUBROUTINE EDIT_LONG_INIT

  SUBROUTINE EDIT_SHORT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Perform a "short edit" by computing and printing global diagnostics.
    !
    !=======================================================================
    use kinds, only: r8
    use constants_module,       only: one_half, zero
    use cutoffs_module,         only: alittle
    use fluid_data_module,      only: fluid_flow, qin, qout, isImmobile
    use fluid_type_module,      only: Div_c
    use kind_module,            only: int_kind, real_kind
    use matl_module,            only: GATHER_VOF
    use mesh_module,            only: Cell
    use nonlinear_solution,     only: NKuser, nonlinear_solutions, DEFAULT_NK_CONTROLS
    use truchas_env,            only: output_file_name
    use parameter_module,       only: ncells, ndim, nmat
    use pgslib_module,          only: PGSLIB_GLOBAL_MAXLOC, PGSLIB_GLOBAL_MAXVAL, &
                                      PGSLIB_GLOBAL_MINLOC, PGSLIB_GLOBAL_MINVAL, &
                                      PGSLIB_GLOBAL_SUM
    use projection_data_module, only: mac_projection_iterations,     &
                                      mac_projection_precond_iter,   &
                                      projection_precond_iterations, &
                                      projection_iterations
    use property_module,        only: Get_User_Material_ID
    use property_data_module,   only: Material_Name
    use property_module,        only: ENTHALPY_DENSITY_MATERIAL, DENSITY_MATERIAL
    use time_step_module,       only: cycle_number, t
    use viscous_data_module,    only: viscous_iterations
    use zone_module,            only: Zone
    use output_utilities,       only: ANNOUNCE_FILE_WRITE
    use solid_mechanics_module, only: STRESS_STRAIN_INVARIANTS
    use solid_mechanics_data,   only: Cell_Mech_Invariant, solid_mechanics
    use gap_output,         only: SET_GAP_ELEMENT_OUTPUT

    implicit none

    ! Argument List

    ! Local Variables
    character(LEN = 128)                        :: string, string2
    integer(KIND = int_kind)                    :: i, m, n, &
                                                   variables = 2*ndim + 6, &
                                                   nmechvar = 4
    integer(KIND = int_kind), dimension(1)      :: MaxLoc_L, MinLoc_L
    real(KIND = real_kind),   dimension(ncells) :: Enthalpy, KE, Mass, &
                                                   Matl_Mass, Tmp, Matl_Vol
    real(KIND = real_kind) :: Temperature
    type(CELL_MECH_INVARIANT), pointer, dimension(:) :: mech_info => NULL()

    real(r8), dimension(nmat) :: Material_Enthalpy, Material_KE, Material_Volume, Material_Mass
    real(r8) :: Material_Momentum(ndim,nmat)
    real(r8) :: Total_Enthalpy, total_KE, total_mass, Total_Momentum(ndim), total_volume

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Write out the header, time, and cycle
    write (string, 10) t, cycle_number
10  format (1x,14('* '),'SHORT EDIT: t = ',1pe13.5,', step = ',i6,14(' *'))
    call TLS_info ('')
    call TLS_info (string)
    call TLS_info ('')

    ! Print out a summary of the linear solver performance.
    call TLS_info ('')
    call TLS_info ('  Linear Solver Performance Summary')
    call TLS_info ('  ---------------------------------')
    write(string,'(4x,i4,1x,a)') mac_projection_iterations, 'MAC projection iterations'
    call TLS_info (string)
    write(string,'(4x,i4,1x,a)') viscous_iterations, 'viscous iterations'
    call TLS_info (string)
    write(string,'(4x,i4,1x,a)') projection_iterations, 'projection iterations'
    call TLS_info (string)
    write(string,'(4x,i4,1x,a)') mac_projection_precond_iter, 'MAC projection preconditioner iterations'
    call TLS_info (string)
    write(string,'(4x,i4,1x,a)') projection_precond_iterations, 'projection preconditioner iterations'
    call TLS_info (string)

    ! Print out a summary of the nonlinear solver performance.
    ! Loop over the NKuser derived type, printing out useful
    ! information if it's present.
    NK_DIAGNOSTICS: do i = 1,DEFAULT_NK_CONTROLS + nonlinear_solutions
    
       ! Cycle if this type contains nothing useful.
       if (NKuser(i)%Newton_tot == 0) cycle NK_DIAGNOSTICS

       ! Print out a summary of this nonlinear solver performance.
       call TLS_info ('')
       call TLS_info (repeat(' ',28) // 'Nonlinear Solver Performance Summary')
       call TLS_info ('')
       call TLS_info ('         Name  Convergence     Eps      Low Damper  High Damper  Iters  Linear Iters')
       call TLS_info ('         ----  -----------     ---      ----------  -----------  -----  ------------')
       write(string,'(9x,a,2x,1pe11.4,3(1x,1pe11.4),4x,i2,8x,i3)') TRIM(NKuser(i)%name), NKuser(i)%tolnewt, & 
           NKuser(i)%eps_NK, NKuser(i)%limit_low, NKuser(i)%limit_high, NKuser(i)%Newton_tot, NKuser(i)%linear_tot
       call TLS_info (string)

       call TLS_info ('')
       call TLS_info ('                    Iter  L2(Residual)  Linf(Residual)  Linf Location')
       call TLS_info ('                    ----  ------------  --------------  -------------')

       do n = 0,NKuser(i)%Newton_tot

          write (string, 18) n, NKuser(i)%L2(n), NKuser(i)%LI(n), &
                                    NKuser(i)%LI_Location(n)
18        format (21x,i2,3x,1pe11.4,4x,1pe11.4,6x,i6)
          call TLS_info (string)

       end do   
    end do NK_DIAGNOSTICS

    ! Compute mass, momentum, kinetic energy, enthalpy, and
    ! total energy for each material. Also compute the sum
    ! over all materials. Print out the results to out_tbrook.
    call TLS_info ('')
    call TLS_info (repeat(' ',32) // 'Mass-Momentum-Energy Summary')
    call TLS_info ('')
    call TLS_info (' Matl    Name     Volume      Mass               (X,Y,Z) Momentum               KE        Enthalpy')
    call TLS_info (' ----    ----     ------      ----               ----------------               --        --------')

    ! Preinitialize
    Mass = zero; Enthalpy = zero; KE = zero
    total_mass = zero; total_momentum = zero
    total_KE = zero; total_enthalpy = zero
 
    ! Zone kinetic energy density
    do n = 1,ndim
       KE = KE + Zone%Vc(n)**2
    end do
    KE = one_half*KE

    ! Loop over each material
    MATERIAL_SUMS: do m = 1,nmat

       ! Gather vof for this material
       call GATHER_VOF (m, Tmp)

       ! Sum material volume.
       Matl_Vol = Cell%Volume * Tmp
       Material_Volume(m) = PGSLIB_GLOBAL_SUM(Matl_Vol)
       if (ABS(Material_Volume(m)) <= alittle) Material_Volume(m) = zero

       ! Sum material mass.
       Matl_Mass = Matl_Vol*DENSITY_MATERIAL(m)
       Material_Mass(m) = PGSLIB_GLOBAL_SUM(Matl_Mass)
       if (ABS(Material_Mass(m)) <= alittle) Material_Mass(m) = zero

       ! Accumulate the total mass.
       Mass = Mass + Matl_Mass
       total_mass = total_mass + Material_Mass(m)

       ! Compute material momentum.
       MOMENTUM: do n = 1,ndim
          Material_Momentum(n,m) = PGSLIB_GLOBAL_SUM(Matl_Mass*Zone%Vc(n))
          if (ABS(Material_Momentum(n,m)) <= alittle .or. isImmobile(m)) Material_Momentum(n,m) = zero
          Total_Momentum(n) = Total_Momentum(n) + Material_Momentum(n,m)
       end do MOMENTUM

       ! Compute material kinetic energy.
       Material_KE(m) = PGSLIB_GLOBAL_SUM(Matl_Mass*KE)
       if (ABS(Material_KE(m)) <= alittle .or. isImmobile(m)) Material_KE(m) = zero
       total_KE = total_KE + Material_KE(m)

       ! Get the material enthalpy.
       ENTHALPY_LOOP: do n=1,ncells
         Tmp(n) = Matl_Vol(n) * ENTHALPY_DENSITY_MATERIAL(m,Zone(n)%Temp)
       end do ENTHALPY_LOOP

       ! Accumulate the material enthalpy.
       Material_Enthalpy(m) = PGSLIB_GLOBAL_SUM(Tmp)
       if (ABS(Material_Enthalpy(m)) <= alittle) Material_Enthalpy(m) = zero

       ! Accumulate the total enthalpy.
       Enthalpy = Enthalpy + Tmp
       total_enthalpy = total_enthalpy + Material_Enthalpy(m)

       ! Write out the summaries for this material.
       write (string, 25) Get_User_Material_ID(m), TRIM(material_name(m)), Material_Volume(m),      &
                                 Material_Mass(m), (Material_Momentum(n,m),n=1,ndim), &
                                 Material_KE(m), Material_Enthalpy(m)
25     format (2x,i2,2x,a8,1x,1pe11.4,1x,1pe11.4,1x,'(',1pe11.4,',',1pe11.4,',', &
               1pe11.4,')',1pe11.4,1x,1pe11.4)
       call TLS_info (string)

    end do MATERIAL_SUMS

    total_volume = PGSLib_Global_SUM(Cell%Volume)

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
       MinLoc_L = PGSLIB_GLOBAL_MINLOC(Tmp)
       MaxLoc_L = PGSLIB_GLOBAL_MAXLOC(Tmp)

       ! Write out the extrema location and values
       write (string2, 40) string, &
                                 PGSLIB_GLOBAL_MINVAL(Tmp), MinLoc_L, &
                                 PGSLIB_GLOBAL_MAXVAL(Tmp), MaxLoc_L
40     format (23x,a11,2x,1pe11.4,1x,i8,2x,1pe11.4,1x,i8)
       call TLS_info (string2)

    end do VARIABLE_EXTREMA

    ! Write out solid mechanics extrema if appropriate
    SOLID_MECHANICS_OUTPUT: if(solid_mechanics) then
       ! Put reasonable values in gap_elements
       call SET_GAP_ELEMENT_OUTPUT()
       ! Calculate data
       mech_info => STRESS_STRAIN_INVARIANTS()
       do i = 1,nmechvar
          ! Select the appropriate variable
          select case (i)
             case (1)
                ! Effective stress
                Tmp = mech_info%mises_stress
                string = 'Eff_Stress'
             case (2)
                ! Effective plastic strain
                Tmp = mech_info%eff_plastic_strain
                string = 'Plast_Strn'
             case (3)
                ! Hydrostatic stress
                Tmp = mech_info%mean_stress
                string = 'Mean_Stress'
             case (4)
                ! Volumetric strain
                Tmp = mech_info%volumetric_strain
                string = 'Vol_Strain'
          end select
          
          ! Find the extrema location
          MinLoc_L = PGSLIB_GLOBAL_MINLOC(Tmp)
          MaxLoc_L = PGSLIB_GLOBAL_MAXLOC(Tmp)
          
          ! Write out the extrema location and values
          write (string2, 40) string, &
                                 PGSLIB_GLOBAL_MINVAL(Tmp), MinLoc_L, &
                                 PGSLIB_GLOBAL_MAXVAL(Tmp), MaxLoc_L
          call TLS_info (string2)

       end do
       DEALLOCATE(mech_info)

    end if SOLID_MECHANICS_OUTPUT


    ! If fluid flow is on, write out velocity divergence norms
    ! and Inflow/Outflow volume if appropriate
    if (fluid_flow) then
       call TLS_info ('')
       call TLS_info ('                            Normalized Velocity Divergence Norms')
       call TLS_info ('')
       call TLS_info ('              Quantity            L1          L2         Linf     Linf Location')
       call TLS_info ('          ----------------        --          --         ----     -------------')
       write(string,45) Div_c%V_f%L1, Div_c%V_f%L2, Div_c%V_f%Linf, Div_c%V_f%Linf_Location
45     format(10x,'(Div)c*Vf*dt/Vol',3x,3(1pe11.4,1x),3x,i6)
       call TLS_info (string)
       call TLS_info ('')

       if (qin /= zero .or. qout /= zero) then
          call TLS_info ('')
          call TLS_info ('                                   Inflow-Outflow Summary')
          call TLS_info ('')
          call TLS_info ('                               Inflow Volume  Outflow Volume')
          call TLS_info ('                               -------------  --------------')
          write(string,50) qin, qout
50        format(32x,1pe11.4,4x,1pe11.4)
          call TLS_info (string)
       end if

    end if

    ! Write out the trailing header.
    write (string, 10) t, cycle_number
    call TLS_info ('')
    call TLS_info (string)
    call TLS_info ('')

    call ANNOUNCE_FILE_WRITE ('Short edit', output_file_name('log'))

    return

  END SUBROUTINE EDIT_SHORT

END MODULE EDIT_MODULE

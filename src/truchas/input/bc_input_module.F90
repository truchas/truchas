!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BC_INPUT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures for the input of boundary condition (BC) parameters.
  !
  !   Public Interface:
  !
  !     * call BC_INPUT ()
  !
  !       Default, read, check, and broadcast variables associated
  !       with the BC namelist.
  !
  ! Contains: BC_CHECK
  !           BC_DEFAULT
  !           BC_INPUT
  !           BC_INPUT_PARALLEL
  !
  ! Author(s): Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use truchas_logging_services
  implicit none
  private

  public :: BC_INPUT

CONTAINS

  SUBROUTINE BC_CHECK (fatal)
    !=======================================================================
    ! Purpose(s):
    !
    !   Check for obvious errors in the BC namelist input variables.
    !
    !=======================================================================
    use input_utilities, only: NULL_I, NULL_C
    use bc_data_module,   only: BC_Type, BC_Variable,                       &
                                Inflow_Material, Inflow_Index,              &
                                BC_Name,                                    &
                                nbc_surfaces, Type_Forms, Variable_Forms,   &
                                Conic_Tolerance,          &
                                Bounding_Box, Conic_Relation, Surface_Name, &
                                BC_Surface_Forms,                           &
                                Surfaces_In_This_BC, Srfmatl_Index,         &
                                Node_Disp_Coords, Mesh_Surface
    use parameter_module, only: bc_forms, maxmat, nvar, mbc_surfaces, mbcsrf
    use legacy_mesh_api, only: ndim
    use utilities_module, only: STRING_COMPARE
    use solid_mechanics_input, only: solid_mechanics
    use material_model_driver, only: matl_model
    
    ! Argument List
    logical, intent(INOUT) :: fatal

    ! Local Variables
    logical :: strings_match
    integer :: l, n, p, v
    character(128) :: errmsg

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Valid surface description.
    CHECK_BC_NAMES: do p = 1, nbc_surfaces-1
       if (BC_Name(p) /= 'Unnamed') then
          do l = p+1, nbc_surfaces
             if(BC_Name(l) == BC_Name(p)) then 
                ! Error found; punt.
                write (errmsg, 637) p, l, BC_Name(p)
637             format ('BCs ', i3, ' and ', i3, ' have the same name: ', a)
                call TLS_fatal (errmsg)
             end if
          end do
       end if
    end do CHECK_BC_NAMES

    CHECK_BC_SURFACES: do p = 1,nbc_surfaces


       ! Check surface name form.
       l = 0
       COUNT_BC_SURFACES: do n = 1,mbcsrf
          if (Surface_Name(n,p) /= 'none') then
             SURFACE_CHECK: do v = 1,mbc_surfaces
                call STRING_COMPARE (Surface_Name(n,p),BC_Surface_Forms(v),strings_match)
                if (strings_match) then

                   ! Now reset the surface name
                   select case(v)
                   case (1:2)
                      Surface_Name(n,p)  = 'conic'
                   case (15:16)
                      Surface_Name(n,p)  = 'from mesh file'
                      l = l + 1
                      Srfmatl_Index(n,p) =  l
                   case (17)
                      Surface_Name(n,p)  = 'node set'
                      l = l + 1
                      Srfmatl_Index(n,p) =  l
                   end select
                   exit SURFACE_CHECK
                end if

             end do SURFACE_CHECK
             if (v == mbc_surfaces + 1) then
                write (errmsg,1) p, TRIM(bc_name(p)),TRIM(Surface_Name(n,p))
1               format('BC surface ',i2,' named ',a,' has an unknown name: ','-',a)
                call TLS_error (errmsg)
                fatal = .true.
                cycle CHECK_BC_SURFACES
             end if
             Surfaces_In_This_BC(p) = Surfaces_In_This_BC(p) + 1 
          end if
       end do COUNT_BC_SURFACES
       if (Surfaces_In_This_BC(p) == 0) then
          write (errmsg,2) p,TRIM(bc_name(p))
2         format('Surfaces for BC namelist #',i2,' named ',a,' have not been specified!')
          call TLS_error (errmsg)
          fatal = .true.
          cycle CHECK_BC_SURFACES
       end if

       ! Check for node list if node set is specified
       do n = 1,Surfaces_In_This_BC(p)
          if (Surface_Name(n,p) == 'node set') then
             if (Node_Disp_Coords(1,1,p) == -1.0d10) then
                write (errmsg,5) p,TRIM(bc_name(p))
5               format('Node coordinates must be specified for BC namelist #',i2,' named ',a)
                call TLS_error (errmsg)
                fatal = .true.
                cycle CHECK_BC_SURFACES
             end if
          end if
       end do

       ! Check BC surface materials array and surfaces read from mesh.
       do n = 1,Surfaces_In_This_BC(p)
          v = Srfmatl_Index(n,p)
          if (v > 0) then
             select case(Surface_Name(n,p))
             case('from mesh file')
                if (Mesh_Surface(v,p) <= 0) then
                   write (errmsg,6) p,TRIM(bc_name(p))
6                  format('Mesh surface number for BC namelist #',i2,' named ',a,' must be positive!')
                   call TLS_error (errmsg)
                   fatal = .true.
                   cycle CHECK_BC_SURFACES
                end if
             end select
          end if
       end do

       ! Check Variable Name Form
       VARIABLE_CHECK: do v = 1,nvar
          do n = 1,bc_forms
             if (BC_Variable(p) /= '') then
                call STRING_COMPARE (BC_Variable(p),Variable_Forms(v,n),strings_match)
                if (strings_match) exit
             end if
          end do
          if (n /= bc_forms + 1) exit VARIABLE_CHECK
       end do VARIABLE_CHECK
       if (v == nvar + 1) then
          write (errmsg,10) p,TRIM(bc_name(p)),TRIM(BC_Variable(p))
10        format('surface ',i2,' named ',a,' has an unknown BC variable: ','-',a)
          call TLS_error (errmsg)
          fatal = .true.
          cycle CHECK_BC_SURFACES
       end if

       ! Check BC Type Form
       do n = 1, bc_forms
          if (BC_Type(p) /= '')then
             call STRING_COMPARE (BC_Type(p),Type_Forms(v,n),strings_match)
             if (strings_match) exit
          end if
       end do
       if (n == bc_forms + 1) then
          write (errmsg,15) p, TRIM(bc_name(p)), TRIM(BC_Type(p))
15        format('surface ',i2,' named ',a,' has an unknown BC type: ','-',a)
          call TLS_error (errmsg)
          fatal = .true.
          cycle CHECK_BC_SURFACES
       end if

       ! Overwrite Input Character Strings
       select case (v)
       case (5) ! Displacement
          if(.not.solid_mechanics) then
             write (errmsg,185) p
185          format('displacement boundary conditions specified', &
                  ' without a solid_mechanics solution. BC namelist #',i2)
             call TLS_warn (errmsg)
             cycle CHECK_BC_SURFACES
          endif
          BC_Variable(p) = 'displacement'
          select case (n)
          case (1)
             BC_Type(p) = 'x-traction'
          case (2)
             BC_Type(p) = 'y-traction'
          case (3)
             BC_Type(p) = 'z-traction'
          case (4)
             BC_Type(p) = 'x-displacement'
          case (5)
             BC_Type(p) = 'y-displacement'
          case (6)
             BC_Type(p) = 'z-displacement'
          case (7)
             BC_Type(p) = 'normal-displacement'
          case (8)
             BC_Type(p) = 'normal-traction'
          case (9)
             BC_Type(p) = 'free-interface'
          case (10)
             BC_Type(p) = 'normal-constraint'
          case (11)
             BC_Type(p) = 'contact'
          end select

          if (n > 11) then
             write (errmsg,187) p, BC_Type(p)
187          format('BC type in BC namelist ',i2,' is not valid',a)
             call TLS_error (errmsg)
             fatal = .true.
             cycle CHECK_BC_SURFACES
          end if

          
       end select
       
       ! Bounding box check; make sure that max >= min
       BOX_CHECK: do n = 1,ndim
          if (Bounding_Box(2,n,p) < Bounding_Box(1,n,p)) then
             write (errmsg,20) n, p
20           format('Bounding box coordinates in dimension ',i1,&
                ' for surface ',i3,' are not monotonically increasing')
             call TLS_error (errmsg)
             fatal = .true.
             cycle CHECK_BC_SURFACES
          end if
       end do BOX_CHECK
       
       ! Check to insure that surface tolerance is positive.
       do n = 1,Surfaces_In_This_BC(p)
          if (Conic_Tolerance(n,p) <= 0) then
             write (errmsg, 40) Conic_Tolerance(n,p)
40           format('Surface tolerance of ',1pe13.5,' cannot be negative!')
             call TLS_error (errmsg)
             fatal = .true.
          end if
       end do
    
       ! Check to insure that the conic relation is valid.
       RELATION_CHECK: do n = 1,Surfaces_In_This_BC(p)
          call STRING_COMPARE (Conic_Relation(n,p),'=',strings_match)
          if (strings_match) cycle RELATION_CHECK
          call STRING_COMPARE (Conic_Relation(n,p),'>',strings_match)
          if (strings_match) cycle RELATION_CHECK
          call STRING_COMPARE (Conic_Relation(n,p),'<',strings_match)
          if (strings_match) cycle RELATION_CHECK
          if (Surface_Name(n,p) == 'conic') then
             write (errmsg, 50) TRIM(Conic_Relation(n,p))
50           format('Conic relation of ',a,' is invalid!')
             call TLS_error (errmsg)
             fatal = .true.
          end if
       end do RELATION_CHECK
       
    end do CHECK_BC_SURFACES
    
  END SUBROUTINE BC_CHECK
  
  SUBROUTINE BC_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default variables in and related to the BC Namelist.
    !
    !=======================================================================
    use bc_data_module,   only: BC_Type, BC_Value, BC_Table, BC_Name,        &
                                BC_Variable, Inflow_Material,                &
                                Inflow_Temperature,                          &
                                Inflow_Index,                                &
                                Variable_Forms,                              &
                                nbc_surfaces, Conic_XX, Conic_YY,            &
                                Conic_ZZ, Conic_XY, Conic_XZ, Conic_YZ,      &
                                Conic_X, Conic_Y, Conic_Z, Conic_Constant,   &
                                Conic_Tolerance, Bounding_Box, Type_Forms,   &
                                Conic_Relation, Surface_Name,                &
                                BC_Surface_Forms,                            &
                                Surfaces_In_This_BC, Srfmatl_Index,          &
                                Node_Disp_Coords
    use input_utilities,  only: NULL_I, NULL_R, NULL_C

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Default BC Variables
    BC_name     = 'Unnamed'
    BC_Type     = 'none'
    BC_Value    = 0
    BC_Table    = NULL_R
    BC_Variable = 'none'

    ! Default Velocity Variables
    Inflow_Material      = NULL_C
    Inflow_Temperature   = NULL_R

    ! Default BC surface constants.
    Surface_Name        = 'none'
    Conic_XX            = 0
    Conic_YY            = 0
    Conic_ZZ            = 0
    Conic_XY            = 0
    Conic_XZ            = 0
    Conic_YZ            = 0
    Conic_X             = 0
    Conic_Y             = 0
    Conic_Z             = 0
    Conic_Constant      = 0
    Conic_Tolerance     = 1.0d-6
    Conic_Relation      = '='
    Bounding_Box(1,:,:) = -huge(1.0d0)
    Bounding_Box(2,:,:) = +huge(1.0d0)
    Node_Disp_Coords(:,:,:)  = -1.0d10

    ! Matl surface number.
    Srfmatl_Index = 0

    ! Inflow surface number.
    Inflow_Index = 0

    ! Number of BC surfaces
    nbc_surfaces        = 0
    Surfaces_In_This_BC = 0

    ! BC surface character strings
    BC_Surface_Forms    = ''

    BC_Surface_Forms(1)  = 'conic'
    BC_Surface_Forms(2)  = 'conic function'
    BC_Surface_Forms(15) = 'from mesh file'
    BC_Surface_Forms(16) = 'in mesh file'
    BC_Surface_Forms(17) = 'node set'

    ! BC Types Character Strings
    Type_Forms = ''

    Type_Forms(5,1) = 'x-traction'     ! 5: Displacement
    Type_Forms(5,2) = 'y-traction'
    Type_Forms(5,3) = 'z-traction'
    Type_Forms(5,4) = 'x-displacement'
    Type_Forms(5,5) = 'y-displacement'
    Type_Forms(5,6) = 'z-displacement'
    Type_Forms(5,7) = 'normal-displacement'
    Type_Forms(5,8) = 'normal-traction'
    Type_Forms(5,9) = 'free-interface'
    Type_Forms(5,10) = 'normal-constraint'
    Type_Forms(5,11) = 'contact'


    ! BC Variable Character Strings
    Variable_Forms = ''

    Variable_Forms(5,1) = 'displacement'  ! 5: Displacement
    Variable_Forms(5,2) = 'disp'

  END SUBROUTINE BC_DEFAULT

  SUBROUTINE BC_INPUT (lun)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read BC namelist. The variables set in this routine are listed
    !   in the BC namelist. If there is no BC namelist in the input deck,
    !   then defaults are used. If there is an error in the BC namelist,
    !   then execution terminates.
    !
    !=======================================================================
    use kinds, only: r8
    use bc_data_module,         only: BC_Table_Array => BC_Table
    use bc_data_module,         only: BC_Type, BC_Value, BC_Name,                  &
                                      BC_Variable, Inflow_Material,                &
                                      Inflow_Temperature,                          &
                                      Conic_XX, Conic_YY,                          &
                                      Conic_ZZ, Conic_XY, Conic_XZ, Conic_YZ,      &
                                      Conic_X, Conic_Y, Conic_Z,                   &
                                      Conic_Constant, Conic_Tolerance,             &
                                      Bounding_Box, nbc_surfaces, Surface_Name,    &
                                      Conic_Relation, Node_Disp_Coords, Mesh_Surface
    use input_utilities,        only: NULL_R
    use input_utilities,        only: seek_to_namelist, NULL_I, NULL_C
    use parallel_info_module,   only: P_Info
    use parameter_module,       only: string_dim, mbc_surfaces
    use string_utilities, only: i_to_c

    integer, intent(in) :: lun

    ! Local Variables
    logical :: fatal, found
    integer :: ios, bcs
    character(128) :: iom
    character(128) :: errmsg
    character(:), allocatable :: label
    real(r8) :: BC_Table(size(BC_Table_Array,1),size(BC_Table_Array,2))

    ! Define BC Namelist
    namelist /BC/ BC_Type, BC_Value, BC_Table, BC_Variable, BC_Name,             &
                  Inflow_Material, Inflow_Temperature,                           &
                   Conic_XX, Conic_YY, Conic_ZZ, Conic_XY,                       &
                  Conic_XZ, Conic_YZ, Conic_X, Conic_Y, Conic_Z, Conic_Constant, &
                  Conic_Tolerance, Bounding_Box, Surface_Name,                   &
                  Node_Disp_Coords, Mesh_Surface

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Error Fatal Flag
    fatal = .false.

    ! Read Notice
    call TLS_info ('Reading BC namelists ...')

    ! Rewind the input deck unit number.
    if (P_Info%IOP) rewind lun

    ! Set BC defaults.
    call BC_DEFAULT ()

    ! Loop over max possible BC surfaces
    BC_NAMELIST_LOOP: do bcs = 1,mbc_surfaces+1

       ! Default the zeroth element of all BC input variable arrays.
       BC_Name(0)               = 'Unnamed'
       BC_Type(0)               = 'none'
       BC_Value(:,0)            = 0
       BC_Table                 = NULL_R
       BC_Variable(0)           = 'none'
       Inflow_Material(0)       = NULL_C
       Inflow_Temperature(0)    = NULL_R
       Surface_Name(:,0)        = 'none'
       !Conic_Relation(:,0)      = 'none'
       Conic_XX(:,0)            = 0
       Conic_YY(:,0)            = 0
       Conic_ZZ(:,0)            = 0
       Conic_XY(:,0)            = 0
       Conic_XZ(:,0)            = 0
       Conic_YZ(:,0)            = 0
       Conic_X(:,0)             = 0
       Conic_Y(:,0)             = 0
       Conic_Z(:,0)             = 0
       Conic_Constant(:,0)      = 0
       Conic_Tolerance(:,0)     = 1.0d-6
       Bounding_Box(1,:,0)      = -huge(1.0d0)
       Bounding_Box(2,:,0)      =  huge(1.0d0)
       Node_Disp_Coords(:,:,0)  = -1.0d10
       Mesh_Surface(:,0)        = 0

       ! Search for and read this BC namelist only on IO PE.
       IO_PE_ONLY: if (P_Info%IOP) then

          ! Find the next BC namelist.
          call seek_to_namelist (lun, 'BC', found)
          if (.not.found) exit

          nbc_surfaces = nbc_surfaces + 1
          label = 'BC[' // i_to_c(nbc_surfaces) // ']'
          
          ! Check that we aren't about to overflow the BC arrays.
          if (nbc_surfaces > mbc_surfaces) then
             write (errmsg, 17) mbc_surfaces
17           format ('too many BC namelists found; limited to MBC_SURFACES=', i0)
             call TLS_error (errmsg)
             fatal = .true.
             exit BC_NAMELIST_LOOP
          end if
          
          ! Read the Namelist.
          read(lun,nml=BC,iostat=ios,iomsg=iom)

          if (ios /= 0) then
             call TLS_error ('error reading ' // label // ' namelist: ' // trim(iom))
             fatal= .true.
             exit BC_NAMELIST_LOOP
          end if   

          ! Copy zeroth element to actual position.
          BC_Name(nbc_surfaces)              = BC_Name(0)
          BC_Type(nbc_surfaces)              = BC_Type(0)
          BC_Value(:,nbc_surfaces)           = BC_Value(:,0)
          BC_Table_Array(:,:,nbc_surfaces)   = BC_Table
          BC_Variable(nbc_surfaces)          = BC_Variable(0)
          Inflow_Material(nbc_surfaces)      = Inflow_Material(0)
          Inflow_Temperature(nbc_surfaces)   = Inflow_Temperature(0)
          Surface_Name(:,nbc_surfaces)       = Surface_Name(:,0)

          Conic_Relation(:,nbc_surfaces)     = '=' !Conic_Relation(:,0)
          Conic_XX(:,nbc_surfaces)           = Conic_XX(:,0)
          Conic_YY(:,nbc_surfaces)           = Conic_YY(:,0)
          Conic_ZZ(:,nbc_surfaces)           = Conic_ZZ(:,0)
          Conic_XY(:,nbc_surfaces)           = Conic_XY(:,0)
          Conic_XZ(:,nbc_surfaces)           = Conic_XZ(:,0)
          Conic_YZ(:,nbc_surfaces)           = Conic_YZ(:,0)
          Conic_X(:,nbc_surfaces)            = Conic_X(:,0)
          Conic_Y(:,nbc_surfaces)            = Conic_Y(:,0)
          Conic_Z(:,nbc_surfaces)            = Conic_Z(:,0)
          Conic_Constant(:,nbc_surfaces)     = Conic_Constant(:,0)
          Conic_Tolerance(:,nbc_surfaces)    = Conic_Tolerance(:,0)
          Bounding_Box(:,:,nbc_surfaces)     = Bounding_Box(:,:,0)
          Node_Disp_Coords(:,:,nbc_surfaces) = Node_Disp_Coords(:,:,0)  
          Mesh_Surface(:,nbc_surfaces)       = Mesh_Surface(:,0)
       end if IO_PE_ONLY

    end do BC_NAMELIST_LOOP

    ! Error Check
    call TLS_fatal_if_any (fatal, 'terminating execution due to previous input errors')

    ! Broadcast namelist data.
    call BC_INPUT_PARALLEL ()

    ! Check Namelist
    call BC_CHECK (fatal)

    ! Error Check
    call TLS_fatal_if_any (fatal, 'terminating execution due to previous input errors')

    select case (nbc_surfaces)
    case (0)
      call TLS_info('  none found; using defaults')
    case (1)
      call TLS_info('  read 1 BC namelist')
    case default
      call TLS_info('  read ' // i_to_c(nbc_surfaces) // ' BC namelists')
    end select

  END SUBROUTINE BC_INPUT

  SUBROUTINE BC_INPUT_PARALLEL ()
    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast variables in the BC namelist to all processors.
    !
    !======================================================================
    use bc_data_module,       only: BC_Type, BC_Value, BC_Table, BC_Name,     &
                                    BC_Variable, Inflow_Material,             &
                                    Inflow_Temperature,                       &
                                    Conic_XX, Conic_YY,                       &
                                    Conic_ZZ, Conic_XY, Conic_XZ, Conic_YZ,   &
                                    Conic_X, Conic_Y, Conic_Z,                &
                                    Conic_Constant, Conic_Tolerance,          &
                                    Bounding_Box, nbc_surfaces, Surface_Name, &
                                    Conic_Relation,                           &
                                    Surfaces_In_This_BC, Node_Disp_Coords,    &
                                    Mesh_Surface
    use parallel_info_module, only: P_Info
    use pgslib_module,        only: PGSLib_BCAST

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Broadcast BC Namelist.
    if (.not. P_Info%UseGlobalServices) then

       call PGSLib_BCAST (BC_name)
       call PGSLib_BCAST (BC_Type)
       call PGSLib_BCAST (BC_Value)
       call PGSLib_BCAST (BC_Table)
       call PGSLib_BCAST (BC_Variable)

       call PGSLib_BCAST (Inflow_Material)
       call PGSLib_BCAST (Inflow_Temperature)

       call PGSLib_BCAST (Conic_XX)
       call PGSLib_BCAST (Conic_YY)
       call PGSLib_BCAST (Conic_ZZ)
       call PGSLib_BCAST (Conic_XY)
       call PGSLib_BCAST (Conic_XZ)
       call PGSLib_BCAST (Conic_YZ)
       call PGSLib_BCAST (Conic_X)
       call PGSLib_BCAST (Conic_Y)
       call PGSLib_BCAST (Conic_Z)
       call PGSLib_BCAST (Conic_Constant)
       call PGSLib_BCAST (Conic_Tolerance)
       call PGSLib_BCAST (Conic_Relation)
       call PGSLib_BCAST (Surface_Name)
       call PGSLib_BCAST (Bounding_Box)
       call PGSLib_BCAST (nbc_surfaces)
       call PGSLib_BCAST (Surfaces_In_This_BC)
       call PGSLib_BCAST (Node_Disp_Coords)
       call PGSLib_BCAST (Mesh_Surface)

    end if

  END SUBROUTINE BC_INPUT_PARALLEL

END MODULE BC_INPUT_MODULE

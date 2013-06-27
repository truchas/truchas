
!!$   Future direction for mesh output
!!$
!!$   Should split mesh into three sections:
!!$  MAKE MESHWRITERVERSION 1.0001
!!$
!!$  TYPE (hex or tet)
!!$
!!$  CELLS
!!$    VERTICES
!!$    VOLUMES
!!$    CENTROIDS
!!$    BOUNDARYFLAGS
!!$    NUMNEIGHBOR
!!$    PERMUTE
!!$    UNPERMUTE
!!$    BLOCKID
!!$
!!$  FACES
!!$    TOTALBOUNDARYFACES
!!$    NUMSIDESETS
!!$    SIDESETDATA
!!$
!!$  VERTICES
!!$    COORDINATES
!!$    RSUMBYRVOL
!!$    PERMUTE
!!$    UNPERMUTE
!!$
!!$  Each of these will have

#include "f90_assert.fpp"

module tbrook_utility
  use kinds, only: r8
  use tbrook_module, only: tbu_make_file_entry, B_Strlen
  use truchas_logging_services
  implicit none
  private

  public :: TBU_WriteBasicData
  public :: TBU_WriteDefaultMesh
  public :: TBU_MeshWriter
  public :: TBU_Zone_Out
  public :: TBU_flow_Out
  public :: TBU_Matl_Out
  public :: TBU_ENSIGHT_ADD_VARIABLE
  public :: TBU_TIMESTEP_OUTPUT
  public :: TBU_LASTSTEP_OUTPUT
  public :: TBU_GMV_WRITE_HEADER
  public :: TBU_GMV_WRITE_FOOTER
  public :: TBU_GMV_ADD_VARIABLE
  character(LEN=6) :: TBU_Version = '1.0000'

  integer, public, parameter :: CELL_MAP = 1
  integer, public, parameter :: FACE_MAP = 2
  integer, public, parameter :: EDGE_MAP = 3
  integer, public, parameter :: NODE_MAP = 4

contains

  SUBROUTINE TBU_ENSIGHT_ADD_VARIABLE(Name, filename, FORMAT, R1, R2, R3, I1, I2, I3, L1, L2, L3, iStatus)
    use tbrook_module
    use file_utility, only: MAKE_FILE_NAME

    !Arguments
    character(len=*), intent(in)    :: Name
    character(len=*), intent(in)    :: filename
    integer,          intent(inout) :: iStatus

    character(len=*), optional, intent(in) :: format
    real(r8), optional, dimension(:),     intent(in) :: R1
    real(r8), optional, dimension(:,:),   intent(in) :: R2
    real(r8), optional, dimension(:,:,:), intent(in) :: R3
    integer,  optional, dimension(:),     intent(in) :: I1
    integer,  optional, dimension(:,:),   intent(in) :: I2
    integer,  optional, dimension(:,:,:), intent(in) :: I3
    logical,  optional, dimension(:),     intent(in) :: L1
    logical,  optional, dimension(:,:),   intent(in) :: L2
    logical,  optional, dimension(:,:,:), intent(in) :: L3

    Type(Brook), Target :: B_E

    ! density
    iStatus=0
    call TBrook_Set(B_E, &
                    FILE=TRIM(filename), &
                    iStatus=iStatus)
    call TBrook_Write(B_E, Variable=TRIM(Name), iStatus=iStatus)
    call TBrook_Write(B_E, Variable=' (Truchas: "', iStatus=iStatus)
    call TBrook_Write(B_E, Variable=TRIM(MAKE_FILE_NAME('inp')), iStatus=iStatus)
    call TBrook_Write(B_E, Variable='")', iStatus=iStatus)
    call TBrook_Endline(B_E, iStatus=iStatus)

    call TBrook_Write(B_E, Variable='part 1', advance=.true., iStatus=iStatus)
    call TBrook_Write(B_E, Variable='hexa8', advance=.true., iStatus=iStatus)

    if (present(r1)) then
       call TBrook_write(B_E, variable=r1, FORMAT=FORMAT, Advance=.true., iStatus=iStatus)
    else if (present(r2)) then
       call TBrook_write(B_E, variable=r2, FORMAT=FORMAT, Advance=.true., iStatus=iStatus)
    else if (present(r3)) then
       call TBrook_write(B_E, variable=r3, FORMAT=FORMAT, Advance=.true., iStatus=iStatus)
    else if (present(i1)) then
       call TBrook_write(B_E, variable=i1, FORMAT=FORMAT, Advance=.true., iStatus=iStatus)
    else if (present(i2)) then
       call TBrook_write(B_E, variable=i2, FORMAT=FORMAT, Advance=.true., iStatus=iStatus)
    else if (present(i3)) then
       call TBrook_write(B_E, variable=i3, FORMAT=FORMAT, Advance=.true., iStatus=iStatus)
    else if (present(l1)) then
       call TBrook_write(B_E, variable=l1, FORMAT=FORMAT, Advance=.true., iStatus=iStatus)
    else if (present(l2)) then
       call TBrook_write(B_E, variable=l2, FORMAT=FORMAT, Advance=.true., iStatus=iStatus)
    else if (present(l3)) then
       call TBrook_write(B_E, variable=l3, FORMAT=FORMAT, Advance=.true., iStatus=iStatus)
    end if
    call TBrook_Close(B_E, iStatus=iStatus)
    call TBrook_Destroy(B_E, iStatus=iStatus)

  END SUBROUTINE TBU_ENSIGHT_ADD_VARIABLE

  SUBROUTINE TBU_GMV_WRITE_HEADER(B_gmv, filename, iStatus, faceBased, mask)
    !---------------------------------------------------------------------------
    ! Purpose:
    !   Write a GMV graphics dump
    !   B_gmv: the brook to which we will write
    !   filename: the name of the file
    !   istatus: return status, inout
    !   faceBased: logical, if true, write face based gmv file
    !   mask: for face based, only write faces that are true, dimensions nfc*ncells
    !---------------------------------------------------------------------------
    use mesh_module,          only: Mesh, Vertex, Face_Vrtx, unpermute_mesh_vector
    use parameter_module,     only: ncells, ndim, nvc, nfc, ncells_tot, nnodes_tot
    use pgslib_module,        only: pgslib_global_sum, &
                                    pgslib_collate,    &
                                    pgslib_bcast
    use tbrook_module,        only: Brook,           &
                                    TBrook_Set,      &
                                    TBrook_Write,    &
                                    TBrook_Endline,  &
                                    ASSIGNMENT(=)
    use parallel_info_module, only: p_info

    ! Argument List
    type(Brook), target :: B_gmv
    character(LEN=*) :: filename
    integer, intent(inout) :: iStatus
    logical, intent(in), optional :: faceBased
    logical, intent(in), optional, dimension(:) :: mask

    ! Local Variables
    integer :: n, f, v, m, out_stat, istat, nf
    logical :: doCells = .true.

    integer, dimension(nvc,ncells) :: Vrtx_Ngbr
    character(LEN=70), dimension(ncells)     :: Vrtx_Ngbr_String
    integer, dimension(:,:), allocatable :: faceInfo
    integer, dimension(nfc*ncells+nfc) :: faceMask
    integer, dimension(:), allocatable :: nfArray
    !---------------------------------------------------------------------------
    if ( present(faceBased)) then
       doCells = .not. faceBased
    else 
       doCells = .true.
    end if
    

    select case (ndim)
    case (3)
    case Default
       INSIST(.false.)
    end select

    ! replace blanks with underscores in name of LINEAR_SOLVER

    ! open the graphics output file
    out_stat=0
    call TBrook_Set(B=B_gmv,             &
                    FILE=TRIM(filename), &
                    FORM="ascii",        &
                    DFORMAT='(7(1es12.5,1x))', &
                    iStatus=out_stat)
    ! Write prologue.
    if(out_stat==0) CALL TBrook_Write(B_gmv, Variable='gmvinput ascii', ADVANCE=.true., iStatus=out_stat)

    ! Write node data.
    if(out_stat==0) CALL TBrook_Write(B_gmv, Variable='nodes ', iStatus=out_stat)
    if(out_stat==0) CALL TBrook_Write(B_gmv, Variable=nnodes_tot,   iStatus=out_stat)
    if(out_stat==0) call TBrook_EndLine(B_gmv, out_stat)
    do n = 1, ndim
       if(out_stat==0) CALL TBrook_Write(B_gmv, Variable=Vertex%Coord(n), iStatus=out_stat)
       if(out_stat==0) call TBrook_EndLine(B_gmv, out_stat)
    end do

    if ( doCells ) then
       ! Write element/vertex connectivity data.
       if(out_stat==0) CALL TBrook_Write(B_gmv, Variable='cells ', iStatus=out_stat)
       if(out_stat==0) CALL TBrook_Write(B_gmv, Variable=ncells_tot,   iStatus=out_stat)
       if(out_stat==0) call TBrook_EndLine(B_gmv, out_stat)
       do v = 1,nvc
          Vrtx_Ngbr(v,:) = Mesh%Ngbr_Vrtx_Orig(v)
       end do
       select case (ndim)
       case (2)
          write (Vrtx_Ngbr_String, 1) Vrtx_Ngbr
1         format (4(i7,1x))
          Vrtx_Ngbr_String = 'hex 4 ' // Vrtx_Ngbr_String
       case (3)
          write (Vrtx_Ngbr_String, 2) Vrtx_Ngbr
2         format (8(i7,1x))
          Vrtx_Ngbr_String = 'hex 8 ' // Vrtx_Ngbr_String
       end select
       if(out_stat==0) CALL TBrook_Write(B_gmv, Variable=Vrtx_Ngbr_String, Advance=.true., iStatus=out_stat)
    else 
       ! Face based data
       allocate(faceInfo(7,ncells*nfc), STAT=istat)
       call TLS_fatal_if_any (istat/=0, 'TBU_gmv_write_header: allocation error')

       if ( .not. present (mask) ) then
          ! write all faces
          faceInfo = 0
          if(out_stat==0) then
             m = 0
             do n = 1, ncells
                do f = 1, nfc 
                   m = m + 1
                   faceInfo(1,m) = 4
                   do v = 1, 4
                      faceInfo(v+1,m) = Mesh(n)%Ngbr_Vrtx_Orig(Face_Vrtx(f,v))
                   end do
                   faceInfo(6,m) = unpermute_mesh_vector(n)
                   faceInfo(7,m) = Mesh(n)%Ngbr_cell_orig(f)
                end do
             end do
          end if
          if(out_stat==0) CALL TBrook_Write(B_gmv, Variable='faces ',     advance=.false., iStatus=out_stat)
          if(out_stat==0) CALL TBrook_Write(B_gmv, Variable=pgslib_global_sum(m), advance=.false., iStatus=out_stat)
          if(out_stat==0) CALL TBrook_Write(B_gmv, Variable=ncells_tot,   iStatus=out_stat)
          if(out_stat==0) call TBrook_EndLine(B_gmv, out_stat)
          if(out_stat==0) CALL TBrook_Write(B_gmv, Variable=faceInfo(:,1:m),   &
               advance=.true., format='(7(i10,1x))', iStatus=out_stat)
       else
          if(out_stat==0) then
             where (mask) 
                faceMask = 1
             elsewhere
                faceMask = 0
             end where

             allocate(nfArray(p_info%nPE), stat = m)
             call TLS_fatal_if_any (m/=0, 'tbu_gmvWRiteHEader: nfArray allocation error')

             nf = sum(facemask)
             call pgslib_collate(nfArray,(/nf/))
             call pgslib_bcast(nfArray)

             nf = 0
             do n = 1, p_info%thisPE-1
                nf = nf + nfArray(n)
             end do
             faceInfo = 0
             m = 0
             do n = 1, ncells
                do f = 1, nfc 
                   if(facemask(n*6+f) == 1) then
                      m = m + 1
                      faceInfo(1,m) = 4
                      do v = 1, 4
                         faceInfo(v+1,m) = Mesh(n)%Ngbr_Vrtx_Orig(Face_Vrtx(f,v))
                      end do
                      faceInfo(6,m) = nf+m
                      ! done by initialization: faceInfo(7,m) = 0
                   end if
                end do
             end do
          end if
          if(out_stat==0) CALL TBrook_Write(B_gmv, Variable='faces ',     advance=.false., iStatus=out_stat)
          if(out_stat==0) CALL TBrook_Write(B_gmv, Variable=sum(nfArray), advance=.false., iStatus=out_stat)
          if(out_stat==0) CALL TBrook_Write(B_gmv, Variable=pgslib_global_sum(m),   iStatus=out_stat)
          if(out_stat==0) call TBrook_EndLine(B_gmv, out_stat)
          if(out_stat==0) CALL TBrook_Write(B_gmv, Variable=faceInfo(:,1:m),   &
               advance=.true., format='(7(i10,1x))', iStatus=out_stat)
          deallocate(faceInfo)
       end if

    end if


    ! Write variable data.
    if(out_stat==0) CALL TBrook_Write(B_gmv, Variable='variable ', Advance=.true., iStatus=out_stat)

    iStatus = out_stat

  END SUBROUTINE TBU_GMV_WRITE_HEADER

  SUBROUTINE TBU_GMV_WRITE_FOOTER(B_Gmv, t, WriteEndVars, iStatus)
    !---------------------------------------------------------------------------
    ! Purpose:
    !   Write a GMV graphics dump
    !---------------------------------------------------------------------------
    use output_utilities,     only: ANNOUNCE_FILE_WRITE
    use tbrook_module,        only: Brook,           &
                                    TBrook_Write,    &
                                    TBrook_Close,    &
                                    TBrook_Destroy,  &
                                    ASSIGNMENT(=)
    use tbrook_module,         only: tBrook_File

    ! Argument List
    type(Brook), target :: B_Gmv
    real(r8), intent(in) :: t
    integer, intent(inout) :: iStatus
    integer :: out_stat
    Logical, optional, intent(in) :: WriteEndVars

    ! Local Variable            
    logical :: endvars
    !---------------------------------------------------------------------------

    if ( present(writeendvars)) then
       endvars = writeendvars
    else 
       endvars=.true.
    end if
    out_stat = 0
    if(out_stat==0 .and. endvars) CALL TBrook_Write(B_gmv, Variable='endvars', Advance=.true., iStatus=out_stat)

    ! Write time and trailer.
    if(out_stat==0) CALL TBrook_Write(B_gmv, Variable='probtime ', Advance=.false., iStatus=out_stat)
    if(out_stat==0) CALL TBrook_Write(B_gmv, Variable=t, Format='(1es13.6)', Advance=.true., iStatus=out_stat)
    if(out_stat==0) CALL TBrook_Write(B_gmv, Variable='endgmv', Advance=.true., iStatus=out_stat)

    call TBrook_Close(B_gmv, iStatus=out_stat)
    call TBrook_Destroy(B_gmv, iStatus=out_stat)


    call ANNOUNCE_FILE_WRITE ('GMV graphics dump into file', TBrook_File(B_Gmv))

    iStatus = out_stat

  END SUBROUTINE TBU_GMV_WRITE_FOOTER

  SUBROUTINE TBU_GMV_ADD_VARIABLE (B_Gmv, Name, Format, RVar, IVar, LVar, iStatus)
    !---------------------------------------------------------------------------
    ! Purpose:
    !   Write a GMV graphics dump
    !---------------------------------------------------------------------------
    use tbrook_module, only: Brook, TBrook_Write, ASSIGNMENT(=)

    ! Argument List
    type(Brook), target    :: B_Gmv
    integer, intent(inout) :: iStatus
    character(LEN=*), intent(in) :: name
    real(r8), dimension(:), intent(in), optional :: RVar
    integer,  dimension(:), intent(in), optional :: IVar
    logical,  dimension(:), intent(in), optional :: LVar
    character(len=*),       intent(in), optional :: Format

    ! Local Variables
    integer :: out_stat
    integer, dimension(:), allocatable :: iLVar

    !---------------------------------------------------------------------------
    out_stat = 0
    if(out_stat==0) CALL TBrook_Write(B_gmv, Variable=name, Advance=.true., iStatus=out_stat)
    
    if (out_stat==0 .and. Present(RVar)) then
       CALL TBrook_Write(B_gmv, Variable=RVar, Format=Format, Advance=.true., iStatus=out_stat)
    else if (out_stat==0 .and. Present(IVar)) then
       CALL TBrook_Write(B_gmv, Variable=IVar, Format=Format, Advance=.true., iStatus=out_stat)
    else if (out_stat==0 .and. Present(LVar)) then
       allocate(iLVar(size(LVar)), STAT=out_stat)
       if (out_stat == 0 ) then
          iLVar = 0
          where(LVar)
             iLVar = 1
          end where
       end if
       CALL TBrook_Write(B_gmv, Variable=iLVar, Format=Format, Advance=.true., iStatus=out_stat)
       if(allocated(iLVar)) then
          deallocate(iLVar)
       end if
    else
       INSIST(.false.) 
    end if 
 
    iStatus = out_stat 
 
  END SUBROUTINE TBU_GMV_ADD_VARIABLE 
 
 
  SUBROUTINE TBU_ABORT(Message, SUBROUTINENAME, error_code, iStatus) 
    use tbrook_module, only: B_Stdout, tBrook_Write 

    character(LEN=*) :: Message 
    character(LEN=*) :: SubroutineName 
    integer          :: error_code 
    integer, optional:: iStatus 
 
    integer :: i 
    i = 0 
    if ( present(iStatus) ) then 
       iStatus = error_code 
       return 
    else 
       call TBrook_Write(B_stdout,       &
                        Variable='ERRORS: '// TRIM(Message) //' in '//TRIM(SUBROUTINENAME), &
                        iStatus=i)
    end if
  END SUBROUTINE TBU_ABORT

  SUBROUTINE TBU_WriteProgramSpecifications (B, iStatus) 
    !--------------------------------------------------------------------------- 
    ! Purpose: 
    ! 
    !   print build time and run time information 
    !--------------------------------------------------------------------------- 
    use tbrook_module,          only: BROOK,               & 
                                       TBrook_WriteXMLTag,  & 
                                       TBrook_OpenXMLTag,   & 
                                       TBrook_CloseXMLTag
 
    use code_module,             only: code_name,          & 
                                       code_version,       & 
                                       libraries,          & 
                                       build_date,         & 
                                       build_architecture, & 
                                       build_flags,        & 
                                       build_host,         & 
                                       run_architecture,   & 
                                       run_host 
    use parameter_module,        only: string_len 
    use parallel_info_module,    only: p_info 
    use utilities_module,        only: TIMESTAMP 
 
    type(Brook), target :: B 
    integer :: iStatus 
 
    ! Local Variables 
    character (LEN=string_len) :: string, istring 
    character (LEN=string_len) :: run_date 
 
    !--------------------------------------------------------------------------- 
    ! collect info we don't already have 
 
 
    run_architecture = "" 
    run_host         = "" 
    call getrunhostinfo (run_architecture, run_host) 
    call TIMESTAMP (run_date) 
 
    if ( iStatus /= 0 ) return 
    if ( .not. p_info%IOP) return 
 
    if ( iStatus == 0 ) call Tbrook_OpenXMLTag(B=B, XMLTag="PROGRAMSPECIFICATIONS",iStatus=iStatus) 
 
1   Format('TSTRLEN="',i10,'"') 
    ! write out code info 
    write(string,1) LEN(code_name) 
    if ( iStatus == 0 ) & 
         call Tbrook_WriteXMLTag(B=B, XMLTag="CODE",XMLAttributes=string,XMLStringData=TRIM(code_name), iStatus=iStatus) 
 
    write(string,1) LEN(code_version) 
    if ( iStatus == 0 ) & 
          call Tbrook_WriteXMLTag(B=B, XMLTag="VERSION",XMLAttributes=string,XMLStringData=TRIM(code_version), iStatus=iStatus) 
 
    write(string,1) LEN(libraries) 
    if ( iStatus == 0 ) & 
          call Tbrook_WriteXMLTag(B=B, XMLTag="LIBRARIES",XMLAttributes=string,XMLStringData=TRIM(libraries), iStatus=iStatus) 
 
    write(string,1) LEN(build_architecture) 
    if ( iStatus == 0 ) & 
          call Tbrook_WriteXMLTag(B=B, XMLTag="BUILDARCHITECTURE",XMLAttributes=string,& 
                                  XMLStringData=TRIM(build_architecture), iStatus=iStatus) 
 
    write(string,1) LEN(build_date) 
    if ( iStatus == 0 ) & 
          call Tbrook_WriteXMLTag(B=B, XMLTag="BUILDDATETIME",XMLAttributes=string,& 
                                  XMLStringData=TRIM(build_date), iStatus=iStatus) 
 
    write(string,1) LEN(build_flags) 
    if ( iStatus == 0 ) & 
          call Tbrook_WriteXMLTag(B=B, XMLTag="BUILDFLAGS",XMLAttributes=string,XMLStringData=TRIM(build_flags), iStatus=iStatus) 
 
    write(string,1) LEN(build_host) 
    if ( iStatus == 0 ) & 
          call Tbrook_WriteXMLTag(B=B, XMLTag="BUILDHOST",XMLAttributes=string,XMLStringData=TRIM(build_host), iStatus=iStatus) 
 
    write(string,1) LEN(run_architecture) 
    if ( iStatus == 0 ) & 
          call Tbrook_WriteXMLTag(B=B, XMLTag="RUNARCHITECTURE",XMLAttributes=string,& 
                                  XMLStringData=TRIM(run_architecture), iStatus=iStatus) 
 
    write(string,1) LEN(run_host) 
    if ( iStatus == 0 ) & 
          call Tbrook_WriteXMLTag(B=B, XMLTag="RUNHOST",XMLAttributes=string,XMLStringData=TRIM(run_host), iStatus=iStatus) 
 
    write(string,1) LEN(run_date) 
    if ( iStatus == 0 ) & 
          call Tbrook_WriteXMLTag(B=B, XMLTag="RUNDATE",XMLAttributes=string,XMLStringData=TRIM(run_date), iStatus=iStatus) 
 
    Write (istring, *) p_info%nPE 
    if ( iStatus == 0 ) & 
          call Tbrook_WriteXMLTag(B=B, XMLTag="RUNPROCESSORS", XMLStringData=TRIM(ADJUSTL(iString)), iStatus=iStatus) 
 
    if (p_info%nPE > 1) then 
       Write (istring, *) p_info%IO_ROOT_PE 
       if ( iStatus == 0 ) & 
          call Tbrook_WriteXMLTag(B=B, XMLTag="IOPROCESSOR",XMLStringData=TRIM(ADJUSTL(iString)), iStatus=iStatus) 
    end if 
 
    ! Close the Program specifications tag 
    if ( iStatus == 0 ) call Tbrook_CloseXMLTag(B=B, XMLTag="PROGRAMSPECIFICATIONS",iStatus=iStatus) 
 
  END SUBROUTINE TBU_WRITEPROGRAMSPECIFICATIONS 


  SUBROUTINE TBU_WRITEPROBES(B, iStatus)

    use parallel_info_module,   only: p_info 
    use probe_module,           only: probes
    use parameter_module,       only: nprobes
    use file_utility,           only: Make_File_Name
    use time_step_module,       only: t, cycle_number
    use tbrook_module,          only: BROOK,               & 
                                      TBrook_Set,          &
                                      TBrook_WriteXMLTag,  & 
                                      TBrook_OpenXMLTag,   & 
                                      TB_SCOPE_LOCAL,      &
                                      B_Strlen,            & 
                                      TBrook_CloseXMLTag 
    use output_data_module,     only: XML_Data_Format 
    use diagnostics_module,     only: PROBES_POSITIONS                   

    type(Brook), target :: B
    integer :: i, iStatus

    ! Local Variables 
    character (LEN=B_Strlen) :: string
    character(LEN=256) :: filename, probename, varname, thesuffix
    integer :: j, count, scalarsize, vectorsize, tensorsize
    real(r8), pointer, dimension(:) :: probe_cycle, probe_cycleV, probe_cycleT

100 FORMAT('NAME="',a,'"')
101 FORMAT('ID="',i8,'" X="',1es12.5,'" Y="',1es12.5,'" Z="',1es12.5, '"' )
102 FORMAT('NAME="',a,'" X="',1es12.5,'" Y="',1es12.5,'" Z="',1es12.5, '"' )
103 FORMAT(i3.3) 

    if ( iStatus /= 0 ) return 

    call PROBES_POSITIONS ()

    do i=1,nprobes

       write(string,102) TRIM(probes(i)%name),probes(i)%coords(1),probes(i)%coords(2),probes(i)%coords(3) 
       if ( iStatus == 0 ) then 
          call TBrook_OpenXMLTag(B = B,            &
               XMLTag              = "PROBE",      &
               XMLAttributes       = TRIM(string), &
               iStatus             = iStatus) 
       end if

       if ( iStatus == 0 ) then 
          call TBrook_WriteXMLTag(B = B,                           &
               XMLTag               = "DESCRIPTION",               &
               XMLStringData        = trim(probes(i)%description), &
               iStatus              = iStatus) 
       end if

       write(string,101)probes(i)%cell%index,probes(i)%cell%coords(1),probes(i)%cell%coords(2),probes(i)%cell%coords(3)

       if ( iStatus == 0 ) then 
          call TBrook_WriteXMLTag(B = B,            &
               XMLTag               = "CELL",       &
               XMLAttributes        = trim(string), &
               iStatus              = iStatus) 
       end if

       write(string,101)probes(i)%node%index,probes(i)%node%coords(1),probes(i)%node%coords(2),probes(i)%node%coords(3)

       if ( iStatus == 0 ) then 
          call TBrook_WriteXMLTag(B = B,            &
               XMLTag               = "NODE",       &
               XMLAttributes        = trim(string), &
               iStatus              = iStatus) 
       end if

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

          varname           = probes(i)%NameLU(count)      
          write(string,100) TRIM(varname) 

          if ( iStatus == 0 ) then 
             call TBrook_OpenXMLTag(B = B,            &
             XMLTag                   = "PROBEVAR",   & 
             XMLAttributes            = TRIM(string), &
             iStatus                  = iStatus) 
          end if
          write(string,103) i

          probename         = TRIM(string) // '_' // TRIM(varname)
          thesuffix         = TRIM(probename) // '.' // TRIM('bin') 
          filename          = MAKE_FILE_NAME (string='pb', suffix=thesuffix)

          call Tbrook_Set(probes(i)%BrookLU(count),    &
               File    = TRIM(filename),               & 
               Form    = "binary",                     &
               iStatus = iStatus) 

          ! Write look aside binary file tag 
          if ( iStatus == 0 ) then 
             call TBrook_WriteXMLTag(B = B,                                                & 
                  XMLTag               = "FILE",                                           & 
                  XMLAttributes        = 'FORMAT="'//trim(ADJUSTL(XML_Data_Format))//'"',  & 
                  XMLStringData        = TRIM(filename),                                   & 
                  iStatus              = iStatus) 
          end if

          probe_cycle(3) = probes(i)%ScalarVarLU(j)%field

          if ( iStatus == 0 .and. p_info%iop) then
             call tbu_make_file_entry(B, probes(i)%BrookLU(count), varname, &
                  probe_cycle, iStatus, .true.,'NoMesh','Attr="none"', TB_SCOPE_LOCAL) 
          end if

          if ( iStatus == 0 ) then 
             call TBrook_CloseXMLTag(B = B,          &
                  XMLTag               = "PROBEVAR", &
                  iStatus              = iStatus)
          end if
          
          count = count + 1

       end do

       do j=1,SIZE(probes(i)%VectorVarLU)

          varname        = probes(i)%NameLU(count)
          write(string,100) TRIM(varname) 
          if ( iStatus == 0 ) then 
             call TBrook_OpenXMLTag(B = B,          &
                  XMLTag              = "PROBEVAR", &
                  XMLAttributes       =TRIM(string),&
                  iStatus=iStatus) 
          end if
          write(string,103) i

          probename           = TRIM(string) // '_' // TRIM(varname)
          thesuffix           = TRIM(probename) // '.' // TRIM('bin') 
          filename            = MAKE_FILE_NAME (string='pb', suffix=thesuffix)

          call Tbrook_Set(probes(i)%BrookLU(count), &
               File   = TRIM(filename),             &
               Form   = "binary",                   &
               iStatus=iStatus) 

          ! Write look aside binary file tag 
          if ( iStatus == 0 ) then 
             call TBrook_WriteXMLTag(B = B,                                               & 
                  XMLTag               = "FILE",                                          & 
                  XMLAttributes        = 'FORMAT="'//trim(ADJUSTL(XML_Data_Format))//'"', & 
                  XMLStringData        = TRIM(filename),                                  & 
                  iStatus              = iStatus) 
          end if

          probe_cycleV(3:vectorsize) = probes(i)%VectorVarLU(j)%field
          
          if ( iStatus == 0 .and. p_info%iop) then
             call tbu_make_file_entry(B, probes(i)%BrookLU(count), varname, &
                  probe_cycleV, iStatus, .true.,'NoMesh','Attr="none"', TB_SCOPE_LOCAL) 
          end if

          if ( iStatus == 0 ) then 
             call TBrook_CloseXMLTag(B = B,          &
                  XMLTag               = "PROBEVAR", &
                  iStatus              = iStatus)
          end if
          
          count = count + 1

       end do

       do j=1,SIZE(probes(i)%TensorVarLU)

          varname        = probes(i)%NameLU(count)
          write(string,100) TRIM(varname) 
          if ( iStatus == 0 ) then 
             call TBrook_OpenXMLTag(B = B,           &
                  XMLTag              = "PROBEVAR",  &
                  XMLAttributes       = TRIM(string),&
                  iStatus             = iStatus) 
          end if
          write(string,103) i

          probename           = TRIM(string) // '_' // TRIM(varname)
          thesuffix           = TRIM(probename) // '.' // TRIM('bin') 
          filename            = MAKE_FILE_NAME (string='pb', suffix=thesuffix)

          call Tbrook_Set(probes(i)%BrookLU(count), File=TRIM(filename), Form="binary", iStatus=iStatus) 

          ! Write look aside binary file tag 
          if ( iStatus == 0 ) then 
             call TBrook_WriteXMLTag(B = B,                                               & 
                  XMLTag               = "FILE",                                          & 
                  XMLAttributes        = 'FORMAT="'//trim(ADJUSTL(XML_Data_Format))//'"', & 
                  XMLStringData        = TRIM(filename),                                  & 
                  iStatus              = iStatus) 
          end if

          probe_cycleT(3:tensorsize) = probes(i)%TensorVarLU(j)%field
          
          if ( iStatus == 0 .and. p_info%iop) then
             call tbu_make_file_entry(B, probes(i)%BrookLU(count), varname, &
                  probe_cycleT, iStatus, .true.,'NoMesh','Attr="none"', TB_SCOPE_LOCAL) 
          end if

          if ( iStatus == 0 ) then 
             call TBrook_CloseXMLTag(B = B,          &
                  XMLTag               = "PROBEVAR", &
                  iStatus              = iStatus)
          end if
          
          count = count + 1

       end do


       if (ASSOCIATED(probe_cycle))  DEALLOCATE(probe_cycle)
       if (ASSOCIATED(probe_cycleV)) DEALLOCATE(probe_cycleV)
       if (ASSOCIATED(probe_cycleT)) DEALLOCATE(probe_cycleT)

       if ( iStatus == 0 ) call TBrook_CloseXMLTag(B=B, XMLTag="PROBE",iStatus=iStatus)
          
    end do

  END SUBROUTINE TBU_WRITEPROBES


  SUBROUTINE TBU_WriteSimulationInformation (B, iStatus) 
    !--------------------------------------------------------------------------- 
    ! Purpose: 
    ! 
    !   print build time and run time information 
    !--------------------------------------------------------------------------- 
    use EM_data_proxy,          only: EM_is_on 
    use fluid_data_module,      only: fluid_flow 
    use solid_mechanics_data,   only: solid_mechanics 
    use output_module,          only: title 
    use tbrook_module,          only: BROOK,               & 
                                      TBrook_WriteXMLTag,  & 
                                      TBrook_OpenXMLTag,   & 
                                      B_Strlen,            & 
                                      TBrook_CloseXMLTag 
    use mesh_input_module,     only:  coordinate_scale_factor
 
    use parallel_info_module,   only: p_info 
    use diffusion_solver_data,  only: ds_enabled, num_species
 
    type(Brook), target :: B 
    integer :: iStatus
 
    ! Local Variables 
    character (LEN=B_Strlen) :: string
    character (LEN=B_Strlen) :: rstart_common, & 
                                rstart_micro,  & 
                                rstart_flow,   & 
                                rstart_ht,     & 
                                rstart_em,     & 
                                rstart_sm 
 
    !--------------------------------------------------------------------------- 
    ! collect info we don't already have 
 
 
1   FORMAT('NAME="',a,'" VALUE="',1L1,'"') 
2   Format('TSTRLEN="',i10,'"') 
    if ( iStatus /= 0 ) return 
 
    if(p_info%IOP) then 
       ! In the rstart lines, a comma indicates a newline 
       rstart_common = 'CODE, VERSION, PHYSICS, LIBRARIES, BUILDDATETIME, BUILDARCHITECTURE, BUILDFLAGS, BUILDHOST, ' //  & 
            'RUNARCHITECTURE, RUNHOST, RUNDATE, TITLE, t, dt, dt_courant, cycle, ncells, nnodes, boundaryfaces, '//    & 
            'MAT_SLOT, NPHASES, NCOMPONENTS, CELLVERTICES, VERTEXCOORDS, RSUMRVOL, '//                                 & 
            'Z_RHO, Z_TEMP, Z_TEMP_OLD, Z_ENTHALPY, Z_ENTHALPYOLD, Z_P, ' //                                 & 
            'Z_VC, Z_VCOLD, Z_VF, '   //                                                                               & 
            'M_ID, M_VOF, M_VOF_OLD ' 
       rstart_micro  = 'LIQUIDUS_TEMP, EUTECTIC_FRACTION, SEC_ARM_SPACING, PHASECONC, PHASECONCOLD' 
       rstart_flow   = ' ' 
       rstart_ht     = ' ' 
       rstart_sm     = 'TOTAL_STRAIN, ELASTIC_STRESS, PLASTIC_STRAIN, PLASTIC_STRAIN_RATE, RHS, epstherm, epspc, DISPLACEMENT' 
       rstart_em     = ' ' 


       if ( iStatus == 0 ) call Tbrook_OpenXMLTag(B=B, XMLTag="SIMULATIONINFORMATION",iStatus=iStatus) 

       write(string,2) LEN(title) 
       if ( iStatus == 0 ) call Tbrook_WriteXMLTag(B=B, XMLTag="TITLE",XMLAttributes=TRIM(string), & 
            XMLStringData=TRIM(title), iStatus=iStatus) 
       if ( iStatus == 0 ) call Tbrook_OpenXMLTag (B=B, XMLTag="RESTARTVARIABLES",iStatus=iStatus) 
       if ( iStatus == 0 ) call TBrook_WriteXMLTag(B=B, XMLTag='RESTART_COMMON', XMLStringData=trim(rstart_common), &
                                                   iStatus=iStatus) 
       if ( iStatus == 0 ) call TBrook_WriteXMLTag(B=B, XMLTag='RESTART_MICRO', XMLStringData=trim(rstart_micro),   &
                                                   iStatus=iStatus) 
       if ( iStatus == 0 ) call TBrook_WriteXMLTag(B=B, XMLTag='RESTART_FLOW', XMLStringData=trim(rstart_flow),     &
                                                   iStatus=iStatus) 
       if ( iStatus == 0 ) call TBrook_WriteXMLTag(B=B, XMLTag='RESTART_HEAT_TRANSFER', XMLStringData=trim(rstart_ht), &
                                                   iStatus=iStatus) 
       if ( iStatus == 0 ) call TBrook_WriteXMLTag(B=B, XMLTag='RESTART_SOLID_MECHANICS', XMLStringData=trim(rstart_sm), &
                                                   iStatus=iStatus) 
       if ( iStatus == 0 ) call TBrook_WriteXMLTag(B=B, XMLTag='RESTART_ELECTROMAGNETICS', XMLStringData=trim(rstart_em), &
                                                   iStatus=iStatus) 
       if ( iStatus == 0 ) call Tbrook_CloseXMLTag(B=B, XMLTag="RESTARTVARIABLES",iStatus=iStatus) 

       !! NNC, 2 May 2012.  Not sure what to do here.  HEAT_CONDUCTION now applies to new solver
       !! but the output parser has understood this to be old solver.  How is this used by the parser?
       !! For know I write out false just to be safe.
       write(string,1) 'HEAT_CONDUCTION', .false. !heat_conduction 
       if ( iStatus == 0 ) call Tbrook_WriteXMLTag(B=B, XMLTag="PHYSICS",XMLAttributes=TRIM(string), iStatus=iStatus) 

       write(string,1) 'FLUID_FLOW', fluid_flow 
       if ( iStatus == 0 ) call Tbrook_WriteXMLTag(B=B, XMLTag="PHYSICS",XMLAttributes=TRIM(string), iStatus=iStatus) 

       write(string,1) 'SOLID_MECHANICS', solid_mechanics 
       if ( iStatus == 0 ) call Tbrook_WriteXMLTag(B=B, XMLTag="PHYSICS",XMLAttributes=TRIM(string), iStatus=iStatus) 

       write(string,1) 'ELECTROMAGNETICS', EM_is_on() 
       if ( iStatus == 0 ) call Tbrook_WriteXMLTag(B=B, XMLTag="PHYSICS",XMLAttributes=TRIM(string), iStatus=iStatus) 

       !! NNC, 2 May 2012. I don't know if I can just delete this.  Does the parser need this?
3      FORMAT('NPHASES="',I10,'" NCOMPONENTS="',I10,'"') 
       write(string,3) -1,-1 
       if ( iStatus == 0 ) call Tbrook_WriteXMLTag(B=B, XMLTag="ALLOYINFORMATION",XMLAttributes=TRIM(string), iStatus=iStatus) 

4      FORMAT('NSPECIES="',I10,'"') 
       if (ds_enabled) then
          write(string,4) num_species
       else
          write(string,4) -1
       end if
       if ( iStatus == 0 ) call Tbrook_WriteXMLTag(B=B, XMLTag="SPECIESINFO",XMLAttributes=TRIM(string), iStatus=iStatus) 

       write(string,*) coordinate_scale_factor
       if ( iStatus == 0 ) call Tbrook_WriteXMLTag(B=B, XMLTag="MESHSCALEFACTOR",XMLStringData=trim(adjustl(string)), &
                                                   iStatus=iStatus) 

       if ( iStatus == 0 ) call Tbrook_CloseXMLTag(B=B, XMLTag="SIMULATIONINFORMATION",iStatus=iStatus) 
    end if 
 
  END SUBROUTINE TBU_WRITESIMULATIONINFORMATION 
 
  SUBROUTINE TBU_WriteBasicData(B, Copyright, Disclaimer, iStatus) 
    ! 
    ! Purpose: Write the basic data block to match that in the design document. 
    ! This includes: Copyright Information 
    !                Disclaimer Information 
    !                Truchas Information 
    !                Timestamp Information 
    !                Machine Information 
    !                Mesh Information 
    ! 
    ! Information is written to brook B 
    ! 
    use file_utility,        only:  Make_File_Name 
    use output_data_module,  only:  retain_last_step
    use tbrook_module,       only:  Brook,               & 
                                    TBrook_Write,        & 
                                    TBrook_WriteXMLTag,  & 
                                    TBrook_OpenXMLTag,   & 
                                    TBrook_CloseXMLTag,  & 
                                    TBrook_Set,          & 
                                    ASSIGNMENT(=),       & 
                                    TB_SCOPE_LOCAL,      & 
                                    TBrook_Destroy 
 
    use parallel_info_module, only: p_info 
 
    type(Brook),      target         :: B 
    character(LEN=*), dimension(:)   :: Copyright 
    character(LEN=*), dimension(:)   :: Disclaimer 
    integer,          intent(INOUT)  :: iStatus 
    character(LEN=70), dimension(22) :: helpstr=(/ & 
                           "<![CDATA[                                                        ", & 
                           "                                                                 ", & 
                           "  In case the run dies in the middle, you can restore this file  ", & 
                           "  by closing any open tags.  Typically, the only tags you will   ", & 
                           "  have to close are CYCLE and TruchasData:                       ", & 
                           "         </CYCLE>                                                ", & 
                           "     </TruchasData>                                              ", & 
                           "                                                                 ", & 
                           "  or TIMESTEP and TruchasData:                                   ", & 
                           "         </TIMESTEP>                                             ", & 
                           "     </TruchasData>                                              ", & 
                           "                                                                 ", & 
                           "  Most other tags write their information in a single pass and   ", & 
                           "  should already be closed.                                      ", & 
                           "                                                                 ", & 
                           "  It is also possible that some of your variable data in the     ", & 
                           "  lookaside files '*.ts.#####.xml' are incomplete, but there is  ", & 
                           "  no way to recover this data since it was never written.  The   ", & 
                           "  OutputParser.py should ignore these files, but if it doesn't,  ", & 
                           "  You can just delete the last timestep tag from this file.      ", & 
                           "                                                                 ", & 
                           "]]>                                                              "  & 
                           /) 
 
    type(Brook), target :: myB
    character(LEN=256)  :: LastStepFile
 
 
    if ( iStatus /= 0 ) return 
 
    ! Copy brook for personal use 
    myB = B 
    call Tbrook_Set(B=myB, form="ascii", iStatus=iStatus) 
 
    ! Open the BasicData Tag 
    if ( iStatus == 0 ) call Tbrook_OpenXMLTag(B=myB, XMLTag="BASICDATA",iStatus=iStatus) 
 
    if ( p_info%IOP) then 
       if ( iStatus == 0 ) call TBrook_Write(B=myB, Variable=helpstr, SCOPE=TB_SCOPE_LOCAL, iStatus=iStatus) 
 
       ! Machine Information, note use of B and not myB 
       if ( iStatus == 0 ) call TBU_WriteProgramSpecifications(B=B, iStatus=iStatus) 
    end if 
 
 
    ! Simulation Information, note use of B and not myB 
    call TBU_WriteSimulationInformation(B=B, iStatus=iStatus) 

    !Write out tags for all probes
    call TBU_WRITEPROBES(B=B,iStatus=iStatus)
 
    !Write out last timestep entry in the *.TBrook.xml file (if retain_last_step = .true.)

    if ((p_info%IOP).and.(retain_last_step)) then
       if ( iStatus == 0 ) then 
          call Tbrook_OpenXMLTag(B = myB,        & 
               XMLTag        = "LASTSTEP",        &
               SCOPE         = TB_SCOPE_LOCAL,   &
               iStatus       = iStatus)
       end if

       LastStepFile = MAKE_FILE_NAME (string=TRIM('LastStep'), suffix='xml') 

       !Write look aside tag
       if ( iStatus == 0 ) then 
          call Tbrook_WriteXMLTag(B = myB,                  & 
               XMLTag        = "FILE",                      & 
               XMLAttributes = 'FORMAT="xml"',              & 
               XMLStringData = TRIM(LastStepFile),          & 
               SCOPE         = TB_SCOPE_LOCAL,              & 
               iStatus       = iStatus) 
       end if

       if ( iStatus == 0 ) then 
          call Tbrook_CloseXMLTag(B=myB, XMLTag="LASTSTEP",SCOPE=TB_SCOPE_LOCAL,iStatus=iStatus) 
       end if

    end if
 
    if ( p_info%IOP) then 
       ! Note that we are writing string arrays here, so we have to set the scope as local 
       ! Copyright information 
       if ( iStatus == 0 ) call Tbrook_OpenXMLTag(B=myB, XMLTag="COPYRIGHT",SCOPE=TB_SCOPE_LOCAL,iStatus=iStatus) 
       if ( iStatus == 0 ) call Tbrook_Write(B=myB, Variable="<![CDATA[",ADVANCE=.true., SCOPE=TB_SCOPE_LOCAL,iStatus=iStatus) 
       if ( iStatus == 0 ) call Tbrook_Write(B=myB, Variable=Copyright,SCOPE=TB_SCOPE_LOCAL,iStatus=iStatus) 
       if ( iStatus == 0 ) call Tbrook_Write(B=myB, Variable="]]>",ADVANCE=.true., SCOPE=TB_SCOPE_LOCAL,iStatus=iStatus) 
       if ( iStatus == 0 ) call Tbrook_CloseXMLTag(B=myB, XMLTag="COPYRIGHT",SCOPE=TB_SCOPE_LOCAL,iStatus=iStatus) 
 
       ! Disclaimer information 
       if ( iStatus == 0 ) call Tbrook_OpenXMLTag(B=myB, XMLTag="DISCLAIMER",SCOPE=TB_SCOPE_LOCAL,iStatus=iStatus) 
       if ( iStatus == 0 ) call Tbrook_Write(B=myB, Variable="<![CDATA[",ADVANCE=.true., SCOPE=TB_SCOPE_LOCAL,iStatus=iStatus) 
       if ( iStatus == 0 ) call Tbrook_Write(B=myB, Variable=Disclaimer,SCOPE=TB_SCOPE_LOCAL,iStatus=iStatus) 
       if ( iStatus == 0 ) call Tbrook_Write(B=myB, Variable="]]>",ADVANCE=.true., SCOPE=TB_SCOPE_LOCAL,iStatus=iStatus) 
       if ( iStatus == 0 ) call Tbrook_CloseXMLTag(B=myB, XMLTag="DISCLAIMER",SCOPE=TB_SCOPE_LOCAL,iStatus=iStatus) 
    end if 
 
    call Tbrook_Set(B=myB, form="XML", iStatus=iStatus) 
 
    ! Close the BasicData Tag 
    if ( iStatus == 0 ) call Tbrook_CloseXMLTag(B=myB, XMLTag="BASICDATA",iStatus=iStatus) 
 
    call TBrook_Destroy(myB, iStatus=iStatus) 

  END SUBROUTINE TBU_WriteBasicData 


  SUBROUTINE TBU_WriteDefaultMesh (B, iStatus) 
    !--------------------------------------------------------------------------- 
    ! Purpose: 
    !   Write a GMV graphics dump 
    !--------------------------------------------------------------------------- 
    use tbrook_module,       only: BROOK, ASSIGNMENT(=)
    use mesh_module,          only: Mesh,                    & 
                                    Vertex,                  & 
                                    Permute_Mesh_Vector,     & 
                                    UnPermute_Mesh_Vector,   & 
                                    Permute_Vertex_Vector,   & 
                                    Unpermute_Vertex_Vector, & 
                                    Cell,                    &
                                    mesh_has_cblockid_data
    use parameter_module,     only: ncells, ndim, nfc, nnodes, nvc, & 
                                    ncells_tot, nnodes_tot, boundary_faces_tot 
    use parallel_info_module, only: p_info 
 
    ! Argument List 
    type(Brook), target, intent(in)  :: B
    integer :: iStatus 
 
    ! Local Variables 
    integer :: n, c, v, f 
 
    integer, dimension(nvc,ncells)   :: Vrtx_Ngbr 
    real(r8), dimension(ndim,nnodes) :: VertexCoords 
    real(r8), dimension(ndim,ncells) :: Centroids 
    integer, dimension(ncells) :: NNeighbors 
    integer, dimension(ncells) :: Boundary_Flag 
    integer, dimension(ncells) :: cpart
    integer, dimension(nnodes) :: vpart
    !--------------------------------------------------------------------------- 
 
    ! our new and improved writer? 
    iStatus = 0 
    do n = 1, ndim 
       Centroids(n,:)    = Cell(:)%Centroid(n) 
       VertexCoords(n,:) = Vertex(:)%Coord(n) 
    end do 
 
    do c = 1, ncells 
       NNeighbors(c) = SIZE(Mesh(c)%Ngbr_Cells_All%v(:)) 
    end do 
    do v = 1,nvc 
       Vrtx_Ngbr(v,:) = Mesh%Ngbr_Vrtx_Orig(v) 
    end do 
 
    Boundary_Flag = 0 
    do f = 1, nfc 
       Where (Mesh%Ngbr_cell(f) == 0) Boundary_Flag = 1 
    end do 
 
    cpart = p_info%thisPE
    vpart = p_info%thisPE

    if (mesh_has_cblockid_data) then
    call TBU_MeshWriter(B,                                             &
                        MLabel              = 'DefaultMesh',           & 
                        MBoundary_faces_tot = boundary_faces_tot,      & 
                        CN_tot              = ncells_tot,              & 
                        CBlockID            = Mesh%CBlockID,           &
                        CPartition          = cpart,                   &
                        CBoundary           = Boundary_Flag,           & 
                        CCentroids          = Centroids,               & 
                        CNNeighbors         = NNeighbors,              & 
                        CPermute            = Permute_Mesh_Vector,     & 
                        CUnpermute          = Unpermute_Mesh_Vector,   & 
                        CVertices           = Vrtx_Ngbr,               & 
                        CVolumes            = Cell(:)%Volume,          & 
                        VN_tot              = nnodes_tot,              & 
                        VCoords             = VertexCoords,            & 
                        VPartition          = vpart,                   &
                        VPermute            = Permute_Vertex_Vector,   & 
                        VRSum_RVol          = Vertex(:)%Rsum_RVol,     & 
                        VUnpermute          = UnPermute_Vertex_Vector, & 
                        iStatus             = iStatus) 
    else
       call TBU_MeshWriter(B,                                             &
                           MLabel              = 'DefaultMesh',           &
                           MBoundary_faces_tot = boundary_faces_tot,      &
                           CN_tot              = ncells_tot,              &
                           CPartition          = cpart,                   &
                           CBoundary           = Boundary_Flag,           &
                           CCentroids          = Centroids,               &
                           CNNeighbors         = NNeighbors,              &
                           CPermute            = Permute_Mesh_Vector,     &
                           CUnpermute          = Unpermute_Mesh_Vector,   &
                           CVertices           = Vrtx_Ngbr,               &
                           CVolumes            = Cell(:)%Volume,          &
                           VN_tot              = nnodes_tot,              &
                           VCoords             = VertexCoords,            &
                           VPartition          = vpart,                   &
                           VPermute            = Permute_Vertex_Vector,   &
                           VRSum_RVol          = Vertex(:)%Rsum_RVol,     &
                           VUnpermute          = UnPermute_Vertex_Vector, &
                           iStatus             = iStatus)
    end if
 
  END SUBROUTINE TBU_WriteDefaultMesh 
 
 
 
  SUBROUTINE TBU_MeshWriter (                    & 
                             B,                  & ! Main xml file brook
                             MLabel,             & !   Mesh label 
                             MBoundary_faces_tot,& ! O Number of boundary faces 
 
                             CN_Tot,             & !   Total number of cells ( across all processors) 
                             CVertices,          & !   Vertices for a given cell 
                             CBoundary,          & ! O Is it a boundary cell 
                             CBlockID,           & ! O Cell block IDS 
                             CPartition,         & ! O Cell Partition 
                             CCentroids,         & ! O Cell Centroids 
                             CNNeighbors,        & ! O Cell Neighbors 
                             CPermute,           & ! O Cell permute from User to T 
                             CUnpermute,         & ! O Cell permute from T to User 
                             CVolumes,           & ! O Cell Volumes 
 
                             FN_Tot,             & ! O Total number of faces 
                             FBoundary,          & ! O Is it a boundary face 
                             FPartition,         & ! O Which partition the face is on 
                             FCentroids,         & ! O Face Centroids 
                             FCells,             & ! O Cells on either side of face 
                             FEdges,             & ! O Edges on face 
                             FVertices,          & ! O Vertices for face 
                             FPermute,           & ! O Face permute from user to T 
                             FUnpermute,         & ! O Face permute from T to User 
 
                             EN_Tot,             & ! O Total number of edges 
                             EPartition,         & ! O WHich partition the edge is on 
                             EBoundary,          & ! O Is it a boundary edge 
                             ECentroids,         & ! O Edge Centroid 
                             ECells,             & ! O Cells the edge borders 
                             EFaces,             & ! O Faces touching edge 
                             EVertices,          & ! O Vertices on edge 
                             EPermute,           & ! O Edge permute from User to T 
                             EUnpermute,         & ! O Edge permute from T to User 
 
                             VN_tot,             & !   Total number of vertices 
                             VCoords,            & !   Vertex Coordinates 
                             VPartition,         & ! O Which partition the vertices are on 
                             VBoundary,          & ! O Is it a boundary vertex 
                             VRSum_RVol,         & ! O Volume associated with vertex 
                             VPermute,           & ! O Permute Vertex from User to T 
                             VUnpermute,         & ! O Permute Vertex from T to User 
 
                             Scope,              & 
                             iStatus) 
    !--------------------------------------------------------------------------- 
    ! Purpose: 
    !   Write a mesh file as a brookxml mesh 
    !--------------------------------------------------------------------------- 
    use file_utility,        only: Make_File_Name 
    use tbrook_module,       only: TBrook_WriteXMLTag,  & 
                                   TBrook_OpenXMLTag,   & 
                                   TBrook_CloseXMLTag,  & 
                                   TBrook_Set,          & 
                                   TBrook_Get_Scope,    & 
                                   BROOK,               & 
                                   TBrook_Close,        & 
                                   TBrook_File,         & 
                                   TBrook_Destroy,      &
                                   TB_SCOPE_LOCAL
    use parallel_info_module, only: p_info 
 
    ! Argument List 
    type(Brook), target, intent(in) :: B  ! Main xml file brook
    
    character(LEN=*),  intent(in) :: MLabel          ! Label of Mesh 
    integer, optional, intent(in) :: Mboundary_faces_tot 
 
    integer,                            intent(in) :: CN_tot       ! Number of cells 
    integer,  dimension(:,:),           intent(in) :: CVertices    ! Vertex IDs for each cell 
    real(r8), dimension(:),   optional, intent(in) :: CVolumes     ! Volumes of cells 
    real(r8), dimension(:,:), optional, intent(in) :: CCentroids   ! Centroids of cells 
    integer,  dimension(:),   optional, intent(in) :: CBlockID     ! Block ID for cell 
    integer,  dimension(:),   optional, intent(in) :: CPartition   ! Partition on which cell sits 
    integer,  dimension(:),   optional, intent(in) :: CBoundary    ! Boudnary cell? 0 = no, 1 = yes 
    integer,  dimension(:),   optional, intent(in) :: CNNeighbors  ! Number of neighbors 
    integer,  dimension(:),   optional, intent(in) :: CPermute     ! Permute Cells from User to T 
    integer,  dimension(:),   optional, intent(in) :: CUnpermute   ! Permute Cells from T to User 
 
    integer,                            intent(in) :: VN_tot       ! Number of vertices 
    real(r8), dimension(:,:),           intent(in) :: VCoords      ! Coordinates of Vertices 
    real(r8), dimension(:),   optional, intent(in) :: VRsum_RVol   ! Sum divided by volume??? 
    integer,  dimension(:),   optional, intent(in) :: VPartition   ! Partition on which Vertex sits 
    integer,  dimension(:),   optional, intent(in) :: VBoundary    ! Boudnary Vertex? 0 = no, 1 = yes 
    integer,  dimension(:),   optional, intent(in) :: VPermute     ! Permute Vertices from User to T 
    integer,  dimension(:),   optional, intent(in) :: VUnpermute   ! Unpermute Vertices from T to User 
 
    integer,                  optional, intent(in) :: FN_tot       ! Number of Faces 
    integer,  dimension(:),   optional, intent(in) :: FBoundary    ! Boudnary Face? 0 = no, 1 = yes 
    integer,  dimension(:),   optional, intent(in) :: FPartition   ! Partition on which Face sits 
    integer,  dimension(:),   optional, intent(in) :: FPermute     ! Permute Faces from User to T 
    integer,  dimension(:),   optional, intent(in) :: FUnpermute   ! Permute Faces from T to User 
    real(r8), dimension(:,:), optional, intent(in) :: FCentroids   ! Centroids of Faces 
    integer,  dimension(:,:), optional, intent(in) :: FVertices    ! Vertex IDs for each Face 
    integer,  dimension(:,:), optional, intent(in) :: FCells       ! Cell IDs for each Face 
    integer,  dimension(:,:), optional, intent(in) :: FEdges       ! Edge IDs for each Face 
 
    integer,                  optional, intent(in) :: EN_tot       ! Number of Edges 
    integer,  dimension(:),   optional, intent(in) :: EBoundary    ! Boudnary Edge? 0 = no, 1 = yes 
    integer,  dimension(:),   optional, intent(in) :: EPartition   ! Partition on which Edge sits 
    integer,  dimension(:),   optional, intent(in) :: EPermute     ! Permute Edges from User to T 
    integer,  dimension(:),   optional, intent(in) :: EUnpermute   ! Permute Edges from T to User 
    real(r8), dimension(:,:), optional, intent(in) :: ECentroids   ! Centroids of Edge 
    integer,  dimension(:,:), optional, intent(in) :: EVertices    ! Vertex IDs for each Edge 
    integer,  dimension(:,:), optional, intent(in) :: ECells       ! Cell IDs for each Edge 
    integer,  dimension(:,:), optional, intent(in) :: EFaces       ! Face IDs for each Edge 
 
 
 
    integer, optional, intent(in)  :: Scope     ! Parallel or serial 
    integer,           intent(out) :: iStatus   ! Did an error occur? 
 
    !Local Variables 
    character(LEN=1024)  :: tmpString 
    character(LEN=1024)  :: ErrString 
    integer              :: error_code 
    integer              :: nvc 
    integer              :: myScope 
    type(Brook), target  :: BAside 
    type(Brook), target  :: BaseBrook 
    character(LEN=B_STRLEN) :: MeshFile
    integer, save :: MeshN = 0
    !--------------------------------------------------------------------------- 
 
    MeshN = MeshN + 1 
 
638 Format('FORMAT="TRUCHASMESH" LABEL="',a,'" MESHWRITERVERSION="',a,'"') 
    myScope = TBrook_Get_Scope(scope=scope) 
 
 
    error_code = 0 

    ! Now write the Mesh file 
    if ( error_code == 0 ) then 
       MeshFile = MAKE_FILE_NAME (string=TRIM(MLabel), number=MeshN, suffix='xml') 
       ErrString = 'Unable to create Mesh XML File' 
       call Tbrook_Set(BaseBrook, File=TRIM(MeshFile), Form="xml", iStatus=error_code) 
    end if 
 
    if ( error_code == 0 ) then 
       ErrString = 'Unable to create Mesh Binary File' 
       tmpString = MAKE_FILE_NAME (string=TRIM(MLabel), number=MeshN, suffix='bin') 
       call Tbrook_Set(BAside, File=TRIM(tmpString), Form="binary", iStatus=error_code) 
    end if 
 
    ! Open Mesh Tag 
    write(tmpString,638) TRIM(MLabel), TRIM(TBU_Version) 
    if ( error_code == 0 ) then 
       ErrString = 'Unable to open mesh tag in lookaside' 
       call tbrook_OpenXMLTag(BaseBrook,                   & 
                              XMLTag="MESH",           & 
                              XMLAttributes=tmpString, & 
                              SCOPE=myScope,           & 
                              iStatus=error_code) 
    end if 
 
    ! Write Look Aside Tag 
    if ( error_code == 0 ) then 
       ErrString = 'Error writing binary mesh file name to xml file' 
       call tbrook_WriteXMLTAG(BaseBrook,                                  & 
                               XMLTag="FILE",                          & 
                               XMLAttributes='FORMAT="binary"',        & 
                               XMLStringData=TRIM(TBrook_File(BAside)), & 
                               SCOPE=myScope,                          & 
                               iStatus=error_code) 
    end if 
    ! Write Mesh Type 
    nvc = size(CVertices,1) 
    if ( nvc == 4 ) then 
       tmpString='TET' 
    else 
       tmpString='HEX'
    end if
    if ( error_code == 0 ) then
       ErrString = 'Error writing mesh type to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'MESHTYPE', trim(tmpString), error_code, &
                                .true., MLabel, 'INFORMATION="MESH"', myScope, 'NONE')
    end if

    if (error_code == 0 ) then
       ErrString = 'Error writing NCELLS type to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'N', CN_Tot, error_code, &
                                .true., MLabel, 'INFORMATION="CELL"', myScope, 'CELLS')
    end if

    if ( error_code == 0 .and. present(CPermute)) then
       ErrString = 'Error Writing PermuteCell'
       call Tbu_make_file_entry(BaseBrook, BAside, 'PERMUTE', CPERMUTE, error_code, &
                                .true., MLabel,  'INFORMATION="CELL"', myScope, 'CELLS')
    end if

    if (error_code == 0 .and. present(CUnpermute)) then
       ErrString = 'Error Writing UnpermuteCell'
       call Tbu_make_file_entry(BaseBrook, BAside, 'UNPERMUTE', CUNPERMUTE, error_code, &
                                .true., MLabel,  'INFORMATION="CELL"', myScope, 'CELLS')
    end if

    if ( error_code == 0 ) then
       ErrSTring = 'Error writing CellData to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'VERTICES', CVertices, error_code, &
                                .true., MLabel,  'INFORMATION="CELL"', myScope, 'CELLS')
     end if


     if (error_code == 0 .and. present(CVolumes)) then
        ErrString = 'Error writing CellVolumes to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'VOLUMES', CVolumes, error_code, &
                                .true., MLabel,  'INFORMATION="CELL"', myScope, 'CELLS')
     end if

     if (error_code == 0 .and. present(CCentroids)) then
        ErrString = 'Error writing CellCentroids to mesh file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'CENTROIDS', CCentroids, error_code, &
                                .true., MLabel,  'INFORMATION="CELL"', myScope, 'CELLS')
     end if

     if (error_code == 0 .and. present(CBoundary)) then
        ErrSTring = 'Error writing CellBoundary to xml file'
        call Tbu_make_file_entry(BaseBrook, BAside, 'BOUNDARY', CBoundary, error_code, &
                                 .true., MLabel,  'INFORMATION="CELL"', myScope, 'CELLS')
     end if

     if (error_code == 0 .and. present(CPartition)) then
        ErrSTring = 'Error writing CELLPARTITION to xml file'
        call Tbu_make_file_entry(BaseBrook, BAside, 'PARTITIONS', CPartition, error_code, &
                                 .true., MLabel,  'INFORMATION="CELL"', myScope, 'CELLS')
     end if

     if (error_code == 0 .and. present(CBlockID)) then
        ErrSTring = 'Error writing CELLBLOCKID to xml file'
        call Tbu_make_file_entry(BaseBrook, BAside, 'BLOCKID', CBlockID, error_code, &
                                 .true., MLabel,  'INFORMATION="CELL"', myScope, 'CELLS')
     end if

     if (error_code == 0 .and. present(CNNeighbors)) then
        ErrString = 'Error writing CellNUMNEIGHBORS to xml file'
        call Tbu_make_file_entry(BaseBrook, BAside, 'NUMNEIGHBORS', CNNeighbors, error_code, &
                                 .true., MLabel,  'INFORMATION="CELL"', myScope, 'CELLS')
     end if

    if (error_code == 0 .and. present(FN_TOT)) then
       ErrString = 'Error writing NFACES type to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'N', FN_tot, error_code, &
                                .true., MLabel,  'INFORMATION="FACE"', myScope, 'FACES')
    end if

    if ( error_code == 0 .and. present(FPermute)) then
       ErrString = 'Error Writing PermuteFace'
       call Tbu_make_file_entry(BaseBrook, BAside, 'PERMUTE', FPERMUTE, error_code, &
                                .true., MLabel,  'INFORMATION="FACE"', myScope, 'FACES')
    end if

    if (error_code == 0 .and. present(FUnpermute)) then
       ErrString = 'Error Writing UnpermuteFace'
       call Tbu_make_file_entry(BaseBrook, BAside, 'UNPERMUTE', FUnPERMUTE, error_code, &
                                .true., MLabel,  'INFORMATION="FACE"', myScope, 'FACES')
    end if

    if ( error_code == 0 .and. present(FVertices) ) then
       ErrSTring = 'Error writing FaceData to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'VERTICES', FVertices, error_code, &
                                .true., MLabel,  'INFORMATION="FACE"', myScope, 'FACES')
     end if

    if ( error_code == 0 .and. present(FCells)) then
       ErrSTring = 'Error writing FaceData to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'CELLS', FCells, error_code, &
                                .true., MLabel,  'INFORMATION="FACE"', myScope, 'FACES')
     end if

    if ( error_code == 0 .and. present(FEdges)) then
       ErrSTring = 'Error writing FaceData to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'EDGES', FEdges, error_code, &
                                .true., MLabel,  'INFORMATION="FACE"', myScope, 'FACES')
    end if


    if (error_code == 0 .and. present(FCentroids)) then
       ErrString = 'Error writing FaceCentroids to mesh file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'CENTROIDS', FCentroids, error_code, &
                                .true., MLabel,  'INFORMATION="FACE"', myScope, 'FACES')
    end if

    if (error_code == 0 .and. present(FBoundary)) then
       ErrSTring = 'Error writing FaceBoundary to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'BOUNDARY', FBoundary, error_code, &
                                .true., MLabel,  'INFORMATION="FACE"', myScope, 'FACES')
    end if
    
    if (error_code == 0 .and. present(FPartition)) then
       ErrSTring = 'Error writing FACEPARTITION to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'PARTITIONS', FPartition, error_code, &
            .true., MLabel,  'INFORMATION="FACE"', myScope, 'FACES')
    end if
    
    if (error_code == 0 .and. present(Mboundary_faces_tot)) then
       ErrString = 'Error writing TOTALBOUNDARYFACES to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'TOTALBOUNDARYFACES', MBoundary_Faces_Tot, error_code, &
                                .true., MLabel,  'INFORMATION="FACE"', myScope, 'FACES')
    end if

    ! Edge based data
    if (error_code == 0 .and. present(EN_TOT)) then
       ErrString = 'Error writing NEDGES type to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'N', EN_Tot, error_code, &
                                .true., MLabel,  'INFORMATION="EDGE"', myScope, 'EDGES')
    end if

    if ( error_code == 0 .and. present(EPermute)) then
       ErrString = 'Error Writing PermuteEdge'
       call Tbu_make_file_entry(BaseBrook, BAside, 'PERMUTE', EPermute, error_code, &
                                .true., MLabel,  'INFORMATION="EDGE"', myScope, 'EDGES')
    end if

    if (error_code == 0 .and. present(EUnpermute)) then
       ErrString = 'Error Writing UnpermuteEdge'
       call Tbu_make_file_entry(BaseBrook, BAside, 'UNPERMUTE', EUnPermute, error_code, &
                                .true., MLabel,  'INFORMATION="EDGE"', myScope, 'EDGES')
    end if

    if ( error_code == 0 .and. present(EVertices)) then
       ErrSTring = 'Error writing EdgeData to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'VERTICES', EVERTICES, error_code, &
                                .true., MLabel,  'INFORMATION="EDGE"', myScope, 'EDGES')
     end if

    if ( error_code == 0 .and. present(ECells)) then
       ErrSTring = 'Error writing EdgeData to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'CELLS', ECELLS, error_code, &
                                .true., MLabel,  'INFORMATION="EDGE"', myScope, 'EDGES')
     end if

    if ( error_code == 0 .and. present(EFaces)) then
       ErrSTring = 'Error writing EdgeData to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'FACES', EFACES, error_code, &
                                .true., MLabel,  'INFORMATION="EDGE"', myScope, 'EDGES')
    end if

    
    if (error_code == 0 .and. present(ECentroids)) then
       ErrString = 'Error writing EdgeCentroids to mesh file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'CENTROIDS', ECentroids, error_code, &
                                .true., MLabel,  'INFORMATION="EDGE"', myScope, 'EDGES')
    end if

    if (error_code == 0 .and. present(EBoundary)) then
       ErrSTring = 'Error writing EdgeBoundary to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'BOUNDARY', EBoundary, error_code, &
                                .true., MLabel,  'INFORMATION="EDGE"', myScope, 'EDGES')
    end if

    if (error_code == 0 .and. present(EPartition)) then
       call Tbu_make_file_entry(BaseBrook, BAside, 'PARTITIONS', EPartition, error_code, &
                                .true., MLabel,  'INFORMATION="EDGE"', myScope, 'EDGES')
    end if

    ! Vertex based data
    if (error_code == 0 ) then
       ErrString = 'Error writing nnodes to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'N', VN_Tot, error_code, &
                                .true., MLabel,  'INFORMATION="VERTEX"', myScope, 'NODES')
    end if

    if ( present(VPermute) .and. error_code == 0) then
       ErrString = 'Error Writing PermuteVertex'
       call Tbu_make_file_entry(BaseBrook, BAside, 'PERMUTE', VPermute, error_code, &
                                .true., MLabel,  'INFORMATION="VERTEX"', myScope, 'NODES')
    end if

    if ( error_code == 0 .and. present(VUnpermute)) then
       ErrString = 'Error Writing UnpermuteVertex'
       call Tbu_make_file_entry(BaseBrook, BAside, 'UNPERMUTE', VUnPermute, error_code, &
                                .true., MLabel,  'INFORMATION="VERTEX"', myScope, 'NODES')
    end if

    if (error_code == 0 ) then
       ErrString = 'Error writing Vertex to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'COORDS', VCoords, error_code, &
                                .true., MLabel,  'INFORMATION="VERTEX"', myScope, 'NODES')
    end if

    if (error_code == 0 .and. present(VPartition)) then
       ErrString = 'Error writing VertexPartition to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'PARTITIONS', VPartition, error_code, &
                                .true., MLabel,  'INFORMATION="VERTEX"', myScope, 'NODES')
    end if
    if (error_code == 0 .and. present(VBoundary)) then
       ErrString = 'Error writing VertexBoundary to xml file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'BOUNDARY', VBoundary, error_code, &
                                .true., MLabel,  'INFORMATION="VERTEX"', myScope, 'NODES')
    end if

    if (error_code == 0 .and. present(VRsum_RVol)) then
       ErrString = 'Error Writing VRsum_RVol to mesh file'
       call Tbu_make_file_entry(BaseBrook, BAside, 'RSUMRVOL', VRsum_RVol, error_code, &
                                .true., MLabel,  'INFORMATION="VERTEX"', myScope, 'NODES')
    end if

    if ( error_code == 0 ) then
       ErrString = 'Error Closing mesh tag in mesh file'
       call tbrook_CloseXMLTag(BaseBrook, XMLTag="MESH", SCOPE=myScope, iStatus=error_code)
    end if


    ! Close and destroy our brooks
    if ( error_code == 0 ) then
       ErrString = 'Error closing Mesh File'
       call Tbrook_Close(BaseBrook, iStatus=error_code)
    end if

    if ( error_code == 0 ) then
       ErrString = 'Error destroying Mesh File'
       call TBrook_Destroy(BaseBrook, error_code)
    end if

    if ( error_code == 0 ) then
       ErrString = 'Error closing Mesh Data File'
       call Tbrook_Close(BAside, iStatus=error_code)
    end if

    if ( error_code == 0 ) then
       ErrString = 'Error destroying Mesh Data File'
       call TBrook_Destroy(BAside, error_code)
    end if

    ! Write the MESH entry in the main xml file
    if ( p_info%IOP  .or. myScope == TB_SCOPE_LOCAL) then 
 
       ! Open Mesh Tag 
       if ( error_code == 0 ) then
          ErrString = 'Error opening mesh tag in xml file'
          call tbrook_OpenXMLTag(B, XMLTag='MESH', &
                                 XMLAttributes='FORMAT="TRUCHASMESH" LABEL="'//trim(MLabel)// &
                                               '" MESHWRITERVERSION="'//trim(tbu_version)//'"', &
                                 SCOPE=myScope, iStatus=error_code)
       end if
       
       ! Write Look Aside Tag 
       if ( error_code == 0 ) then
          ErrString = 'Error writing mesh file name to xml file'
          call tbrook_WriteXMLTAG(B,                                & 
                                  XMLTag='FILE',                    & 
                                  XMLAttributes='FORMAT="xml"',     & 
                                  XMLStringData=TRIM(MeshFile),     & 
                                  SCOPE=myScope,                    & 
                                  iStatus=error_code) 
       end if
 
       ! Close tag in main file 
       if ( error_code == 0 ) then
          ErrString = 'Error Closing mesh tag in xml file'
          call tbrook_CloseXMLTag(B, XMLTag='MESH', SCOPE=myScope, iStatus=error_code) 
       end if
    end if 

    if ( error_code /= 0 ) then
       call TBU_Abort(Message=ErrString,                &
            SUBROUTINENAME='TBU_MeshWriter',  &
            error_code=error_code,            &
            iStatus=iStatus)
    end if

  END SUBROUTINE TBU_MeshWriter


  SUBROUTINE TBU_TIMESTEP_OUTPUT (iStatus)
    !---------------------------------------------------------------------------
    ! Purpose:
    !   Write out a Truchas timestep to the *.TBrook.xml file (BaseBrook)
    !---------------------------------------------------------------------------
    use gap_output, only: set_gap_element_output
    use file_utility,         only: Make_File_Name
    use output_utilities,     only: ANNOUNCE_FILE_WRITE
    use output_data_module,   only: XML_Data_Format,     &
                                    XML_IFormat
    use time_step_module,     only: t, dt, cycle_number, dt_courant
    use tbrook_module,        only: Brook,               &
                                    BaseBrook,           &
                                    TBrook_WriteXMLTag,  &
                                    TBrook_OpenXMLTag,   & 
                                    TBrook_CloseXMLTag,  & 
                                    TBrook_Set,          & 
                                    ASSIGNMENT(=),       & 
                                    B_IFORM_BINARY,      & 
                                    TBrook_Destroy,      & 
                                    TBrook_Close 
#ifdef DEBUG_TBU 
#define UTILIZE use 
    UTILIZE Tbrook_module,     only: Brook_Unit,         & 
                                     Brook_File,         & 
                                     Brook_Form 
#undef UTILIZE 
#endif 
 
    ! Argument List 
    integer, optional :: iStatus 
    ! Local Variables 
    type(Brook), target :: B_Tmp 

    character(LEN=256) :: tmpString 
    character(LEN=256) :: filename 
    integer, save :: file_count = 0 
    integer :: error_code 
 
    !--------------------------------------------------------------------------- 
 
    ! Set gap element data, if necessary (provides non-spurious field values for cells near gaps) 
    call SET_GAP_ELEMENT_OUTPUT() 
 
    error_code = 0 
 
    ! Create Timestep Tag Attribute before incrementing file_count 
637 FORMAT('CYCLE="',i20,'" t="',1e30.20,'" dt="',1e30.20,'" dt_courant="',1e30.20,'" ID="',i10,'"') 
    write(tmpString,637) cycle_number, t, dt, dt_courant, file_count 
 
    ! Create lookaside binary or xml file *.00*.ts.bin or *.00*.ts.xml (B_Tmp)
    if (XML_IFormat == B_IFORM_BINARY) then 
       filename = MAKE_FILE_NAME (string ='ts', number = file_count, suffix = 'bin') 
    else 
       filename = MAKE_FILE_NAME (string ='ts', number = file_count, suffix = 'xml') 
    end if 

    if ( error_code == 0 ) then 
       call Tbrook_Set(B  = B_Tmp,                  &
            file          = TRIM(filename),         &
            form          = TRIM(XML_data_format),  &
            iStatus       = error_code) 
    end if

    file_count = file_count + 1 
 
    ! Open Timestep Tag 
    if ( error_code == 0 ) then 
       call Tbrook_OpenXMLTag(B             = BaseBrook,  &
                              XMLTag        = "TIMESTEP", & 
                              XMLAttributes = tmpString,  &
                              iStatus       = error_code) 
    end if
 
    ! Write Look Aside Tag 
    if ( error_code == 0 ) then 
       call Tbrook_WriteXMLTag(B             = BaseBrook,                                       & 
                               XMLTag        = "FILE",                                          & 
                               XMLAttributes = 'FORMAT="'//trim(ADJUSTL(XML_Data_Format))//'"', & 
                               XMLStringData = TRIM(filename),                                  & 
                               iStatus       = error_code) 

    end if

    ! We will add variable information to the file as needed 
 
#ifdef DEBUG_TBU 
    write(*,*) 'B_Tmp diagnostic error: ',  error_code, & 
               '__Unit: ',Brook_Unit(B_Tmp),        & 
               '__File:  ',TRIM(Brook_File(B_Tmp)), & 
               '__Form:  ',TRIM(Brook_Form(B_Tmp)) 
    write(*,*) 'Compare file: ', TRIM(filename) 
#endif 

    ! Write timestep tag to temporary brook 
    if ( error_code == 0 ) then 
       call Tbrook_OpenXMLTag(B             = B_Tmp,      &
                              XMLTag        = "TIMESTEP", &
                              XMLAttributes = tmpString,  &
                              iStatus       = error_code) 
    end if

    call TBU_WriteTimeStepData(BBase        = BaseBrook,  &
                               BAside       = B_Tmp,      &
                               iStatus      = error_code)
 
    ! Close tag in main file 
    if ( error_code == 0 ) then 
       call Tbrook_CloseXMLTag(B            = BaseBrook,  &
                               XMLTag       = "TIMESTEP", &
                               iStatus      = error_code) 
    end if
 
#ifdef DEBUG_TBU 
    write(*,*) 'B_Tmp diagnostic error: ',  error_code, & 
               '__Unit: ', Brook_Unit(B_Tmp),           & 
               '__File:  ',TRIM(Brook_File(B_Tmp)),     & 
               '__Form:  ',TRIM(Brook_Form(B_Tmp)) 
#endif 
    call Tbrook_CloseXMLTag(B = B_Tmp,                &
         XMLTag               = "TIMESTEP",           & 
         iStatus              = error_code) 

    call TBrook_Close(B_Tmp,   & 
           iStatus = error_code) 

    call TBrook_Destroy(B_Tmp, &
         iStatus = error_code ) 

    call ANNOUNCE_FILE_WRITE('Timestep dump', filename) 

    if ( present(iStatus)) then 
       iStatus = error_code 
    end if 

  END SUBROUTINE TBU_TIMESTEP_OUTPUT 
 

  SUBROUTINE TBU_LASTSTEP_OUTPUT (iStatus)
    !--------------------------------------------------------------------------- 
    ! Purpose: 
    ! 
    !  Dumps very last successful timestep data, to a lookaside xml file
    !  *.LastStep.TBrook.xml (BLastStep) in order to restart and diagnose
    !  from this final successful timestep
    !
    !--------------------------------------------------------------------------- 

    use file_utility,      only: MAKE_FILE_NAME
    use time_step_module,  only: t, dt, cycle_number, dt_courant
    use tbrook_module,     only: Brook,    &
         BLastStep,                        &
         TBrook_WriteXMLTag,               &
         TBrook_OpenXMLTag,                & 
         TBrook_CloseXMLTag,               & 
         TBrook_Set,                       & 
         ASSIGNMENT(=),                    & 
         TBrook_Destroy,                   & 
         TBrook_Close 

    ! Argument List
    integer, optional :: iStatus 

    ! Local Variables
    character(LEN=256) :: tmpString, filename, LastStepFile 
    type(Brook), target :: BAside
    integer :: error_code, file_count

    error_code = 0
    file_count = 0

    ! Create Timestep Tag Attribute before incrementing file_count 
638 FORMAT('CYCLE="',i20,'" t="',1e30.20,'" dt="',1e30.20,'" dt_courant="',1e30.20,'" ID="',i10,'"') 
    write(tmpString,638) cycle_number, t, dt, dt_courant, file_count 

    if (error_code == 0) then 

       !(Re)create lookaside xml file *.LastStep.xml (=BLastStep)
       LastStepFile = MAKE_FILE_NAME (string=TRIM('LastStep'), suffix='xml') 

       call Tbrook_Set(B = BLastStep,          &
            File         = TRIM(LastStepFile), & 
            Form         = "xml",              &
            iStatus      = error_code) 
    end if

    if (error_code == 0) then 

       ! Create lookaside binary file *.lts.bin (=BAside)
       filename = MAKE_FILE_NAME (string='lts', suffix='bin') 

       call Tbrook_Set(B = BAside,           & 
            file         = TRIM(filename),   &
            Form         = "binary",         &
            iStatus      = error_code) 
 
    end if
 
    ! Open Timestep Tag 
    if ( error_code == 0 ) then 
       call Tbrook_OpenXMLTag(B = BLastStep,  &
            XMLTag              = "LASTSTEP", &
            XMLAttributes       = tmpString,  &
            iStatus             = error_code) 
 
    end if

    ! Write Look Aside Tag 
    if ( error_code == 0 ) then 
       call Tbrook_WriteXMLTag(B = BLastStep,            & 
            XMLTag               = "FILE",               & 
            XMLAttributes        = 'FORMAT="binary"',    & 
            XMLStringData        = TRIM(filename),       & 
            iStatus              = error_code) 
    end if
       
    call TBU_WriteTimeStepData(BBase = BLastStep, &
         BAside                      = BAside,    &
         iStatus                     = error_code)

       
    ! close TimeStep tag
    if ( error_code == 0 ) then 
       call Tbrook_CloseXMLTag(B = BLastStep,  &
            XMLTag               = "LASTSTEP", &
            iStatus              = error_code) 
    end if

    call TBrook_Close(BAside,   & 
         iStatus = error_code) 

    call TBrook_Destroy(BAside, &
         iStatus = error_code ) 

    call TBrook_Close(BLastStep,   & 
         iStatus = error_code) 

    call TBrook_Destroy(BlastStep, &
         iStatus = error_code ) 

    if ( present(iStatus)) then 
       iStatus = error_code 
    end if 

  END SUBROUTINE TBU_LASTSTEP_OUTPUT

  SUBROUTINE TBU_WriteTimeStepData (BBase, BAside, iStatus)
    !--------------------------------------------------------------------------- 
    ! Purpose: 
    ! 
    !   for a given base XML file (BBase, e.g *.TBrook.xml) and 
    !   for a given binary lookaside file (BAside, e.g *.ts.*.bin)
    !   write out timestep data
    !--------------------------------------------------------------------------- 

    use fluid_data_module,    only: fluid_flow, boussinesq_approximation, courant
    use matl_module,          only: GATHER_VOF, Matl
    use parameter_module,     only: ncells, ndim, nmat, mat_slot
    use diagnostics_module,   only: DIVERGENCE
    use property_module,      only: get_density, get_density_delta, get_user_material_id
    use zone_module,          only: Zone
    use node_operator_module, only: nipc
    use solid_mechanics_data, only: solid_mechanics, &
                                    SMech_Cell,      &
                                    SMECH_IP,        &
                                    RHS,             &
                                    Displacement,    &
                                    Thermal_Strain,  &
                                    PC_Strain,       &
                                    Node_Gap,        &
                                    Node_Norm_Trac,  &
                                    Rotation_Magnitude
    use mech_bc_data_module,  only: Interface_List
    use pgslib_module,         only: PGSLib_GLOBAL_ANY
    use EM_data_proxy,        only: joule_power_density, EM_is_on

    use tbrook_module,        only: Brook,           &
                                    ASSIGNMENT(=),   & 
                                    TB_SCOPE_GLOBAL
    use time_step_module,          only: dt
    use diffusion_solver_data,     only: ds_enabled, num_species, heat_eqn
    use diffusion_solver,          only: ds_get_phi, ds_get_temp_grad
    use string_utilities,          only: i_to_c
    use gap_output, only: set_gap_element_output
    use physics_module, only: heat_transport, heat_species_transport
 
    ! Argument List 
    integer, optional   :: iStatus 
    type(Brook), target :: BBase, BAside 

    ! Local Variables 
    integer :: n
    real(r8), dimension(ncells) :: Vof, fluidDeltaRho, materialRho 
    real(r8), pointer, dimension(:,:) :: tmp_2dArray 
    character(LEN=32) :: vof_name
    real(r8), pointer, dimension(:,:) :: gradT

    tmp_2dArray => NULL()

    !Allocate temporary array
    Allocate(tmp_2dArray(ndim, ncells))
    
    ! Write variable data. 
    if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside, "MATSLOT", mat_slot, iStatus) 

    ! If gap elements are present, the data in them is bogus and needs to
    ! be set to the data for an adjacent cell.
    call SET_GAP_ELEMENT_OUTPUT ()
 
    ! Write out derived types describing the thermodynamic state of the problem. 
    ! Note use of dummy variable in the file entry 
    if ( iStatus == 0 ) call tbu_zone_out(BAside, Zone, iStatus, TB_SCOPE_GLOBAL, BBase) 

     ! Density - calculate from current VOF data and material density instead of Zone%Rho or 
     ! phase densities.  Dave Korzekwa 02-04-09 
     ! This uses Vof as a temporary ncells array
     call get_density (Zone%Temp, materialRho)
     call tbu_make_file_entry(BBase, BAside,"Z_RHO", materialRho, iStatus) 

    if ( iStatus == 0 ) call tbu_flow_out(BAside, iStatus, TB_SCOPE_GLOBAL, BBase) 

    if ( iStatus == 0 ) call tbu_matl_out(BAside, Matl, iStatus, TB_SCOPE_GLOBAL, BBase) 
 
    ! Variables relevant to electromagnetic induction heating. 
    if (EM_is_on()) then 
       ! Joule power density 
       Vof = joule_power_density() 
       if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside,"Joule_P", Vof, iStatus) 
    end if 
    ! Variables relevant to heat transfer. 
    if (heat_transport .or. heat_species_transport) then 
       ! Instantaneous temperature change. 
       vof = (Zone%Temp - Zone%Temp_Old)/dt 
       if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside,"dTdt", Vof, iStatus) 
    end if 
 
    if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside, "COURANT", courant, iStatus)

    ! Material volume fraction(s). 
    if (nmat > 1) then 
       do n = 1, nmat 
          call GATHER_VOF (n,Vof) 
3         format('VOF',i4.4) 
          write (vof_name,3) get_user_material_id(n) 
          if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside, vof_name, Vof, iStatus) 
       end do 
    end if 
 
    ! Delta Density 
    if ( (iStatus == 0 ) .and. fluid_flow .and. boussinesq_approximation) then 
      call get_density_delta (Zone%Temp, fluidDeltaRho)
      call tbu_make_file_entry(BBase, BAside,"del-rho", fluidDeltaRho, iStatus) 
    endif 
 
    ! divergence (Vof) (only if fluid flow is on). 
    if (iStatus == 0 .and. fluid_flow) then 
       call DIVERGENCE(Vof) 
       call tbu_make_file_entry(BBase, BAside,"Volume_Error", Vof, iStatus) 
    end if 
 
    ! Stresses, strains and displacements 
    if (iStatus == 0 .and. solid_mechanics) then 
 
       ! Cell centroid - useful for postprocessing stresses and strains 
11     FORMAT(a,i2.2) 
       do n = 1, nipc 
          write(vof_name,11) 'TOTAL_STRAIN_',n 
          if (iStatus == 0 ) call tbu_make_file_entry(BBase, BAside, vof_name,           & 
                                                         SMech_IP(n)%Total_Strain, iStatus) 
          write(vof_name,11) 'ELASTIC_STRESS_',n 
          if (iStatus == 0 ) call tbu_make_file_entry(BBase, BAside, vof_name,         & 
                                                         SMech_IP(n)%Elastic_Stress, iStatus) 
          write(vof_name,11) 'PLASTIC_STRAIN_',n 
          if (iStatus == 0 ) call tbu_make_file_entry(BBase, BAside, vof_name,         & 
                                                         SMech_IP(n)%Plastic_Strain, iStatus) 
          write(vof_name,11) 'PLASTIC_STRAIN_RATE_',n 
          if (iStatus == 0 ) call tbu_make_file_entry(BBase, BAside, vof_name,    & 
                                                         SMech_IP(n)%Plastic_Strain_Rate, iStatus) 
       end do 
 
       if (iStatus == 0 ) call tbu_make_file_entry(BBase, BAside, 'RHS', RHS, iStatus, MAP="NODE") 
 
       if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside,'sigma',     SMech_Cell%Elastic_Stress, iStatus) 
       if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside,'epsilon',   SMech_Cell%Total_Strain,   iStatus) 
       if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside,'e_plastic', SMech_Cell%Plastic_Strain, iStatus) 
 
       if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside,'epstherm', Thermal_Strain, iStatus) 
       if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside,'epspc', PC_Strain, iStatus) 
 
       Vof = SMech_Cell%Plastic_Strain_Rate(:) 
       if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside,"epsdot",       Vof,          iStatus) 
       if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside,"Rotation", Rotation_Magnitude, iStatus) 
       if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside,"Displacement", Displacement, iStatus, MAP="NODE") 

       ! Gap Displacements and scaled forces
       do n = 1, SIZE(Node_Gap,2)
          if (PGSLib_GLOBAL_ANY((Node_Gap(:,n) /= 0.0).or.(Node_Norm_Trac(:,n) /= 0.0))) then
             write(vof_name,11) 'GAP_',Interface_List(n) 
             if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside,vof_name, &
                  Node_Gap(:,n), iStatus, MAP="NODE") 
             write(vof_name,11) 'NTRAC_',Interface_List(n) 
             if ( iStatus == 0 ) call tbu_make_file_entry(BBase, BAside,vof_name, &
                  Node_Norm_Trac(:,n), iStatus, MAP="NODE") 
          end if
       end do

    end if

    ! Diffusion solver solution fields not otherwise known.
    if (ds_enabled) then
      do n = 1, num_species
        call ds_get_phi (n, vof)  ! use VOF as a ncells-length temporary
        if (iStatus == 0) call tbu_make_file_entry (BBase, BAside, "phi"//i_to_c(n), vof, iStatus)
      end do

      if (heat_eqn) then
         allocate(gradT(ndim,ncells))
         call ds_get_temp_grad(gradT)
         if (iStatus == 0) call tbu_make_file_entry (BBase, BAside, "Grad_T", gradT, iStatus)
         
         vof(:) = 0.0_r8
         do n=1,ndim
            vof(:) = vof(:) + gradT(n,:)**2
         enddo
         vof(:) = sqrt(vof(:))

         if (iStatus == 0) call tbu_make_file_entry (BBase, BAside, "||Grad_T||", vof, iStatus)
         deallocate(gradT)
      endif
    end if

    ! Destroy work arrays. 
    if (ASSOCIATED(Tmp_2dArray)) Deallocate(Tmp_2dArray) 

  END SUBROUTINE TBU_WriteTimeStepData

 
  subroutine tbu_zone_out(b, Variable, istatus, scope, b_base) 
    !================================================================== 
    ! Purpose(s): 
    !   Output zone to stream 
    !================================================================== 
    use tbrook_module, only: brook, assignment(=) 
    use zone_module, only: CELL_AVG 
    use parameter_module, only: ndim

    type(brook), intent(in) :: b 
    type(brook), optional, target, intent(in) :: b_base 
    type(CELL_AVG),  target, dimension(:), intent(IN) :: Variable 
    integer, intent(in), optional :: scope 
    integer, intent(inout) :: istatus 
 
    ! Arguments 
    ! Local variables 
    type(CELL_AVG), pointer, dimension(:) :: Zone 
    real(r8), dimension(:,:), allocatable :: tmp_2DArray 
    integer :: n, out_stat 
    type(Brook), pointer :: B_lbase 
    type(Brook), target  :: B_dummy 
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
 
    if(present(B_base)) then 
       b_lbase => B_base 
    else 
       b_lbase => b_dummy 
    end if 
 
    out_stat = 0 
    Zone => Variable 
 
! Write out density calculated from current volume fractions and material densities (not phase namelist densities)
!  instead of these.  Remove old density completely from the output to avoid confusion.  Dave Korzekwa 2-04-09
!    call Tbu_make_file_entry(B_lbase, B, 'Z_RHO',         Zone(:)%Rho,          out_stat, .true., 'DefaultMesh', ' ', Scope) 
!    call Tbu_make_file_entry(B_lbase, B, 'Z_RHOOLD',      Zone(:)%Rho_Old,      out_stat, .true., 'DefaultMesh', ' ', Scope) 
    call Tbu_make_file_entry(B_lbase, B, 'Z_TEMP',        Zone(:)%Temp,         out_stat, .true., 'DefaultMesh', ' ', Scope) 
    call Tbu_make_file_entry(B_lbase, B, 'Z_TEMP_OLD',    Zone(:)%Temp_Old,     out_stat, .true., 'DefaultMesh', ' ', Scope) 
    call Tbu_make_file_entry(B_lbase, B, 'Z_ENTHALPY',    Zone(:)%Enthalpy,     out_stat, .true., 'DefaultMesh', ' ', Scope) 
    call Tbu_make_file_entry(B_lbase, B, 'Z_ENTHALPYOLD', Zone(:)%Enthalpy_Old, out_stat, .true., 'DefaultMesh', ' ', Scope) 
    call Tbu_make_file_entry(B_lbase, B, 'Z_P',           Zone(:)%P,            out_stat, .true., 'DefaultMesh', ' ', Scope) 
 
 
    allocate(tmp_2DArray(ndim, size(zone,1))) 
    do n = 1,ndim 
       tmp_2DArray(n, :)  = Zone(:)%Vc(n) 
    end do 
    call Tbu_make_file_entry(B_lbase, B, 'Z_VC', tmp_2DArray, out_stat, .true., 'DefaultMesh', 'Attr="none"', Scope) 
 
    do n = 1,ndim 
       tmp_2DArray(n, :)  = Zone(:)%Vc_Old(n) 
    end do 
    call Tbu_make_file_entry(B_lbase, B, 'Z_VCOLD', tmp_2DArray, out_stat, .true., 'DefaultMesh', 'Attr="none"' , Scope) 
    deallocate(tmp_2DArray) 
    
    istatus = out_stat
 
  END subroutine tbu_zone_out 
 
  subroutine tbu_flow_out(b, istatus, scope, b_base) 
    !================================================================== 
    ! Purpose(s): 
    !   Output zone to stream 
    !================================================================== 
    use tbrook_module, only: brook, assignment(=) 
    use fluid_data_module, only: Fluxing_Velocity

    type(brook), intent(in) :: b 
    type(brook), optional, target, intent(in) :: b_base 
    integer, intent(in), optional :: scope 
    integer, intent(inout) :: istatus 
 
    ! Arguments 
    ! Local variables 
    integer :: out_stat 
    type(Brook), pointer :: B_lbase 
    type(Brook), target  :: B_dummy 
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
 
    if(present(B_base)) then 
       b_lbase => B_base 
    else 
       b_lbase => b_dummy 
    end if 
 
    out_stat = 0 
 
    call Tbu_make_file_entry(B_lbase, B, 'Face_Vel', Fluxing_Velocity, out_stat, .true., 'DefaultMesh', 'Attr="none"' , Scope) 
 
    istatus = out_stat

  END subroutine tbu_flow_out 
 
 
  subroutine tbu_MATL_out(b, Variable, istatus, scope, b_base) 
    !================================================================== 
    ! Purpose(s): 
    !   Output zone to stream 
    !================================================================== 
    use matl_module,   only: MATL_SLOT 
    use tbrook_module, only: brook, assignment(=) 
    use parameter_module, only: mat_slot 
 
    type(brook), intent(in) :: b 
    type(brook), optional, target, intent(in) :: b_base 
    type(MATL_SLOT), target, dimension(:), intent(IN) :: Variable 
    integer, intent(in), optional :: scope 
    integer, intent(inout) :: istatus 
 
    ! Arguments 
    ! Local variables 
    integer :: s, out_stat 
    type(MATL_SLOT), pointer, dimension(:) :: Matl 
    type(Brook), pointer :: B_lbase 
    type(Brook), target  :: B_dummy 
    character(LEN=256) :: myAttr 
 
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> 
 
    if(present(B_base)) then 
       b_lbase => B_base 
    else 
       b_lbase => b_dummy 
    end if 
 
    out_stat = 0 
    Matl => Variable 
 
    ! The assumption is that mat_slot is common across all Processors 
1   Format('SLOT="',I10,'"') 
    do s = 1,mat_slot 
      write(myAttr,1) s 
      call Tbu_make_file_entry(B_lbase, B, 'M_ID',      Matl(s)%Cell(:)%ID,      out_stat, .true., 'DefaultMesh', myAttr, Scope) 
      call Tbu_make_file_entry(B_lbase, B, 'M_VOF',     Matl(s)%Cell(:)%Vof,     out_stat, .true., 'DefaultMesh', myAttr, Scope) 
      call Tbu_make_file_entry(B_lbase, B, 'M_VOF_OLD', Matl(s)%Cell(:)%Vof_OLD, out_stat, .true., 'DefaultMesh', myAttr, Scope) 
    end do 
 
    iStatus=out_stat 

  END subroutine tbu_MATL_out 
 
end module tbrook_utility 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
!! in the LICENSE file found in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BODY_INPUT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures for the input of body parameters.
  !
  ! Contains: BODY_CHECK
  !           BODY_DEFAULT
  !           BODY_INPUT
  !           BODY_INPUT_PARALLEL
  !           TAB_DATA_INPUT
  !
  ! Author(s): Jerry S. Brock, LANL (jsbrock@lanl.gov)
  !            Douglas B. Kothe, LANL (dbk@lanl.gov)
  !            Larry J. Cox, LANL (ljcox@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  public :: BODY_CHECK, BODY_DEFAULT, BODY_INPUT

CONTAINS

  SUBROUTINE BODY_CHECK (fatal)
    !=======================================================================
    ! Purpose(s):
    !
    !   Check BODY namelist.
    !
    !=======================================================================
    use input_utilities,   only: NULL_I, NULL_R
    use interfaces_module, only: Axis,                                 &
                                 File_Name, Fill, Height, Isrftype,    &
                                 Matnum, nbody, Nsurf, Offset,         &
                                 Radius, Rotangl, Rotpt, Rtab,         &
                                 Rtheta_Tabular_Pt, RZ_Tabular_Pt,     &
                                 Sgeom, Surface_Name, Tab_Type,        &
                                 Length, Ztab, Mesh_Matnum
    use parameter_module,  only: msurf, mtab, mtype, nmat
    use legacy_mesh_api,   only: ndim
    use property_module,   only: Get_User_Material_ID

    ! Argument List
    logical :: fatal

    ! Local Variables
    integer :: is, j, l, m, n
    integer, parameter :: maxforms = 16
    real(r8) :: degree, rtheta_flag

    character(32), dimension(maxforms) :: Axis_Form,    &
                                          Fill_Inside,  &
                                          Fill_Outside, &
                                          Tab_Type_Form
    character(32), dimension(maxforms,0:mtype) :: Srf_Forms

    character, dimension(ndim) :: Axis_Label
    character(128) :: message

    ! Axis labels
    data Axis_Label / 'x', 'y', 'z' /

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Initialize character string forms.
    Fill_Inside = ''
    Fill_Inside(1) = 'fill inside'
    Fill_Inside(2) = 'inside'
    Fill_Inside(3) = '+'
    Fill_Inside(4) = 'right'
    Fill_Inside(5) = 'above'
    Fill_Inside(6) = 'top'
    Fill_Inside(7) = 'front'

    Fill_Outside    = ''
    Fill_Outside(1) = 'fill outside'
    Fill_Outside(2) = 'outside'
    Fill_Outside(3) = '-'
    Fill_Outside(4) = 'left'
    Fill_Outside(5) = 'below'
    Fill_Outside(6) = 'bottom'
    Fill_Outside(7) = 'back'

    Axis_Form    = ''
    Axis_Form(1) = 'x'
    Axis_Form(2) = 'x-axis'
    Axis_Form(3) = '1'
    Axis_Form(4) = 'y'
    Axis_Form(5) = 'y-axis'
    Axis_Form(6) = '2'
    Axis_Form(7) = 'z'
    Axis_Form(8) = 'z-axis'
    Axis_Form(9) = '3'

    Tab_Type_Form    = ''
    Tab_Type_Form(1) = 'rotation'
    Tab_Type_Form(2) = 'rotate'
    Tab_Type_Form(3) = 'rot'
    Tab_Type_Form(4) = 'translation'
    Tab_Type_Form(5) = 'translate'
    Tab_Type_Form(6) = 'trans'

    Srf_Forms      = ''
    ! Background
    Srf_Forms(1,0) = 'fill'
    Srf_Forms(2,0) = 'background'
    Srf_Forms(3,0) = 'back'
    Srf_Forms(4,0) = 'bck'
    ! (2) Planar surface
    Srf_Forms(1,1) = 'pla'
    Srf_Forms(2,1) = 'planar'
    Srf_Forms(3,1) = 'plane'
    ! (3) Box
    Srf_Forms(1,2) = 'box'
    Srf_Forms(2,2) = 'parallelopiped'
    Srf_Forms(3,2) = 'cube'
    ! (4) Sphere
    Srf_Forms(1,3) = 'sph'
    Srf_Forms(2,3) = 'sphere'
    ! (5) Ellipsoid
    Srf_Forms(1,4) = 'ell'
    Srf_Forms(2,4) = 'ellipse'
    Srf_Forms(3,4) = 'ellipsoid'
    ! (6) Cylinder
    Srf_Forms(1,5) = 'cyl'
    Srf_Forms(2,5) = 'cylinder'
    ! (7) Cone
    Srf_Forms(1,6) = 'cone'
    Srf_Forms(2,6) = 'con'
    Srf_Forms(3,6) = 'conical'
    ! (8) Tabular
    Srf_Forms(1,7) = 'tab'
    Srf_Forms(2,7) = 'tabular'
    ! (9) Read in from the mesh file.
    Srf_Forms(1,8) = 'mesh'
    Srf_Forms(2,8) = 'mesh file'
    Srf_Forms(3,8) = 'from mesh file'
    Srf_Forms(4,8) = 'get from mesh file'
    Srf_Forms(5,8) = 'in mesh file'

    ! Initialize variables
    degree = ACOS(-1.0_r8)/180.0
    fatal = .false.

    ! Check material number
    m = Matnum(nbody)
    if (m < 1 .or. m > nmat) then
       write (message, 10) nbody, m, nmat
10     format('Body ',i2,' material number = ',i2,' < 1 or > nmat = ',i2,'!')
       call TLS_error (message)
       fatal = .true.
       goto 300
    end if
    if (m == NULL_I) then
       write (message, 15) nbody
15     format('Body ',i2,' material number not defined!')
       call TLS_error (message)
       fatal = .true.
       goto 300
    end if

    ! Identify surface type
    SURFACE: do m = 1, msurf
       if (Surface_Name(m) /= 'none') then
          do n=0,mtype
             do l=1,maxforms
                if (Surface_Name(m) == Srf_Forms(l,n)) go to 45
             end do
          end do

          write (message, 40) m, nbody, TRIM(Surface_Name(m))
40        format('Surface ',i2,' of body ',i2,' has an unknown surface type ', a)
          call TLS_error (message)
          fatal = .true.
          cycle SURFACE

45        continue               ! Found surface type
          Isrftype(m,nbody) = n  ! Save surface type id

          ! Fill on the "outside" of the surface?
          do l=1,maxforms
             if (Fill(m) == Fill_Outside(l)) Isrftype(m,nbody) = -Isrftype(m,nbody)
          end do

          if (Isrftype(m,nbody) == 0) then

             ! Background material
             Nsurf(nbody) = 0

          else if (ABS(Isrftype(m,nbody)) == 1) then

             ! Plane
             do l=1,maxforms
                if (Axis(m)(1:3) == Axis_Form(l)(1:3)) goto 55
             end do

             write (message, 50) m,nbody
50           format('No axis specified for plane surface ',i3, ' of body ',i3)
             call TLS_error (message)
             fatal = .true.
             cycle SURFACE

55           continue
             if (l >= 1 .and. l <= 3) then
                Sgeom(1,m,nbody) = 1
             else if (l >= 4 .and. l <= 6) then
                Sgeom(1,m,nbody) = 2
             else if (l >= 7 .and. l <= 9) then
                Sgeom(1,m,nbody) = 3
             end if

          else if (ABS(Isrftype(m,nbody)) == 2) then

             ! Box
             do n = 1,3
                Sgeom(n,m,nbody) = Length(n,m)
                if (Length(n,m) == NULL_R) then
                   write (message,60) Axis_label(n),m,nbody
60                 format('No ',a1,'-Length specified for box surface ',i3,' of body ',i3)
                   call TLS_error (message)
                   fatal = .true.
                end if
             end do

          else if (ABS(Isrftype(m,nbody)) == 3) then

             ! Sphere
             Sgeom(1,m,nbody) = ABS(Radius(1,m))
             if (Radius(1,m) == NULL_R) then
                write (message,65) m,nbody
65              format('No radius specified for sphere surface ',i3,' of body ',i3)
                call TLS_error (message)
                fatal = .true.
             end if

          else if (ABS(Isrftype(m,nbody)) == 4) then

             ! Ellipsoid
             do n = 1,3
                Sgeom(n,m,nbody) = ABS(Radius(n,m))
                if (Radius(n,m) == NULL_R) then
                   write (message,70) n,m,nbody
70                 format('No radius(',i1,') specified for ellipsoid surface ',i3,' of body ',i3)
                   call TLS_error (message)
                   fatal = .true.
                end if
             end do

          else if (ABS(Isrftype(m,nbody)) == 5) then

             ! Cylinder
             do l = 1,maxforms
                if (Axis(m)(1:3) == Axis_Form(l)(1:3)) goto 80
             end do

             write (message,75) m,nbody
75           format('No axis specified for cylinder surface ',i3,' of body ',i3)
             call TLS_error (message)
             fatal = .true.
             cycle SURFACE

80           continue
             if (l >= 1 .and. l <= 3) then
                Sgeom(1,m,nbody) = 1
             else if (l >= 4 .and. l <= 6) then
                Sgeom(1,m,nbody) = 2
             else if (l >= 7 .and. l <= 9) then
                Sgeom(1,m,nbody) = 3
             end if

             Sgeom(2,m,nbody) = ABS(Radius(1,m))
             if (Radius(1,m) == NULL_R) then
                write (message,85) m,nbody
85              format('No radius specified for cylinder surface ',i3,' of body ',i3)
                call TLS_error (message)
                fatal = .true.
             end if

             Sgeom(3,m,nbody) = Height(m)
             if (Height(m) == NULL_R) then
                write (message,90) m,nbody
90              format('No height specified for cylinder surface ',i3,' of body ',i3)
                call TLS_error (message)
                fatal = .true.
             end if

          else if (ABS(Isrftype(m,nbody)) == 6) then

             ! Cone
             do l = 1,maxforms
                if (Axis(m)(1:3) == Axis_Form(l)(1:3)) goto 100
             end do

             write (message,95) m,nbody
95           format('No axis specified for cone surface ',i3,' of body ',i3)
             call TLS_error (message)
             fatal = .true.
             cycle SURFACE

100          continue
             if (l >= 1 .and. l <= 3) then
                Sgeom(1,m,nbody) = 1
             else if (l >= 4 .and. l <= 6) then
                Sgeom(1,m,nbody) = 2
             else if (l >= 7 .and. l <= 9) then
                Sgeom(1,m,nbody) = 3
             end if

             Sgeom(2,m,nbody) = Height(m)
             if (Height(m) == NULL_R) then
                write (message,105) m,nbody
105             format('No height specified for cone surface ',i3,' of body ',i3)
                call TLS_error (message)
                fatal = .true.
             end if

             Sgeom(3,m,nbody) = ABS(Radius(1,m))
             if (Radius(1,m) == NULL_R) then
                write (message,110) m,nbody
110             format('No radius specified for cone surface ',i3,' of body ',i3)
                call TLS_error (message)
                fatal = .true.
             end if

          else if (ABS(Isrftype(m,nbody)) == 7) then

             ! Tabular surface

             ! Tabular type
             do l=1,6
                if (Tab_Type(m)(1:3) == Tab_Type_Form(l)(1:3)) goto 120
             end do

             write (message,115) m,nbody
115          format('No Tab_Type specified for tabular surface ',i3,' of body ',i3)
             call TLS_error (message)
             fatal = .true.
             cycle SURFACE

120          continue
             if (l <= 6 .and. l >= 4) Sgeom(1,m,nbody) = 1

             ! Axis
             do l = 1,maxforms
                if (Axis(m)(1:3) == Axis_Form(l)(1:3)) goto 130
             end do

             write (message,125) m,nbody
125          format('No axis specified for tabular surface ',i3,' of body ',i3)
             call TLS_error (message)
             fatal = .true.
             cycle SURFACE

130          continue
             if (l >= 1 .and. l <= 3) then
                Sgeom(2,m,nbody) = 1
             else if (l >= 4 .and. l <= 6) then
                Sgeom(2,m,nbody) = 2
             else if (l >= 7 .and. l <= 9) then
                Sgeom(2,m,nbody) = 3
             end if

             ! Set Up Data ...
             ! From a file
             if (File_Name(m) /= ' ' .and. File_Name(m) /= 'none') then
                call TAB_DATA_INPUT (m)
             end if

             ! (r,theta) data
             rtheta_flag = 0
             do n = 1,mtab
                if(Rtheta_Tabular_Pt(1,n) /= NULL_R) then
                   rtheta_flag = 1
                   do j = 1,mtab
                      RZ_Tabular_Pt(1,j) = Rtheta_Tabular_Pt(1,j)
                      RZ_Tabular_Pt(2,j) = Rtheta_Tabular_Pt(2,j)
                   end do
                   exit
                end if
             end do

             ! Count points and load
             do n = mtab,1,-1
                if(RZ_Tabular_Pt(1,n) /= NULL_R) then
                   Sgeom(3,m,nbody) = n
                   do l=1,n
                      if (Rtheta_flag /= 0) then
                         Rtab(l,nbody) = RZ_Tabular_Pt(1,l) &
                                       * COS(degree*RZ_Tabular_Pt(2,l))
                         Ztab(l,nbody) = RZ_Tabular_Pt(1,l) &
                                       * SIN(degree*RZ_Tabular_Pt(2,l))
                      else
                         Rtab(l,nbody) = RZ_Tabular_Pt(1,l)
                         Ztab(l,nbody) = RZ_Tabular_Pt(2,l)
                      end if
                   end do
                   cycle SURFACE
                end if
             end do

          else if (ABS(Isrftype(m,nbody)) == 8) then

             ! Read body from mesh file; make sure mesh_matnum has been assigned.
             if (Mesh_Matnum(nbody) == NULL_I) then
                write (message,126) nbody
126             format('mesh_material_number for body ',i2,' has not been assigned')
                call TLS_error (message)
                fatal = .true.
             end if
             cycle SURFACE

          end if

       end if
    end do SURFACE

    ! Print body data
300 continue

    if (nbody == 1) then
       call TLS_info ('')
       call TLS_info ('                                   Geometry Data')
       call TLS_info ('')
       call TLS_info ('          Body  Material  Surface  Translation  Rotation  Rotation  Surface')
       call TLS_info ('                 Number     Type      Point       Point     Angle  Parameters')
       call TLS_info ('          ----  --------  -------  -----------  --------  -------- ----------')
    end if

    do is = 1,Nsurf(nbody)
       if (is == 1) then
          l = 8
          write (message,310) nbody, Get_User_Material_ID(Matnum(nbody)), Surface_Name(is)(1:8), &
                                    Offset(1,is,nbody), Rotpt(1,is,nbody), &
                                    Rotangl(1,is,nbody), Sgeom(1,is,nbody)
310       format (11x,i2,6x,i2,4x,a,2x,1pe10.3,2x,1pe10.3,2x,0pf5.1,3x,1pe10.3)
          call TLS_info (message)
          write (message,315) Fill(is)(1:7), Offset(2,is,nbody), &
                                    Rotpt(2,is,nbody), Rotangl(2,is,nbody), &
                                    Sgeom(2,is,nbody)
315       format (25x,'(',a,')',1x,1pe10.3,2x,1pe10.3,2x,0pf5.1,3x,1pe10.3)
          call TLS_info (message)
          write (message,320) Offset(3,is,nbody), Rotpt(3,is,nbody), &
                                    Rotangl(3,is,nbody), Sgeom(3,is,nbody)
320       format (35x,1pe10.3,2x,1pe10.3,2x,0pf5.1,3x,1pe10.3)
          call TLS_info (message)
       else
          write (message,325) Surface_Name(is)(1:8), &
                                    Offset(1,is,nbody), Rotpt(1,is,nbody), &
                                    Rotangl(1,is,nbody), Sgeom(1,is,nbody)
325       format (25x,a,2x,1pe10.3,2x,1pe10.3,2x,0pf5.1,3x,1pe10.3)
          call TLS_info (message)
          write (message,315) Fill(is)(1:7), Offset(2,is,nbody), &
                                    Rotpt(2,is,nbody), Rotangl(2,is,nbody), &
                                    Sgeom(2,is,nbody)
          call TLS_info (message)
          write (message,320) Offset(3,is,nbody), Rotpt(3,is,nbody), &
                                    Rotangl(3,is,nbody), Sgeom(3,is,nbody)
          call TLS_info (message)
       end if

       if (ABS(Isrftype(is,nbody)) == 7) then
          write (message,330) Rtab(1,nbody),Ztab(1,nbody)
330       format(25x,'tabular points:',2x,1pe10.3,2x,1pe10.3)
          call TLS_info ('')
          call TLS_info (message)
          do n = 2,NINT(Sgeom(3,is,nbody))
             write (message,335) Rtab(n,nbody),Ztab(n,nbody)
335          format (42x,1pe10.3,2x,1pe10.3)
             call TLS_info (message)
          end do
          call TLS_info ('')
       end if
    end do

    if (Nsurf(nbody) == 0) then
       write (message,310) nbody, Get_User_Material_ID(Matnum(nbody)), Surface_Name(1)(1:10)
       call TLS_info (message)
    end if

  END SUBROUTINE BODY_CHECK

  SUBROUTINE BODY_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default BODY namelist.
    !
    !=======================================================================
    use interfaces_module,   only: Ar, Body_mass, Body_vel, Cosa, Offset, &
                                   old_Body_mass, Rotangl, Rotpt, Rtab,   &
                                   Sgeom, Sina, Surface_Name, Ztab

    ! Constants
    Ar = 0

    Body_mass     = 0
    Body_Vel      = 0
    Old_Body_Mass = 0
    Surface_Name  = ' '

    Cosa   = 0
    Offset = 0

    Rotangl = 0
    Rotpt   = 0
    Sgeom   = 0
    Sina    = 0

    Rtab = 0
    Ztab = 0

  END SUBROUTINE BODY_DEFAULT

  SUBROUTINE BODY_INPUT (lun, body_namelist)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read BODY namelist, put read data into place. This routine is
    !   called from GEN_INPUT At each call, data for body nbody is read.
    !   If a namelist is found, then a namelist is read and data is put
    !   into body data arrays at location nbody. On succesful read,
    !   then nbody -> nbody + 1.
    !
    !=======================================================================
    use interfaces_module,      only: Axis, body_phi, body_surfaces,            &
                                      body_temp, body_vel, File_Name, Fill,     &
                                      Height, material_number, Matnum, nbody,   &
                                      Nsurf, Offset, Rtheta_Tabular_Pt,         &
                                      RZ_Tabular_Pt, Radius, Rotpt, Rotangl,    &
                                      surface, Surface_Name, Tab_Type,          &
                                      Length, Rotation_Angle, Rotation_Pt,      &
                                      Translation_Pt, Mesh_Matnum,              &
                                      mesh_material_number
    use input_utilities,        only: seek_to_namelist, NULL_I, NULL_R, NULL_C
    use parallel_info_module,   only: p_info
    use parameter_module,       only: mbody, msurf, mtab, string_dim, nrot, mphi
    use legacy_mesh_api,        only: ndim
    use property_module,        only: Get_Truchas_Material_ID

    use scalar_func_factories,  only: alloc_const_scalar_func
    use scalar_func_table,      only: lookup_func

    ! Argument List
    integer, intent(in) :: lun
    logical :: body_namelist

    ! Local Variables
    character(80),  dimension(msurf) :: tabular_type
    character(120) :: fatal_error_string
    logical :: fatal
    integer :: ioerror, is, n
    real(r8) :: phi(mphi), temperature
    character(31) :: temperature_function
    real(r8), dimension(ndim) :: Velocity
    real(r8) :: scratch
    character(128) :: message

    ! Define BODY Namelist
    namelist /BODY/ Axis, Fill, Height, Radius, material_number, tabular_type,     &
                    phi, temperature, temperature_function, Velocity, File_Name,   &
                    Surface_Name, Rtheta_Tabular_Pt, RZ_Tabular_Pt, Length,        &
                    Rotation_Angle, Rotation_Pt, Translation_Pt, mesh_material_number

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    do n = 1, msurf
       tabular_type(n) = ''
    end do

    ! Initialize for error checking
    fatal = .false.
    fatal_error_string = 'BODY namelist input error!'

    ! Initialize the surface input
    body_surfaces        = 0
    material_number      = NULL_I
    mesh_material_number = NULL_I
    temperature          = NULL_R
    temperature_function = NULL_C
    phi                  = 0
    Velocity             = 0

    do is = 1,msurf
       Height(is)       = NULL_R
       Surface_Name(is) = 'none'
       Axis(is)         = 'none'
       Fill(is)         = 'none'
       File_Name(is)    = 'none'
       Tab_Type(is)     = 'none'

       do n = 1,nrot
          Rotation_Angle(n,is) = 0
       end do
       do n = 1,ndim
          Translation_Pt(n,is) = 0
          Rotation_Pt(n,is)    = 0
          Length(n,is)         = NULL_R
          Radius(n,is)         = NULL_R
       end do
    end do

    do is = 1,mtab
       do n = 1,2
          RZ_Tabular_Pt(n,is)     = NULL_R
          Rtheta_Tabular_Pt(n,is) = NULL_R
       end do
    end do

    ! Set error detection stuff
    IO_PE_ONLY: if (p_info%IOP) then
       ! Find namelist
       call seek_to_namelist (lun, 'BODY', found=body_namelist)

       ! Read namelist if found one
       if (body_namelist) then
          read (lun, NML = body, IOSTAT = ioerror)
          if (ioerror /= 0) then ! If read error, then didn't read namelist
             body_namelist = .false.
             fatal         = .true.
          end if
       end if

    end if IO_PE_ONLY

    ! Broadcast data just read in to all PE's.
    call BODY_INPUT_PARALLEL (body_namelist, phi, tabular_type, &
                              temperature, temperature_function, Velocity)

    ! Continue only if we found the namelist and read it without error.
    BODY_NML: if (body_namelist) then

       nbody = nbody + 1

       ! Read notice
       write (message, 15) nbody
15     format (' Reading BODY Namelist #',i2,' ...')
       call TLS_info ('')
       call TLS_info (message)

       ! Too many bodies is a fatal error
       if (nbody > mbody) then
          fatal = .true.
          write (fatal_error_string, 20) mbody
20        format('exceeded maximum number: mbody = ',i5)
       end if

       ! Assuming that we don't have too many bodies,
       ! Count surfaces
       if (.not. fatal) then
          do is = msurf,1,-1
             fatal = .false.
             if (Surface_Name(is) /= 'none') exit
             fatal = .true.
          end do
       end if

       ! If we got through msurf without finding a surface, then fatal == .true.
       if (fatal) then
          write (fatal_error_string, 25) nbody
25        format ('no surfaces for body ',i2)
       else
          body_surfaces = is
       end if

       ! Check that either the temperature or a function name was read, and not both.
       if (temperature == NULL_R .eqv. temperature_function == NULL_C) then
          fatal_error_string = 'either TEMPERATURE or TEMPERATURE_FUNCTION must be specified'
          fatal = .true.
       end if

       ! Save body data if we still don't have any error
       if (.not. fatal) then
         Nsurf(nbody)       = body_surfaces
         Matnum(nbody)      = Get_Truchas_Material_ID(material_number)
         Mesh_Matnum(nbody) = mesh_material_number
         
         ! Generate and store the body temperature function.
         if (temperature /= NULL_R) then
           call alloc_const_scalar_func (body_temp(nbody)%f, temperature)
         else
           call lookup_func (temperature_function, body_temp(nbody)%f)
           if (.not.allocated(body_temp(nbody)%f)) then
             fatal_error_string =  'unknown function name: ' // trim(temperature_function)
             fatal = .true.
           end if
         end if

         Body_Phi(nbody,:)  = phi
         do n = 1,ndim
           Body_Vel(n,nbody) = Velocity(n)
         end do
         
         do is = 1,Nsurf(nbody)
           Tab_Type(is)           = Tabular_Type(is)
           do n = 1,nrot
             Rotangl(n,is,nbody) = Rotation_Angle(n,is)
           end do
           do n = 1,ndim
             Offset(n,is,nbody)  = Translation_Pt(n,is)
             Rotpt(n,is,nbody)   = Rotation_Pt(n,is)
           end do
           Surface(is,nbody)      = Surface_Name(is)
         end do
         
       end if

    end if BODY_NML

    ! Now need to check for a fatal error
    call TLS_fatal_if_any (fatal, fatal_error_string)

  END SUBROUTINE BODY_INPUT

  SUBROUTINE BODY_INPUT_PARALLEL (body_namelist, phi, tabular_type, &
                                  temperature, temperature_function, velocity)
    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast items in BODY namelist from IO PE to all other PE's.
    !   Also broadcast body_namelist flag.
    !
    !======================================================================
    use interfaces_module,    only: Axis, File_Name, Fill, Height,        &
                                    material_number, Radius,              &
                                    Rtheta_Tabular_Pt, RZ_Tabular_Pt,     &
                                    Surface_Name, Length, Rotation_Angle, &
                                    Translation_Pt, Rotation_Pt,          &
                                    mesh_material_number
    use parameter_module,     only: msurf
    use legacy_mesh_api,      only: ndim
    use parallel_info_module, only: p_info
    use pgslib_module,        only: PGSLIB_BCAST

    ! Argument List
    character(80), dimension(msurf) :: tabular_type
    logical :: body_namelist
    real(r8) :: phi(:)
    real(r8) :: temperature
    character(*), intent(inout) :: temperature_function
    real(r8), dimension(ndim) :: Velocity
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Broadcast Data
    if (.NOT. p_info%UseGlobalServices) then
       call PGSLIB_BCAST (body_namelist)
       call PGSLIB_BCAST (phi)
       call PGSLIB_BCAST (tabular_type)
       call PGSLIB_BCAST (temperature)
       call PGSLIB_BCAST (temperature_function)
       call PGSLIB_BCAST (Velocity)

       call PGSLIB_BCAST (Axis)
       call PGSLIB_BCAST (Fill)
       call PGSLIB_BCAST (Height)
       call PGSLIB_BCAST (Radius)
       call PGSLIB_BCAST (material_number)
       call PGSLIB_BCAST (mesh_material_number)
       call PGSLIB_BCAST (File_Name)
       call PGSLIB_BCAST (Surface_Name)
       call PGSLIB_BCAST (Rtheta_Tabular_Pt)
       call PGSLIB_BCAST (RZ_Tabular_Pt)
       call PGSLIB_BCAST (Length)
       call PGSLIB_BCAST (Rotation_Angle)
       call PGSLIB_BCAST (Translation_Pt)
       call PGSLIB_BCAST (Rotation_Pt)
    end if

  END SUBROUTINE BODY_INPUT_PARALLEL


  SUBROUTINE TAB_DATA_INPUT (m)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read tabular data from TABULAR_DATA namelist.
    !
    !=======================================================================
    use interfaces_module, only: File_Name, RZ_Tabular_Pt, Rtheta_Tabular_Pt
    use input_utilities, only: seek_to_namelist

    ! Argument List
    integer :: m

    ! Local Variables
    logical :: file_exist, found
    integer :: l, tab_lun
    character(128) :: message

    ! Define TABULAR_DATA namelist
    namelist /TABULAR_DATA/ Rtheta_Tabular_Pt, RZ_Tabular_Pt

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Find File
10  l = LEN_TRIM(File_Name(m))
    inquire (FILE = File_Name(m)(1:l) , EXIST = file_exist)

    ! Read Data
    if (.not.file_exist) then
       call TLS_fatal ('file does not exist: ' // trim(file_name(m)))
    else
       ! Read Data
       open (NEWUNIT = tab_lun, FILE = File_Name(m)(1:l), STATUS = 'old')
       open (UNIT = tab_lun, FILE = File_Name(m)(1:l), STATUS = 'old')
       call seek_to_namelist (tab_lun, 'TABULAR_DATA', found)
       if (.not.found) go to 20
       read (tab_lun, tabular_data)
       close (tab_lun)
    end if

20  continue

  END SUBROUTINE TAB_DATA_INPUT

END MODULE BODY_INPUT_MODULE

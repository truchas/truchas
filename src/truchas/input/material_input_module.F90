!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MATERIAL_INPUT_MODULE
  !=======================================================================
  ! Purpose(s):
  !
  !   Define procedures for the input of material parameters.
  !
  !   Public Interface:
  !
  !     * call MATERIAL_INPUT ()
  !
  !       Defaults, reads, checks, and broadcasts input variables
  !       in the MATERIAL namelist(s).
  !
  !     * call MATERIAL_SIZES ()
  !
  !       Set the relevant problem material dimensions.
  !
  ! Contains: MATERIAL_CHECK
  !           MATERIAL_DEFAULT
  !           MATERIAL_INPUT
  !           MATERIAL_INPUT_PARALLEL
  !           MATERIAL_SIZES
  !
  ! Author(s): Jerry S. Brock (jsbrock@lanl.gov)
  !            Douglas B. Kothe (dbk@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use truchas_logging_services
  implicit none
  private

  public :: MATERIAL_INPUT, MATERIAL_SIZES

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

CONTAINS

  SUBROUTINE MATERIAL_CHECK (fatal)
    !=======================================================================
    ! Purpose(s):
    !
    !   Check MATERIAL namelist.
    !
    !=======================================================================
    use input_utilities,      only: NULL_I, NULL_R
    use matl_module,          only: relation_forms
    use parameter_module,     only: maxcon, maxmat, max_Relation_forms , nmat
    use property_data_module, only: background_material,                          &
                                    density, material_name, matpri, priority,     &
                                    Permeability_Constant,                        &
                                    Material_Feature,                             &
                                    Sound_Speed
    use utilities_module,     only: STRING_COMPARE
    use property_module,      only: Get_User_Material_ID

    ! Argument List
    logical :: fatal

    ! Local Variables
    character(80) :: initialized_string
    logical :: strings_match
    integer :: l, m, n
    character(128) :: message

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    fatal = .false.

    write (message,'(a,i0,a)') 'Identified ', nmat, ' material(s)'
    call TLS_info (message)

    MATERIALFEATURE: do m = 1,nmat

       if (Material_Feature(m) /= 'normal') then

          ! Background material.
          call STRING_COMPARE (Material_Feature(m), 'background', strings_match)
          if(strings_match) then
             background_material = m
             cycle MATERIALFEATURE 
          end if

          ! Invalid material.
          write (message,140) TRIM(Material_Feature(m)), Get_User_Material_ID(m)
140       format('Material_Feature = "',a,'" is invalid for material ',i2)
          call TLS_error (message)
          fatal = .true.

       end if

    end do MATERIALFEATURE

    ! Check background material number; punt and inform user if
    ! background material has not been set.
    if (background_material == NULL_I) then
       write (message, 14) 
14     format('Background material number is not defined; must have Material_Feature = ''background'' for one of the materials')
       call TLS_error (message)
       fatal = .true.
    else if (background_material <= 0 .or. background_material > maxmat) then
       write (message, 15) Get_User_Material_ID(background_material), maxmat
15     format('Background material number of ',i2,' is not between 0 and ',i2)
       call TLS_error (message)
       fatal = .true.
    end if


    ! Check material interface advection priority input
    matpri = 0  ! Clear list of materials in priority order

    do m = 1,nmat
       n = priority(m)
       if (n == NULL_I) then
          write (message, 20) Get_User_Material_ID(m),m
20        format('Material ',i2,' priority not initialized; set to ',i2)
          call TLS_warn (message)
          priority(m) = m
       else if (n < 1) then
          write (message, 25) Get_User_material_ID(m),n
25        format('Material priority(',i2,') = ',i2,' < 1')
          call TLS_error (message)
          fatal = .true.
       else if (n > nmat) then
          write (message, 30) Get_User_Material_ID(m),n,nmat
30        format('Material priority(',i2,') = ',i2,' > nmat = ',i2)
          call TLS_error (message)
          fatal = .true.
       end if
    end do

    ! Check for unique priority numbers; set matpri
    do m = 1,nmat
       do l = 1,nmat
          if (l == m) cycle
          if (priority(l) == priority(m)) then
             write (message, 35) Get_User_Material_ID(m), Get_User_Material_ID(l)
35           format('Materials ',i2,' and ',i2,' have the same priorities')
             call TLS_error (message)
             fatal = .true.
             exit
          end if
       end do
       if (fatal) cycle
       matpri(priority(m)) = m
    end do

    ! Print material priorities
    call TLS_info ('')
    call TLS_info ('               Material Priorities')
    call TLS_info ('')
    call TLS_info ('         Material     Name     Priority')
    call TLS_info ('         --------     ----     --------')

    do m = 1,nmat
       write (message, 45) Get_User_Material_ID(m), trim(material_name(m)), priority(m)
45     format(12x,i2,5x,a8,7x,i2)
       call TLS_info (message)
    end do

    ! Check property relation forms
    initialized_string = 'none'
    PROPERTY_RELATIONS: do m = 1,nmat
       ! Count the number of constants for each material
       ! and set appropriate flags

       !   User must specifiy the "reference" density
       if (density(m) == NULL_R) then
          write (message,69) Get_User_Material_ID(m)
69        format('Density must be specified for Material ',i2)
          call TLS_error (message)
          fatal = .true.
       end if

    end do PROPERTY_RELATIONS

    ! Make sure permeability constants are positive.
    PERMEABILITY: do m = 1,nmat

       ! Check constants.
       if (any(Permeability_Constant(:,m) < 0)) then
          write (message, 130) Get_User_Material_ID(m)
130       format('Material ',i2,' has a negative permeability constant; must be positive!')
          call TLS_error (message)
          fatal = .true.
       end if

    end do PERMEABILITY

    ! Make sure SOUND SPEEDS are positive.
    SOUNDSPEED: do m = 1,nmat

       ! Check constants.
       if (Sound_Speed(m) < 0) then
          write (message, 132) Get_User_Material_ID(m)
132       format('Material ',i2,' has a negative sound speed; must be positive!')
          call TLS_error (message)
          fatal = .true.
       end if

       if (Sound_Speed(m) > 0 .and. Density(m) /= 0) then
          write (message, 133) m
133       format('Material ',i2,' has a positive sound speed and is not a void')
          call TLS_warn (message)
       end if

    end do SOUNDSPEED

  END SUBROUTINE MATERIAL_CHECK

  SUBROUTINE MATERIAL_DEFAULT ()
    !=======================================================================
    ! Purpose(s):
    !
    !   Default MATERIAL namelist.
    !
    !=======================================================================
    use parameter_module,     only: nmat
    use input_utilities,      only: NULL_I, NULL_R
    use matl_module,          only: relation_forms
    use property_data_module, only: background_material,                          &
                                    density, material_name, matpri,               &
                                    Permeability_Constant, Material_Feature,      &
                                    Sound_Speed, Void_Temperature, isImmobile

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    nmat = 0

    ! Mobility flag
    isImmobile = .false.

    ! Character string defaults
    material_name                  = 'none'         
    Material_Feature               = 'normal'

    ! Material numbers
    background_material = NULL_I

    density                         = NULL_R
    Permeability_Constant           = 0

    Void_Temperature             = 0
    Sound_Speed                  = 0
    
    ! Valid character strings for property relations
    relation_forms = ''
    relation_forms(1) = 'constant'  
    relation_forms(2) = 'const'
    relation_forms(3) = 'temperature polynomial'
    relation_forms(4) = 'temperature_polynomial'
    relation_forms(5) = 'temp polynomial'

  END SUBROUTINE MATERIAL_DEFAULT

  SUBROUTINE MATERIAL_INPUT (lun)
    !=======================================================================
    ! Purpose(s):
    !
    !   Read MATERIAL namelist.
    !
    !=======================================================================
    use input_utilities,        only: seek_to_namelist, NULL_I, NULL_R
    use parallel_info_module,   only: p_info
    use parameter_module,       only: nmat
    use property_module,        only: Set_User_Material_ID
    use property_data_module,   only: density, material_name, priority,               &
                                      Permeability_Constant,                          &
                                      Material_Feature,                               &
                                      Sound_Speed, Void_Temperature, isImmobile
  
    integer, intent(in) :: lun

    ! Local Variables
    logical :: fatal, found, read_done
    integer :: material_number
    integer :: ioerror
    logical :: Immobile
    character(128) :: message

    ! Define MATERIAL Namelist
    namelist /material/ density, Immobile,                                &
                        material_name, material_number, priority,         &
                        Permeability_Constant,                            &
                        Material_Feature,                                 &
                        Sound_Speed, Void_Temperature
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    read_done = .false.

    ! Default namelist
    call MATERIAL_DEFAULT ()

    ! Input is done only on the IO PE.  Results of read are broadcast to all PE's
    ! For initialization, need to rewind input file.  Only do this once
    REWIND_IO_PE_ONLY: if (p_info%IOP) then
       rewind lun
    end if REWIND_IO_PE_ONLY

    ! Read namelist (reads multiple namelist block)
    ! This loops until either end of namelist or a fatal error is encountered
    nmat = 0
    READ_MATERIALS: do

       ! Re-initialize critical input values
       material_number           = NULL_I
       Immobile                  = .false.
       priority(0)               = NULL_I

       material_name(0)                  = 'none'
       Material_Feature(0)               = 'normal'

       density(0)                           = NULL_R
       Permeability_Constant(:,0)           = 0
       
       Void_Temperature(0)             = 0
       Sound_Speed(0)                  = 0

       ! Read next MATERIAL
       READ_IO_PE_ONLY: if (p_info%IOP) then

          call seek_to_namelist (lun, 'MATERIAL', found)
          read_done = .not.found
          if (read_done) then
             call TLS_info ('This is the last MATERIAL Namelist.')
          else
             ioerror = 0
             ! Inform usr that a material namelist is being read
             write (message, 10) nmat + 1
10           format (' Reading MATERIAL Namelist #',i2,' ...')
             call TLS_info ('')
             call TLS_info (message)

             ! Read the namelist
             read (lun, material, IOSTAT = ioerror)
          end if

       end if READ_IO_PE_ONLY

       ! Broadcast data just read to all PE's. Also broadcast the
       ! flags read_done and ioerror. If we didn't actually do a
       ! read, then this bcast is unnecessary (since default values
       ! are bcast). The code is cleaner, though, so do the extra work.
       ! Within material_input_parallel, only the 0th element of each
       ! array is broadcast, since that is the only element just read.
       call MATERIAL_INPUT_PARALLEL (ioerror, material_number, immobile, read_done)
       if (read_done) exit READ_MATERIALS

       if (ioerror /= 0) call TLS_fatal ('error reading MATERIAL namelist')

       ! New abbreviated input w/o indexing...temp, do in
       ! separate routine with checking & error msgs, etc. later

       nmat = nmat + 1
       MATERIAL_ASSIGNED: if (material_number == NULL_I) then
          material_number = nmat
          write (message, 15) nmat, nmat
15        format ('Assigned material_numer = ',i2,' to MATERIAL namelist ',i2)
          call TLS_info (message)
       end if MATERIAL_ASSIGNED

       ! Put input variables in proper material slot
       call Set_User_Material_ID(nmat, material_number)

       material_name(nmat)                     = material_name(0)
       priority(nmat)                          = priority(0)
       isImmobile(nmat)                        = immobile

       density(nmat)                           = density(0)

       Permeability_Constant(:,nmat)           = Permeability_Constant(:,0)

       Material_Feature(nmat)                  = Material_Feature(0)

       Void_Temperature(nmat)                  = Void_Temperature(0)
       Sound_Speed(nmat)                       = Sound_Speed(0)

       if (material_name(nmat) == 'none') then
          write (message, 25) nmat
25        format ('Material name in MATERIAL namelist ',i2,' not assigned')
          call TLS_fatal (message)
       end if

    end do READ_MATERIALS

    ! Check MATERIAL namelist
    call MATERIAL_CHECK (fatal)
    call TLS_fatal_if_any (fatal, 'terminating execution due to previous input errors')

  END SUBROUTINE MATERIAL_INPUT

  SUBROUTINE MATERIAL_INPUT_PARALLEL (ioerror, material_number, immobile, read_done)
    !======================================================================
    ! Purpose(s):
    !
    !   Broadcast all items in the materials namelist. Only 0th element
    !   in each array is read in, only bcast that element. Also broadcast
    !   the two flags read_done and ioerror.
    !
    !======================================================================
    use parallel_info_module, only: p_info
    use pgslib_module,        only: PGSLIB_BCAST
    use property_data_module, only: background_material,                            &
                                    density, material_name, priority,               &
                                    Permeability_Constant,                          &
                                    Material_Feature,                               &
                                    Sound_Speed, Void_Temperature

    ! Argument List
    logical, intent(inout) :: read_done
    logical, intent(inout) :: immobile
    integer, intent(inout) :: material_number
    integer, intent(inout) :: ioerror

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Broadcast Data
    BROADCAST_VARIABLES: if (.not. p_info%UseGlobalServices) then

       call PGSLIB_BCAST (background_material)

       call PGSLIB_BCAST (density(0))

       call PGSLIB_BCAST (ioerror)

       call PGSLIB_BCAST (immobile)
       call PGSLIB_BCAST (material_name(0))
       call PGSLIB_BCAST (material_number)
       call PGSLIB_BCAST (priority(0))

       call PGSLIB_BCAST (read_done)

       call PGSLIB_BCAST (Permeability_Constant(:,0))

       call PGSLIB_BCAST (Material_Feature(0))

       call PGSLIB_BCAST (Void_Temperature(0))
       call PGSLIB_BCAST (Sound_Speed(0))
       
    end if BROADCAST_VARIABLES
    
  END SUBROUTINE MATERIAL_INPUT_PARALLEL
  
  SUBROUTINE MATERIAL_SIZES ()
     !=======================================================================
     ! Purpose(s):
     !
     !   Set the relevant problem material dimensions.
     !
     !=======================================================================
     use parameter_module, only: mat_slot_new, mmat, nmat

     ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

     ! Set max number of materials.
     mmat = nmat

     ! Use 1 slot if we only have 1 material (= nmat); otherwise use 2
     ! slots. More slots are automatically allocated as needed.
     mat_slot_new = merge(1, 2, nmat <= 1)

  END SUBROUTINE MATERIAL_SIZES

end module MATERIAL_INPUT_MODULE




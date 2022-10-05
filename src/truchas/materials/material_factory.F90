!!
!! MATERIAL_FACTORY
!!
!! This module provides several procedures for creating MATERIAL class objects:
!!
!! * ALLOC_MATERIAL allocates a MATERIAL class object that is specified by a
!!   parameter list of a prescribed form; and
!!
!! * LOAD_MATERIAL_DATABASE allocates a parameter list-specified collection
!!   of MATERIAL class objects, storing them in a MATERIAL_DATABASE object.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! December 2019
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module material_factory

  use material_class
  use material_database_type
  use parameter_list_type
  use scalar_func_factories, only: alloc_scalar_func
  implicit none
  private

  public :: load_material_database, alloc_material

contains

  !! Load the materials specified by PARAMS into the material database MATL_DB.
  !! This can be called multiple times, though it is an error to attempt to
  !! load a material with the same name as one that already exists in MATL_DB.

  subroutine load_material_database(matl_db, params, stat, errmsg)

    class(material_database), intent(inout) :: matl_db
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
    type(parameter_list_iterator) :: piter
    class(material), allocatable :: matl

    piter = parameter_list_iterator(params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      call alloc_material(matl, piter%name(), plist, stat, errmsg)
      if (stat /= 0) return
      if (matl_db%has_matl(matl%name)) then
        stat = 1
        errmsg = 'attempt to overwrite existing material: ' // piter%name()
        return
      end if
      !! Write a solid fraction plot file for each phase change
      block
        use multiphase_matl_type
        use parallel_communication, only: is_IOP
        use truchas_env, only: output_dir
        integer :: n, i, ios
        character(:), allocatable :: filename
        select type (matl)
        class is (multiphase_matl)
          do n = 1, matl%num_phase() - 1
            filename = trim(output_dir) // matl%phase_name(n) // '-frac.dat'
            do ! until all blanks replaced with underscores
              i = scan(filename, ' ')
              if (i == 0) exit
              filename(i:i) = '_'
            end do
            if (is_IOP) call matl%write_solid_frac_plotfile(n, filename, digits=6, npoints=100, iostat=ios)
          end do
        end select
      end block
      call matl_db%add_matl(matl)
      call piter%next
    end do

  end subroutine load_material_database

  !! Allocate the MATERIAL class object specified by PARAMS.

  subroutine alloc_material(this, name, params, stat, errmsg)

    class(material), allocatable, intent(out) :: this
    character(*), intent(in) :: name
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    if (params%is_sublist('phases')) then
      call alloc_multiphase(this, name, params, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'ill-defined multi-phase material: ' // errmsg
        return
      end if
    else
      call alloc_single_phase(this, name, params, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'ill-defined single-phase material: ' // errmsg
        return
      end if
    end if

  end subroutine alloc_material

  !! Allocate the SINGLE_PHASE_MATL object specified by PARAMS and assign it
  !! the name NAME.

  subroutine alloc_single_phase(this, name, params, stat, errmsg)

    use single_phase_matl_type

    class(material), allocatable, intent(out) :: this
    character(*), intent(in) :: name
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(single_phase_matl), allocatable, target :: matl
    type(parameter_list), pointer :: plist

    allocate(matl)
    matl%name = name
    if (params%is_sublist('properties')) then
      plist => params%sublist('properties')
      call add_phase_properties(matl, plist, stat, errmsg)
      if (stat /= 0) return !TODO: refine errmsg
    else
      stat = 1
      errmsg = 'missing "properties" sublist parameter'
      return
    end if

    call move_alloc(matl, this)

  end subroutine alloc_single_phase

  !! Allocate the MULTIPHASE_MATL object specified by PARAMS and assign it
  !! the name NAME.

  subroutine alloc_multiphase(this, name, params, stat, errmsg)

    use multiphase_matl_type
    use phase_change_factory

    class(material), allocatable, intent(out) :: this
    character(*), intent(in) :: name
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n, nphase
    type(multiphase_matl), allocatable, target :: matl
    type(parameter_list), pointer :: plist
    type(parameter_list_iterator) :: piter
    type(parameter_list) :: phase_list
    character(:), allocatable :: string, pname

    allocate(matl)
    matl%name = name

    !! Properties spanning all phases
    if (params%is_sublist('properties')) then
      plist => params%sublist('properties')
      call add_phase_properties(matl, plist, stat, errmsg)
      if (stat /= 0) return !TODO: refine errmsg
    end if

    if (.not.params%is_sublist('phases')) then
      stat = 1
      errmsg = 'missing "phases" sublist parameter'
      return
    end if

    !! An initial sanity check over phase names, collecting the phase names in
    !! a temporary PHASE_LIST parameter list used below to order the phases.
    plist => params%sublist("phases")
    piter = parameter_list_iterator(plist, sublists_only=.true.)
    nphase = piter%count()
    if (nphase < 2) then
      stat = 1
      errmsg = '"phases" must define at least two phases'
      return
    end if
    do while (.not.piter%at_end())
      if (phase_list%is_sublist(piter%name())) then
        stat = 1
        errmsg = 'duplicate definition of phase ' // piter%name()
        return
      else
        plist => phase_list%sublist(piter%name())  ! create sublist for temp phase data
      end if
      call piter%next
    end do

    if (.not.params%is_sublist('phase-changes')) then
      stat = 1
      errmsg = 'missing "phase-changes" sublist parameter'
      return
    end if

    plist => params%sublist('phase-changes')
    piter = parameter_list_iterator(plist, sublists_only=.true.)
    if (piter%count() /= nphase-1) then
      stat = 1
      errmsg = '"phase-changes" not a sequence of phase changes connecting all phases'
      return
    end if

    !! Create a linked-list structure of phases, checking that the referenced
    !! phases are defined
    do while (.not.piter%at_end())
      string = piter%name()
      n = index(string, ':')
#ifdef INTEL_BUG20191228
      block
        character(:), allocatable :: low_phase, high_phase
        low_phase  = string(:n-1)
        high_phase = string(n+1:)
#else
      associate (low_phase => string(:n-1), high_phase => string(n+1:))
#endif
        if (phase_list%is_sublist(low_phase)) then
          plist => phase_list%sublist(low_phase)
          if (plist%is_parameter('next')) then
            stat = 1
            errmsg = '"phase-changes" not a sequence of phase changes'
            return
          else
            call plist%set('next', high_phase)
          end if
        else
          stat = 1
          errmsg = 'unknown phase in "phase-changes": ' // low_phase
          return
        end if
        if (phase_list%is_sublist(high_phase)) then
          plist => phase_list%sublist(high_phase)
          if (plist%is_parameter('prev')) then
            stat = 1
            errmsg = '"phase-changes" not a sequence of phase changes'
            return
          else
            call plist%set('prev', low_phase)
          end if
        else
          stat = 1
          errmsg = 'unknown phase in "phase-changes": ' // high_phase
          return
        end if
#ifdef INTEL_BUG20191228
      end block
#else
      end associate
#endif
      call piter%next
    end do

    !! Number the phases in order from low to high temperature.
    piter = parameter_list_iterator(phase_list, sublists_only=.true.)
    pname = piter%name() ! seed phase
    do ! walk back the sequence of phases to the beginning
      plist => phase_list%sublist(pname)
      if (.not.plist%is_parameter('prev')) exit
      call plist%get('prev', pname)
    end do
    n = 0
    do ! walk down the sequence of phases to the end and assign indexes
      n = n + 1
      call plist%set('index', n)
      if (.not.plist%is_parameter('next')) exit
      call plist%get('next', pname)
      plist => phase_list%sublist(pname)
    end do

    !! Initialize the PHASE objects.
    allocate(matl%phi(nphase))
    plist => params%sublist("phases")
    piter = parameter_list_iterator(plist, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => phase_list%sublist(piter%name())
      call plist%get('index', n)
      matl%phi(n)%name = piter%name()
      matl%phi(n)%matl => matl%phase
      plist => piter%sublist()
      call add_phase_properties(matl%phi(n), plist, stat, errmsg)
      if (stat /= 0) return !TODO: refine errmsg
      call piter%next
    end do

    !! Initialize the PHASE_CHANGE objects.
    allocate(matl%pc_seq(nphase-1))
    plist => params%sublist('phase-changes')
    piter = parameter_list_iterator(plist, sublists_only=.true.)
    do while (.not.piter%at_end())
      string = piter%name()
      n = index(string, ':')
      plist => phase_list%sublist(string(:n-1))
      call plist%get('index', n)
      plist => piter%sublist()
      call alloc_phase_change(matl%pc_seq(n)%pc, plist, stat, errmsg)
      if (stat /= 0) then
        errmsg = 'processing ' // plist%name() // ': ' // errmsg
        return
      end if
      call piter%next
    end do

    !! Check the phase changes do not overlap.
    do n = 1, size(matl%pc_seq) - 1
      if (matl%pc_seq(n)%pc%liquidus_temp() > matl%pc_seq(n+1)%pc%solidus_temp()) then
        stat = 1
        errmsg = 'processing ' // params%name() // ': overlapping phase changes'
        return
      end if
    end do

    call move_alloc(matl, this)

  end subroutine alloc_multiphase

  !! This auxiliary subroutine adds the properties and attributes specified by
  !! PARAMS to the PHASE class object.

  subroutine add_phase_properties(this, params, stat, errmsg)

    use scalar_func_class

    class(phase), intent(inout) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    character(:), allocatable :: pname
    class(scalar_func), allocatable :: f
    logical :: flag

    piter = parameter_list_iterator(params)
    do while (.not.piter%at_end())
      pname = piter%name()
#ifdef GNU_PR93762
      block
        character(:), allocatable :: dummy
        call params%get(pname, flag, stat=stat, errmsg=dummy)
      end block
#else
      call params%get(pname, flag, stat=stat)
#endif
      if (stat == 0) then ! this is an attribute
        if (flag) call this%add_attr(pname)
      else  ! this is a property
        call alloc_scalar_func(params, pname, f, stat, errmsg)
        if (stat /= 0) then
          errmsg = pname // ' property value error: ' // errmsg
          return
        end if
        call this%add_prop(pname, f)  !NB: this overwrites an existing function
      end if
      call piter%next
    end do

    stat = 0

  end subroutine add_phase_properties

end module material_factory

module thes_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use bndry_cfunc1_class
  implicit none
  private

  type, public :: thes_bc
    class(bndry_cfunc1), allocatable :: dirichlet!, neumann !, robin, robin_lhs
  contains
    procedure :: init
    procedure :: compute
  end type

  !! Interface for the boundary condition call back subroutine
  abstract interface
    subroutine bc_cb(plist, setids, stat, errmsg)
      use parameter_list_type
      implicit none
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
    end subroutine
  end interface

contains

  subroutine compute(this, t)
    class(thes_bc), intent(inout) :: this
    real(r8), intent(in) :: t
    call this%dirichlet%compute(t)
  end subroutine

  subroutine init(this, mesh, params, stat, errmsg)

    use simpl_mesh_type
    use parameter_list_type
    use bndry_node_cfunc_type

    class(thes_bc), intent(out) :: this
    type(simpl_mesh), intent(in), target :: mesh
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(bndry_node_cfunc), allocatable :: dir
    
    call iterate_list(params, 'potential', dir_proc, stat, errmsg)
    if (allocated(dir)) then
      call dir%add_complete
      call move_alloc(dir, this%dirichlet)
    end if
    
  contains
  
    subroutine dir_proc(plist, setids, stat, errmsg)
      use complex_scalar_func_factories
      type(parameter_list), intent(inout) :: plist
      integer, intent(in) :: setids(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      class(complex_scalar_func), allocatable :: f
      call alloc_complex_scalar_func(plist, 'phi', f, stat, errmsg)
      if (stat /= 0) return
      if (.not.allocated(dir)) then
        allocate(dir)
        call dir%init(mesh, compute_type=1, no_overlap=.true.)
      end if
      call dir%add(f, setids, stat, errmsg)
    end subroutine

  end subroutine

  !! This auxiliary subroutine iterates over the parameter list and for each BC
  !! sublist that matches the given TYPE, it calls the supplied subroutine
  !! PROC, which is expected to incrementally construct the BC using that data.

  subroutine iterate_list(params, type, proc, stat, errmsg)

    use parameter_list_type
    use string_utilities, only: lower_case

    type(parameter_list), intent(inout) :: params
    character(*), intent(in) :: type
    procedure(bc_cb) :: proc
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    integer, allocatable :: setids(:)
    character(:), allocatable :: this_type, context

    stat = 0
    piter = parameter_list_iterator(params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      context = 'processing ' // plist%path() // ': '
      call plist%get('type', this_type, stat, errmsg)
      if (stat /= 0) exit
      if (lower_case(this_type) == lower_case(type)) then  ! use this sublist
        call plist%get('face-set-ids', setids, stat, errmsg)
        if (stat /= 0) exit
        call proc(plist, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do
    if (stat /= 0) errmsg = context // errmsg

  end subroutine iterate_list

end module thes_bc_type

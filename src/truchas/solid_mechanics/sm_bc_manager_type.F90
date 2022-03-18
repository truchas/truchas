!!
!! SM_BC_TYPE
!!
!! TODO
!!
!! Zach Jibben <zjibben@lanl.gov>
!! March 2021
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_manager_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use timer_tree_type
  use truchas_logging_services
  use unstr_mesh_type
  use bndry_face_func_type
  use bndry_ip_func_type
  use sm_normal_traction_bc_type
  use sm_bc_list_type
  use sm_bc_face_list_type
  use sm_bc_node_list_type
  use sm_bc_class
  implicit none
  private

  type, public :: sm_bc_manager
    private
    logical, public :: contact_active

    real(r8) :: contact_penalty

    type(unstr_mesh), pointer :: mesh ! unowned reference

    type(sm_bc_list), pointer :: list => null()
    type(sm_bc_face_list) :: face_list
    type(sm_bc_node_list) :: node_list

    type(bndry_ip_func) :: traction(3)
    type(sm_normal_traction_bc) :: tractionn
    type(sm_bc_box), allocatable :: bcs(:)
  contains
    procedure :: init
    procedure :: apply
    procedure :: apply_traction
    procedure :: apply_nontraction
    procedure :: apply_deriv_diagonal
    procedure :: compute_viz_fields
    final :: delete_sm_bc_manager
  end type sm_bc_manager

contains

  subroutine delete_sm_bc_manager(this)
    type(sm_bc_manager), intent(inout) :: this
    if (associated(this%list)) deallocate(this%list)
  end subroutine delete_sm_bc_manager


  subroutine init(this, params, mesh, ig, contact_penalty, contact_distance, contact_traction, &
      stat, errmsg)

    use parallel_communication, only: global_sum
    use parameter_list_type
    use integration_geometry_type
    use sm_bc_c0d1_type
    use sm_bc_c0d2_type
    use sm_bc_c0d3_type
    use sm_bc_c1d0_type
    use sm_bc_c1d1_type
    use sm_bc_c1d2_type
    use sm_bc_c2d0_type
    use sm_bc_c2d1_type
    use sm_bc_c3d0_type

    class(sm_bc_manager), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    type(unstr_mesh), intent(in), target :: mesh
    type(integration_geometry), intent(in), target :: ig
    real(r8), intent(in) :: contact_penalty, contact_distance, contact_traction
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: b

    stat = 0
    this%mesh => mesh
    this%contact_penalty = contact_penalty

    allocate(this%list)
    call this%list%init(params, stat, errmsg)
    if (stat /= 0) return
    call this%face_list%init(mesh, this%list)
    call this%node_list%init(mesh, ig, this%list, this%face_list)

    allocate(this%bcs(9))
    allocate(sm_bc_c0d1 :: this%bcs(1)%p)
    allocate(sm_bc_c0d2 :: this%bcs(2)%p)
    allocate(sm_bc_c0d3 :: this%bcs(3)%p)
    allocate(sm_bc_c1d0 :: this%bcs(4)%p)
    allocate(sm_bc_c1d1 :: this%bcs(5)%p)
    allocate(sm_bc_c1d2 :: this%bcs(6)%p)
    allocate(sm_bc_c2d0 :: this%bcs(7)%p)
    allocate(sm_bc_c2d1 :: this%bcs(8)%p)
    allocate(sm_bc_c3d0 :: this%bcs(9)%p)

    this%contact_active = .false.
    do b = 1, size(this%bcs)
      call this%bcs(b)%p%init(mesh, this%node_list, this%list, &
          contact_penalty, contact_distance, contact_traction)
      if (b > 3) this%contact_active = this%contact_active .or. this%bcs(b)%p%enabled
    end do

    call alloc_bc('traction', this%traction)
    if (stat /= 0) return
    call alloc_tractionn_bc
    if (stat /= 0) return

    ! sanity check
    block
      character(128) :: msg
      integer :: n1, n2
      n1 = count(this%node_list%node <= mesh%nnode_onP)
      n1 = global_sum(n1)
      n2 = 0
      do b = 1, size(this%bcs)
        n2 = n2 + count(this%bcs(b)%p%index <= mesh%nnode_onP)
      end do
      n2 = global_sum(n2)
      write(msg,'("Nodes with requested BCs: ",i6,"    Nodes with applied BCs: ",i6)') n1, n2
      call TLS_info(trim(msg))
      if (n1 /= n2) call TLS_fatal("ERROR mismatching node BCs")
    end block

  contains

    subroutine alloc_bc(prefix, bc)

      character(*), intent(in) :: prefix
      type(bndry_ip_func), intent(out) :: bc(:)

      character(1), parameter :: dirstr(3) = ['x','y','z']
      integer :: d
      type(bndry_face_func), allocatable :: bff

      do d = 1, 3
        allocate(bff)
        call bff%init(mesh, bndry_only=.false.)
        call iterate_list(params, prefix//'-'//dirstr(d), prefix, bff, stat, errmsg)
        if (stat /= 0) return
        call bff%add_complete
        call bc(d)%init(mesh, ig, bff)
      end do

    end subroutine alloc_bc

    subroutine alloc_tractionn_bc
      type(bndry_face_func), allocatable :: bff
      allocate(bff)
      call bff%init(mesh, bndry_only=.false.)
      call iterate_list(params, 'traction-n', 'traction', bff, stat, errmsg)
      if (stat /= 0) return
      call bff%add_complete
      call this%tractionn%init(mesh, ig, bff)
    end subroutine alloc_tractionn_bc

  end subroutine init


  !! This auxiliary subroutine iterates over the parameter list and for each
  !! BC sublist that matches the given TYPE, it calls the supplied subroutine
  !! PROC, which is expected to incrementally construct the BC using that data.
  subroutine iterate_list(params, type, data_label, bff, stat, errmsg)

    use scalar_func_class
    use scalar_func_factories, only: alloc_scalar_func
    use string_utilities, only: lower_case

    type(parameter_list), intent(inout) :: params
    character(*), intent(in) :: type, data_label
    type(bndry_face_func), intent(inout) :: bff
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: plist
    integer, allocatable :: setids(:)
    character(:), allocatable :: this_type
    class(scalar_func), allocatable :: f

    stat = 0
    piter = parameter_list_iterator(params, sublists_only=.true.)
    do while (.not.piter%at_end())
      plist => piter%sublist()
      call plist%get('type', this_type, stat=stat, errmsg=errmsg)
      if (stat /= 0) exit
      if (lower_case(this_type) == type) then  ! use this sublist
        call TLS_info('  using SM_BC[' // piter%name() // ']')
        call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
        if (stat /= 0) exit
        call alloc_scalar_func(plist, data_label, f, stat, errmsg)
        if (stat /= 0) exit
        call bff%add(f, setids, stat, errmsg)
        if (stat /= 0) exit
      end if
      call piter%next
    end do
    if (stat /= 0) errmsg = 'SM_BC[' // piter%name() // ']: ' // errmsg

  end subroutine iterate_list


  subroutine apply(this, t, scaling_factor, displ, r)
    class(sm_bc_manager), intent(inout) :: this
    real(r8), intent(in) :: t, scaling_factor(:), displ(:,:)
    real(r8), intent(inout) :: r(:,:)
    ASSERT(size(displ,dim=2) == this%mesh%nnode .and. size(r,dim=2) == this%mesh%nnode)
    call this%apply_traction(t, r)
    call this%mesh%node_imap%gather_offp(r)
    call this%apply_nontraction(t, scaling_factor, displ, r)
  end subroutine apply


  !! Compute & apply boundary conditions to the residual.
  subroutine apply_nontraction(this, t, scaling_factor, displ, r)

    class(sm_bc_manager), intent(inout) :: this
    real(r8), intent(in) :: t, scaling_factor(:), displ(:,:)
    real(r8), intent(inout) :: r(:,:)

    integer :: b
    real(r8), allocatable :: ftot(:,:)

    call start_timer("BCs")

    ASSERT(size(displ,dim=2) == this%mesh%nnode .and. size(r,dim=2) == this%mesh%nnode)

    ! Make a copy of the sum of forces (residual without displacement BCs
    ! applied). This is needed to essentially avoid a race condition. Nodes in
    ! linked pairs need the forces on their sibling node, and they need to
    ! update their own residuals to account for gap forces.
    ftot = r

    do b = 1, size(this%bcs)
      call this%bcs(b)%p%apply(t, displ, ftot, scaling_factor, r)
    end do

    call stop_timer("BCs")

  end subroutine apply_nontraction


  subroutine apply_traction(this, t, r)

    class(sm_bc_manager), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(inout) :: r(:,:)

    integer :: d, i, n

    call start_timer("BCs")

    ! At traction BCs the stress is defined along the face. Each face has
    ! traction discretized at the integration points, one for each node. The
    ! values array provides the equivalent of traction = tensor_dot(stress,
    ! n) where n is normal of the boundary with magnitude equal to the area
    ! of this integration region (associated with node and face center, not
    ! the whole face).
    do d = 1, 3
      call this%traction(d)%compute(t)
      associate (nodes => this%traction(d)%index, values => this%traction(d)%value, &
          area => this%traction(d)%factor)
        do i = 1, size(nodes)
          n = nodes(i)
          if (n > this%mesh%nnode_onP) cycle
          r(d,n) = r(d,n) + values(i) * area(i)
        end do
      end associate
    end do

    ! Normal traction BCs apply traction in the normal direction, but zero
    ! traction in tangential directions. Each element of this BC contains a
    ! vector traction, which corresponds to the accumulated area-multiplied
    ! traction contribution to a node from its neighboring boundary
    ! integration points. The area is already factored into the values.
    call this%tractionn%compute(t)
    associate (nodes => this%tractionn%index, values => this%tractionn%value)
      do i = 1, size(nodes)
        n = nodes(i)
        if (n > this%mesh%nnode_onP) cycle
        r(:,n) = r(:,n) + values(:,i)
      end do
    end associate

    call stop_timer("BCs")

  end subroutine apply_traction


  !! Compute boundary conditions and modify the residual matrix diagonal
  !! accordingly.
  !!
  !! NB: Multi-contact and contact + displacement combo-BCs haven't yet had
  !! their contact-part preconditioner contributions implemented. These are
  !! only applying the displacement-part of the term.
  subroutine apply_deriv_diagonal(this, t, scaling_factor, displ, force, diag, F)

    class(sm_bc_manager), intent(inout) :: this
    real(r8), intent(in) :: t, scaling_factor(:), displ(:,:), force(:,:), F(:,:,:)
    real(r8), intent(inout) :: diag(:,:)

    integer :: b
    integer :: i, n

    ASSERT(size(displ,dim=2) == this%mesh%nnode .and. size(force,dim=2) == this%mesh%nnode)
    do b = 1, size(this%bcs)
      call this%bcs(b)%p%apply_deriv(t, displ, force, scaling_factor, F, diag)
    end do

  end subroutine apply_deriv_diagonal


  !! Copy gap displacement and gap traction fields into arrays for
  !! visualization.
  !!
  !! Right now this follows a rudimentary design: a single scalar
  !! displacement and traction at each gap node. The underlying model
  !! actually uses one displacement and one traction for each
  !! independent condition at a given node. So where two sidesets
  !! intersect, there will be two scalar displacements and two scalar
  !! tractions. Thus only one gap condition will be plotted at a given
  !! node. This is expected to be consistent, so the last gap listed in
  !! an input file is what gets plotted at collisions.
  !!
  !! If something more sophisticated is needed, we ought to implement a
  !! way of plotting fields over subsets of the mesh. This should be
  !! considered in any conversion to VTK. Otherwise we will need to dump
  !! a node-field for every gap sideset, which is wasteful.
  subroutine compute_viz_fields(this, gap_displacement, gap_traction)

    use,intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

    class(sm_bc_manager), intent(in) :: this
    real(r8), intent(out), allocatable :: gap_displacement(:), gap_traction(:)

    integer :: i, n, b

    allocate(gap_displacement(this%mesh%nnode_onP), gap_traction(this%mesh%nnode_onP))
    gap_displacement = ieee_value(0.0_r8, ieee_quiet_nan)
    gap_traction = ieee_value(0.0_r8, ieee_quiet_nan)

    do b = 1, size(this%bcs)
      ! only gap-contact conditions contribute something here
      if (.not.allocated(this%bcs(b)%p%displacement) &
          .or. .not.allocated(this%bcs(b)%p%traction)) cycle

      associate (nodes => this%bcs(b)%p%index, &
          displacement => this%bcs(b)%p%displacement, &
          traction => this%bcs(b)%p%traction)
        do i = 1, size(nodes)
          n = nodes(i)
          if (n > this%mesh%nnode_onP) cycle
          gap_displacement(n) = displacement(i)
          gap_traction(n) = traction(i)
        end do
      end associate
    end do

  end subroutine compute_viz_fields

end module sm_bc_manager_type

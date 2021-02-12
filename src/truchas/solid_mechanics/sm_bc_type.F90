!!
!! Zach Jibben <zjibben@lanl.gov>
!! September 2020
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module sm_bc_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parameter_list_type
  use truchas_logging_services
  use unstr_mesh_type
  use bndry_face_func_type
  use bndry_ip_func_type
  use sm_normal_displacement_bc_type
  use sm_normal_traction_bc_type
  use sm_gap_contact_bc_type
  implicit none
  private

  type, public :: sm_bc
    private
    logical, public :: gap_contact_active
    real(r8) :: contact_penalty

    type(unstr_mesh), pointer :: mesh ! unowned reference
    type(bndry_ip_func) :: displacement(3), traction(3)
    type(sm_normal_displacement_bc) :: displacementn
    type(sm_normal_traction_bc) :: tractionn
    type(sm_gap_contact_bc) :: gap_contact
  contains
    procedure :: init
    procedure :: apply
    procedure :: apply_traction
    procedure :: apply_deriv_diagonal
  end type sm_bc

contains

  subroutine init(this, params, mesh, ig, contact_penalty, contact_distance, contact_traction, &
      stat, errmsg)

    use parameter_list_type
    use integration_geometry_type

    class(sm_bc), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    type(unstr_mesh), intent(in), target :: mesh
    type(integration_geometry), intent(in), target :: ig
    real(r8), intent(in) :: contact_penalty, contact_distance, contact_traction
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    this%mesh => mesh
    this%contact_penalty = contact_penalty

    stat = 0
    call alloc_bc('displacement', this%displacement)
    if (stat /= 0) return
    call alloc_bc('traction', this%traction)
    if (stat /= 0) return
    call alloc_displacementn_bc
    if (stat /= 0) return
    call alloc_tractionn_bc
    if (stat /= 0) return
    call alloc_gap_contact_bc
    if (stat /= 0) return

    this%gap_contact_active = this%gap_contact%enabled ! TODO clean this up

    ! TODO: ensure consistency

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

    subroutine alloc_displacementn_bc
      type(bndry_face_func), allocatable :: bff
      allocate(bff)
      call bff%init(mesh, bndry_only=.false.)
      call iterate_list(params, 'displacement-n', 'displacement', bff, stat, errmsg)
      if (stat /= 0) return
      call bff%add_complete
      call this%displacementn%init(mesh, ig, bff)
    end subroutine alloc_displacementn_bc

    subroutine alloc_tractionn_bc
      type(bndry_face_func), allocatable :: bff
      allocate(bff)
      call bff%init(mesh, bndry_only=.false.)
      call iterate_list(params, 'traction-n', 'traction', bff, stat, errmsg)
      if (stat /= 0) return
      call bff%add_complete
      call this%tractionn%init(mesh, ig, bff)
    end subroutine alloc_tractionn_bc

    subroutine alloc_gap_contact_bc

      use scalar_func_class
      use scalar_func_factories, only: alloc_scalar_func
      use string_utilities, only: lower_case

      type(parameter_list_iterator) :: piter
      type(parameter_list), pointer :: plist
      integer, allocatable :: setids(:)
      character(:), allocatable :: this_type
      class(scalar_func), allocatable :: f

      call this%gap_contact%init(mesh, ig, contact_penalty, contact_distance, contact_traction)

      stat = 0
      piter = parameter_list_iterator(params, sublists_only=.true.)
      do while (.not.piter%at_end())
        plist => piter%sublist()
        call plist%get('type', this_type, stat=stat, errmsg=errmsg)
        if (stat /= 0) exit
        if (lower_case(this_type) == 'gap-contact') then  ! use this sublist
          call TLS_info('  using SM_BC[' // piter%name() // ']')
          call plist%get('face-set-ids', setids, stat=stat, errmsg=errmsg)
          if (stat /= 0) exit
          call alloc_scalar_func(plist, 'displacement', f, stat, errmsg)
          if (stat /= 0) exit
          call this%gap_contact%add(f, setids, stat, errmsg)
          if (stat /= 0) exit
        end if
        call piter%next
      end do
      if (stat /= 0) then
        errmsg = 'SM_BC[' // piter%name() // ']: ' // errmsg
        return
      end if

      call this%gap_contact%add_complete

    end subroutine alloc_gap_contact_bc

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


  !! Compute & apply boundary conditions to the residual.
  subroutine apply(this, t, scaling_factor, displ, r)

    class(sm_bc), intent(inout) :: this
    real(r8), intent(in) :: t, scaling_factor(:), displ(:,:)
    real(r8), intent(inout) :: r(:,:)

    call this%apply_traction(t, r)
    call gap_conditions
    call displacement_bcs

  contains

    ! Normal displacement gap condition. Enforces a given separation between two
    ! sides of a gap face, and possibly a shear force on each side.
    subroutine gap_conditions()

      integer :: i, n1, n2

      call this%gap_contact%compute(t, displ, r, scaling_factor)

      associate (link => this%gap_contact%index, values => this%gap_contact%value)

        if (this%gap_contact_active) then
          ! TODO-WARN: Need halo node displacements and stresses/residuals? For now just get
          !            everything working in serial.
#ifndef NDEBUG
          do i = 1, size(link, dim=2)
            n1 = link(1,i)
            n2 = link(2,i)
            if (n1 <= this%mesh%nnode_onP .or. n2 <= this%mesh%nnode_onP) then
              ASSERT(n1 <= this%mesh%nnode_onP .and. n2 <= this%mesh%nnode_onP)
            end if
          end do
#endif
        end if

        do i = 1, size(link, dim=2)
          n1 = link(1,i)
          n2 = link(2,i)
          if (n1 <= this%mesh%nnode_onP) r(:,n1) = r(:,n1) + values(:,1,i)
          if (n2 <= this%mesh%nnode_onP) r(:,n2) = r(:,n2) + values(:,2,i)
        end do

      end associate

    end subroutine gap_conditions

    subroutine displacement_bcs()

      integer :: d, i, n
      real(r8) :: displn

      ! Displacement BCs transform the matrix to a diagonal and the right hand
      ! side to the desired solution.
      do d = 1, 3
        call this%displacement(d)%compute(t)
        associate (nodes => this%displacement(d)%index, values => this%displacement(d)%value)
          do i = 1, size(nodes)
            n = nodes(i)
            r(d,n) = this%contact_penalty * (displ(d,n) - values(i)) * scaling_factor(n)
          end do
        end associate
      end do

      ! Normal displacement BCs enforce a given displacement in the normal
      ! direction, and apply zerotraction in tangential directions. To do this
      ! we must rotate the residual components components to a surface-aligned
      ! coordinate space, and set the normal component of the residual to the
      ! appropriate value. The node normal vector is computed from an
      ! area-weighted average of the surrounding integration points. After
      ! setting the appropriate coordinate-aligned value, it's important to
      ! rotate back into the usual coordinate space, for consistency with the
      ! preconditioner.
      call this%displacementn%compute(t)
      associate (nodes => this%displacementn%index, values => this%displacementn%value, &
          normal => this%displacementn%normal)
        do i = 1, size(nodes)
          n = nodes(i)
          displn = dot_product(displ(:,n), normal(:,i))
          r(:,n) = r(:,n) - dot_product(r(:,n), normal(:,i)) * normal(:,i)
          r(:,n) = r(:,n) - this%contact_penalty * (displn - values(i)) * normal(:,i) * scaling_factor(n)
        end do
      end associate

    end subroutine displacement_bcs

  end subroutine apply


  subroutine apply_traction(this, t, r)

    class(sm_bc), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(inout) :: r(:,:)

    integer :: d, i, n

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
        r(:,n) = r(:,n) + values(:,i)
      end do
    end associate

  end subroutine apply_traction


  !! Compute boundary conditions and modify the residual matrix diagonal
  !! accordingly.
  subroutine apply_deriv_diagonal(this, t, scaling_factor, displ, force, diag)

    class(sm_bc), intent(inout) :: this
    real(r8), intent(in) :: t, scaling_factor(:), displ(:,:), force(:,:)
    real(r8), intent(inout) :: diag(:,:)

    integer :: i, n, d

    call gap_conditions

    ! Dirichlet BCs
    do d = 1, 3
      call this%displacement(d)%compute(t)
      associate (nodes => this%displacement(d)%index)
        do i = 1, size(nodes)
          n = nodes(i)
          diag(d,n) = this%contact_penalty * scaling_factor(n)
        end do
      end associate
    end do

    ! The residual passed in will have, in the surface-aligned reference frame,
    ! the difference between the test normal displacement and the BC normal
    ! displacement. So we want the application of the preconditioner (inverse of
    ! diagonal) to result in this difference in that rotated reverence frame.
    call this%displacementn%compute(t)
    associate (nodes => this%displacementn%index, normal => this%displacementn%normal)
      do i = 1, size(nodes)
        n = nodes(i)
        diag(:,n) = diag(:,n) - diag(:,n) * normal(:,i)**2
        diag(:,n) = diag(:,n) - this%contact_penalty * normal(:,i)**2 * scaling_factor(n)
      end do
    end associate

  contains

    subroutine gap_conditions()

      integer :: n1, n2

      associate (link => this%gap_contact%index, dvalues => this%gap_contact%dvalue)

        if (this%gap_contact%enabled) then
          ! TODO-WARN: Need halo node displacements and stresses/residuals? For now just get
          !            everything working in serial.
          call this%gap_contact%compute_deriv(t, displ, force, scaling_factor, diag)
#ifndef NDEBUG
          do i = 1, size(link, dim=2)
            n1 = link(1,i)
            n2 = link(2,i)
            if (n1 <= this%mesh%nnode_onP .or. n2 <= this%mesh%nnode_onP) then
              ASSERT(n1 <= this%mesh%nnode_onP .and. n2 <= this%mesh%nnode_onP)
            end if
          end do
#endif
        end if

        do i = 1, size(link, dim=2)
          n1 = link(1,i)
          n2 = link(2,i)
          if (n1 <= this%mesh%nnode_onP) diag(:,n1) = diag(:,n1) + dvalues(:,1,i)
          if (n2 <= this%mesh%nnode_onP) diag(:,n2) = diag(:,n2) + dvalues(:,2,i)
        end do

      end associate

    end subroutine gap_conditions

  end subroutine apply_deriv_diagonal

end module sm_bc_type

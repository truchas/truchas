#include "f90_assert.fpp"

module fdme_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use simpl_mesh_type
  use material_model_type
  use avg_phase_prop_type
  use fdme_solver_type
  implicit none
  private

  type, public :: fdme_sim
    private
    type(simpl_mesh), pointer :: mesh => null()
    type(material_model) :: matl_model
    real(r8), allocatable :: vol_frac(:,:)
    type(avg_phase_prop) :: eps_prop, dlt_prop, mu_prop, sigma_prop
    real(r8), allocatable :: eps_re(:), eps_im(:), mu(:), sigma(:)
    real(r8) :: omega
    type(fdme_solver) :: solver
  contains
    procedure :: init
    procedure :: run
  end type

contains

  subroutine init(this, params, stat, errmsg)

    use parameter_list_type
    use simpl_mesh_factory, only: new_simpl_mesh
    use material_database_type
    use material_factory, only: load_material_database
    use material_utilities, only: define_property_default

    class(fdme_sim), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist, bodies_plist
    type(material_database) :: matl_db
    character(:), allocatable :: context, matl_names(:)
    real(r8) :: freq
    integer :: j

    !! Create the mesh
    if (params%is_sublist('mesh')) then
      plist => params%sublist('mesh')
      this%mesh => new_simpl_mesh(plist, stat, errmsg)
      context = 'processing ' // plist%path() // ': '
      if (stat /= 0) errmsg = context // errmsg
    else
      stat = 1
      errmsg = 'missing "mesh" sublist parameter'
    end if
    if (stat /= 0) return

    !! Load the material database
    if (params%is_sublist('materials')) then
      plist => params%sublist('materials')
      context = 'processing ' // plist%path() // ': '
      call load_material_database(matl_db, plist, stat, errmsg)
      if (stat /= 0) errmsg = context // errmsg
    else
      stat = 1
      errmsg = 'missing "materials" sublist parameter'
    end if
    if (stat /= 0) return

    if (params%is_sublist('bodies')) then
      bodies_plist => params%sublist('bodies')
    else
      stat = 1
      errmsg = 'missing "bodies" sublist parameter'
      if (stat /= 0) return
    end if

    !! Extract the material names from the body sublists and initialize the material model
    call get_matl_names(bodies_plist, matl_names, stat, errmsg, unique=.true.)
    if (stat /= 0) return
    call this%matl_model%init(matl_names, matl_db, stat, errmsg)
    if (stat /= 0) return

    !! Only single-phase materials are currently supported
    do j = 1, this%matl_model%nmatl_real
      if (this%matl_model%num_matl_phase(j) /= 1) then
        stat = 1
        errmsg = 'found unsupported multiphase material: ' // this%matl_model%matl_name(j)
        return
      end if
    end do

    !! Initialize the material volume fraction array
    call get_matl_vol_frac(this%mesh, this%matl_model, bodies_plist, this%vol_frac, stat, errmsg)

    !! Set the default property values if not specified by input
    call define_property_default(this%matl_model, 'relative-permittivity',   1.0_r8)
    call define_property_default(this%matl_model, 'relative-permeability',   1.0_r8)
    call define_property_default(this%matl_model, 'dielectric-loss-tangent', 0.0_r8)
    call define_property_default(this%matl_model, 'sigma', 0.0_r8)

    !! Create the property objects
    call this%eps_prop%init('relative-permittivity', this%matl_model, stat, errmsg, void_value=1.0_r8)
    if (stat /= 0) return

    call this%mu_prop%init('relative-permeability', this%matl_model, stat, errmsg, void_value=1.0_r8)
    if (stat /= 0) return

    call this%dlt_prop%init('dielectric-loss-tangent', this%matl_model, stat, errmsg, void_value=0.0_r8)
    if (stat /= 0) return

    call this%sigma_prop%init('sigma', this%matl_model, stat, errmsg, void_value=0.0_r8)
    if (stat /= 0) return

    !! Compute the cell-based property arrays
    block
      real(r8) :: state(0) ! state variables e.g. (T, x, y, z) would go here
      associate (n => this%mesh%ncell)
        allocate(this%eps_re(n), this%eps_im(n), this%mu(n), this%sigma(n))
      end associate
      do j = 1, this%mesh%ncell
        call this%eps_prop%compute_value(this%vol_frac(:,j), state, this%eps_re(j))
        call this%dlt_prop%compute_value(this%vol_frac(:,j), state, this%eps_im(j))
        this%eps_im(j) = this%eps_im(j) * this%eps_re(j)
        call this%mu_prop%compute_value(this%vol_frac(:,j), state, this%mu(j))
        call this%sigma_prop%compute_value(this%vol_frac(:,j), state, this%sigma(j))
      end do
    end block

    !! Initialize the solver
    if (params%is_sublist('solver')) then
      plist => params%sublist('solver')
      call this%solver%init(this%mesh, plist, stat, errmsg)
    else
      stat = 1
      errmsg = 'missing "solver" sublist parameter'
    end if
    if (stat /= 0) return

    call plist%get('frequency', freq, stat, errmsg)
    if (stat /= 0) return
    if (freq <= 0) then
      stat = 1
      errmsg = 'non-positive frequency'
      return
    end if
    this%omega = 8*atan(1.0_r8)*freq

    call this%solver%setup(this%omega, this%eps_re, this%eps_im, this%mu, this%sigma, stat, errmsg)
    if (stat /= 0) return

  contains

    !! Return the array of body material names appearing in the bodies parameter
    !! list and an array of the unique material names.

    subroutine get_matl_names(bodies, matl_names, stat, errmsg, unique)

      type(parameter_list), intent(inout) :: bodies
      character(:), allocatable, intent(out) :: matl_names(:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg
      logical, intent(in), optional :: unique

      type(parameter_list_iterator) :: piter
      type(parameter_list), pointer :: plist
      type(parameter_list) :: matl_list
      character(:), allocatable :: name
      logical :: unique_
      integer :: n

      unique_ = .false.
      if (present(unique)) unique_= unique

      !! Scan the material names found in the body sublists and verify that
      !! they appear in the material database. Also find the length of the
      !! longest name, which will be used to allocate the arrays. If returning
      !! an array of the unique names, we use a temporary parameter list to
      !! generate a list of the unique names.

      piter = parameter_list_iterator(bodies, sublists_only=.true.)
      n = 0 ! max name length
      do while (.not.piter%at_end())
        plist => piter%sublist()
        context = 'processing ' // plist%path() // ': '
        call plist%get('material', name, stat, errmsg)
        if (stat /= 0) then
          errmsg = context // errmsg
          return
        else if (name /= 'VOID' .and. .not.matl_db%has_matl(name)) then
          stat = 1
          errmsg = context // 'unknown "material": ' // name
          return
        end if
        n = max(n, len(name))
        if (unique_) call matl_list%set(name, 1)
        call piter%next
      end do

      if (unique_) then ! copy the list of unique material names into the array

        piter = parameter_list_iterator(matl_list)
        allocate(character(n) :: matl_names(piter%count()))
        n = 0
        do while (.not.piter%at_end())
          n = n + 1
          matl_names(n) = piter%name()
          call piter%next
        end do

      else ! copy the body material names in order into the array

        piter = parameter_list_iterator(bodies, sublists_only=.true.)
        allocate(character(n) :: matl_names(piter%count()))
        n = 0
        do while (.not.piter%at_end())
          plist => piter%sublist()
          call plist%get('material', name)
          n = n + 1
          matl_names(n) = name
          call piter%next
        end do

      end if

    end subroutine get_matl_names

    !! Return the cell material volume fraction array that is defined by the
    !! body sublists of the bodies parameter list.

    subroutine get_matl_vol_frac(mesh, matl_model, bodies, vol_frac, stat, errmsg)

      use compute_body_volumes_proc, only: compute_body_volumes

      type(simpl_mesh), intent(in) :: mesh
      type(material_model), intent(in) :: matl_model
      type(parameter_list), intent(inout) :: bodies
      real(r8), allocatable :: vol_frac(:,:)
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      integer :: num_matl, num_body, j, m
      real(r8), allocatable :: body_vol(:,:)
      integer, allocatable :: matl_index(:)
      logical, allocatable :: body_matl_mask(:,:)
      character(:), allocatable :: matl_names(:)

      call compute_body_volumes(this%mesh, 3, bodies, body_vol, stat, errmsg)
      context = 'processing ' // bodies%path() // ': '
      if (stat /= 0) then
        errmsg = context // errmsg
        return
      end if

      num_matl = matl_model%nmatl
      num_body = size(body_vol,dim=1)

      !! BODY_MATL_MASK(b,m) is true if body b is material m.
      call get_matl_names(bodies, matl_names, stat, errmsg)
      if (stat /= 0) return
      matl_index = matl_model%matl_index(matl_names)
      allocate(body_matl_mask(num_body,num_matl), source=.false.)
      do j = 1, num_body
        body_matl_mask(j, matl_index(j)) = .true.
      end do

      !! Initialize the material volume fraction array.
      allocate(vol_frac(num_matl,mesh%ncell))
      do j = 1, mesh%ncell_onP
        do m = 1, num_matl
          vol_frac(m,j) = sum(body_vol(:,j), mask=body_matl_mask(:,m)) / abs(mesh%volume(j))
        end do
      end do
      call mesh%cell_imap%gather_offp(vol_frac)

    end subroutine get_matl_vol_frac

  end subroutine init

  subroutine run(this, stat, errmsg)

    use fdme_vtk_graphics_proc

    class(fdme_sim), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    call this%solver%solve(stat, errmsg)
    if (stat /= 0) return

    call fdme_vtk_graphics(this%solver, 'out.vtkhdf', stat, errmsg)

  end subroutine run

end module fdme_sim_type

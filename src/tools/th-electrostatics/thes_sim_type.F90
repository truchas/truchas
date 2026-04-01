#include "f90_assert.fpp"

module thes_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use simpl_mesh_type
  use material_model_type
  use avg_phase_prop_type
  use thes_solver_type
  implicit none
  private

  type, public :: thes_sim
    private
    type(simpl_mesh), pointer :: mesh => null()
    type(material_model) :: matl_model
    real(r8), allocatable :: vol_frac(:,:)
    type(avg_phase_prop) :: eps_prop, eps_im_prop
    complex(r8), allocatable :: eps(:), phi(:)
    type(thes_solver) :: solver
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
    use thes_bc_type

    class(thes_sim), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist, bodies_plist
    type(material_database) :: matl_db
    character(:), allocatable :: context, matl_names(:)
    type(thes_bc) :: bc
    real(r8) :: eps0
    integer :: j

    !! Read physical constants (optional)
    plist => params%sublist('physical-constants', stat, errmsg)
    if (stat /= 0) return
    call plist%get('vacuum-permittivity', eps0, stat, errmsg, default= 8.854188e-12_r8)

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

    !!
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

    !! Initialize the cell-based permittivity array

    !! Set the default permittivity value if not specified by input
    call define_property_default(this%matl_model, 'relative-permittivity', 1.0_r8)
    call define_property_default(this%matl_model, 'relative-permittivity-im', 0.0_r8)

    call this%eps_prop%init('relative-permittivity', this%matl_model, stat, errmsg, void_value=1.0_r8)
    if (stat /= 0) return

    call this%eps_im_prop%init('relative-permittivity-im', this%matl_model, stat, errmsg)
    if (stat /= 0) return

    block
      real(r8) :: state(0) ! state variables e.g. (T, x, y, z) would go here
      allocate(this%eps(this%mesh%ncell))
#if defined(INTEL_BUG20250605) || defined(GNU_PR119986)
      associate (eps => this%eps) ! works around an ifx 2025.1 bug
        do j = 1, this%mesh%ncell
          call this%eps_prop%compute_value(this%vol_frac(:,j), state, eps(j)%re)
          call this%eps_im_prop%compute_value(this%vol_frac(:,j), state, eps(j)%im)
          eps(j) = eps0 * eps(j)
        end do
      end associate
#else
      do j = 1, this%mesh%ncell
        call this%eps_prop%compute_value(this%vol_frac(:,j), state, this%eps(j)%re)
        call this%eps_im_prop%compute_value(this%vol_frac(:,j), state, this%eps(j)%im)
        this%eps(j) = eps0 * this%eps(j)
      end do
#endif
    end block

    !! Initialize the boundary condition data
    if (params%is_sublist('bc')) then
      plist => params%sublist('bc')
      call bc%init(this%mesh, plist, stat, errmsg)
      context = 'processing ' // plist%path() // ': '
      if (stat /= 0) errmsg = context // errmsg
      call bc%compute(t=0.0_r8)
    else
      stat = 1
      errmsg = 'missing "bc" sublist parameter'
    end if
    if (stat /= 0) return

    !! Initialize the solver
    if (params%is_sublist('solver')) then
      plist => params%sublist('solver')
      call this%solver%init(this%mesh, this%eps, bc, plist, stat, errmsg)
    else
      stat = 1
      errmsg = 'missing "solver" sublist parameter'
    end if
    if (stat /= 0) return

    allocate(this%phi(this%mesh%nnode))

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

    class(thes_sim), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    call this%solver%solve(this%phi, stat, errmsg)
    if (stat /= 0) return

    call this%mesh%node_imap%gather_offp(this%phi)
    call write_vtk_graphics(this, 'out.vtkhdf', stat, errmsg)

  end subroutine run


  subroutine write_vtk_graphics(this, filename, stat, errmsg)

    use,intrinsic :: iso_fortran_env, only: int8
    use vtkhdf_mb_file_type
    use vtkhdf_vtk_cell_types, only: VTK_TETRA

    class(thes_sim), intent(in) :: this
    character(*), intent(in) :: filename
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(vtkhdf_mb_file) :: viz_file
    type(vtkhdf_block_handle) :: handle
    real(r8), allocatable :: x(:,:), vector(:,:)
    complex(r8), allocatable :: zscalar(:), zvector(:,:)
    complex(r8) :: grad_phi(3,this%mesh%ncell)
    integer, allocatable :: xcnode(:), cnode(:), block_cells(:), block_nodes(:)
    integer :: n, bitmask
    integer(int8), allocatable :: types(:)

    call viz_file%create(filename, comm, stat, errmsg)
    if (stat /= 0) return

    block ! Cell averaged electric field phasor
      use mimetic_discretization, only: grad, w1_vector_on_cells
      real(r8) :: w1(this%mesh%nedge)
#ifdef GNU_PR119986
      call grad(this%mesh, (this%phi%re), w1)
#else
      call grad(this%mesh, this%phi%re, w1)
#endif
      grad_phi%re = w1_vector_on_cells(this%mesh, w1)
#ifdef GNU_PR119986
      call grad(this%mesh, (this%phi%im), w1)
#else
      call grad(this%mesh, this%phi%im, w1)
#endif
      grad_phi%im = w1_vector_on_cells(this%mesh, w1)
    end block

    do n = 1, size(this%mesh%cell_set_name)
      associate (name => this%mesh%cell_set_name(n)%s)
        handle = viz_file%add_block(name, mode=UG_STATIC)
        bitmask = ibset(0,pos=n)
        call this%mesh%get_mesh_block(bitmask, x, xcnode, cnode, block_cells, block_nodes)
        types = spread(VTK_TETRA, dim=1, ncopies=size(xcnode)-1)
        call viz_file%write_mesh(handle, x, cnode, xcnode, types)

        !! Cell-based material volume fractions
        vector = this%vol_frac(:,block_cells)
        call viz_file%write_cell_data(handle, 'vol-frac', vector)

        !! Cell-based complex permittivities
        zscalar = this%eps(block_cells)
        call viz_file%write_cell_data(handle, 'eps_re', zscalar%re)
#ifdef GNU_PR117774
        call viz_file%write_cell_data(handle, 'eps_im', [zscalar%im])
#else
        call viz_file%write_cell_data(handle, 'eps_im', zscalar%im)
#endif

        !! Node-based complex electric potential phasor
        zscalar = this%phi(block_nodes)
        call viz_file%write_point_data(handle, 'phi_re', zscalar%re)
#ifdef GNU_PR117774
        call viz_file%write_point_data(handle, 'phi_im', [zscalar%im])
#else
        call viz_file%write_point_data(handle, 'phi_im', zscalar%im)
#endif
        call viz_file%write_point_data(handle, '|phi|', abs(zscalar))

        !! Cell-averaged electric field, grad(phi)
        zvector = grad_phi(:,block_cells)
        call viz_file%write_cell_data(handle, 'grad phi_re', zvector%re)
#ifdef GNU_PR117774
        call viz_file%write_cell_data(handle, 'grad phi_im', reshape([zvector%im],shape(zvector)))
#else
        call viz_file%write_cell_data(handle, 'grad phi_im', zvector%im)
#endif

        call viz_file%write_cell_data(handle, 'MPI rank', spread(this_PE,dim=1,ncopies=size(block_cells)))
      end associate
    end do

    call viz_file%close

  end subroutine write_vtk_graphics

end module thes_sim_type

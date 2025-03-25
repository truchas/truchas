#include "f90_assert.fpp"

module th_electrostatics_sim_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use parallel_communication
  use simpl_mesh_type
  use material_model_type
  implicit none
  private

  type, public :: th_electrostatics_sim
    private
    type(simpl_mesh), pointer :: mesh => null()
    type(material_model) :: matl_model
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
    use compute_body_volumes_proc, only: compute_body_volumes

    class(th_electrostatics_sim), intent(out) :: this
    type(parameter_list), intent(inout) :: params
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(parameter_list), pointer :: plist
    type(material_database) :: matl_db
    character(:), allocatable :: context
    real(r8), allocatable :: vof(:,:)

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

    !! Initialize the material model
    !TODO: input name instead of hardwiring it
    call this%matl_model%init(['default'], matl_db, stat, errmsg)
    if (stat /= 0) errmsg = context // errmsg

    if (params%is_sublist('bodies')) then
      plist => params%sublist('bodies')
      call compute_body_volumes(this%mesh, 5, plist, vof, stat, errmsg)
      context = 'processing ' // plist%path() // ': '
      if (stat /= 0) errmsg = context // errmsg
    else
      stat = 1
      errmsg = 'missing "bodies" sublist parameter'
    end if
    if (stat /= 0) return

  end subroutine init

  subroutine run(this, stat, errmsg)
  
    class(th_electrostatics_sim), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
  
    call write_vtk_graphics('out.vtkhdf', this%mesh, stat, errmsg)

  end subroutine run
  

  subroutine write_vtk_graphics(filename, mesh, stat, errmsg)
  
    use vtkhdf_file_type
    
    character(*), intent(in) :: filename
    type(simpl_mesh), intent(in) :: mesh
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    
    type(vtkhdf_file) :: viz_file
    
    if (is_IOP) call viz_file%create(filename, stat, errmsg)
    call broadcast(stat)
    if (stat /= 0) then
      call broadcast(errmsg)
      return
    end if
    
    call write_mesh
    
  contains
  
    subroutine write_mesh

      use,intrinsic :: iso_fortran_env, only: int8

      integer, allocatable, target :: cnode(:,:)
      integer, allocatable :: xcnode(:)
      integer(int8), allocatable :: types(:)
      real(r8), allocatable :: x(:,:)
      integer, pointer :: connectivity(:)
      integer :: j, stat
      character(:), allocatable :: errmsg

      !! Collate the mesh data structure onto the IO process
      call mesh%get_global_cnode_array(cnode)
      call mesh%get_global_x_array(x)

      if (is_IOP) then
        xcnode = [(1+4*j, j=0, size(cnode,dim=2))]
        connectivity(1:size(cnode)) => cnode ! flattened view
        types = spread(VTK_TETRA, dim=1, ncopies=size(cnode,dim=2))
        call viz_file%write_mesh(x, connectivity, xcnode, types, stat, errmsg)
      end if
      call broadcast(stat)
      INSIST(stat == 0)

    end subroutine write_mesh

  end subroutine write_vtk_graphics
    
end module th_electrostatics_sim_type

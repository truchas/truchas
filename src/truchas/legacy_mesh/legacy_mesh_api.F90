!!
!! LEGACY_MESH_API
!!
!! This module provides the legacy API to the main Truchas mesh instantiated as
!! an UNSTR_MESH object.  This is intended to serve as a temporary bridge while
!! client code is transitioned to use the native API of the next generation
!! UNSTR_MESH mesh data structure.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Ondrej Certik <certik@lanl.gov>
!! Matt Calef <mcalef@lanl.gov>
!! February 2016
!!

module legacy_mesh_api

  use common_impl, only: unpermute_mesh_vector, unpermute_vertex_vector
  use common_impl, only: ncells, nnodes, nnodes_tot, ncells_tot
  use ee_gather_impl, only: ee_gather, gather_boundarydata
  use en_gather_impl, only: en_gather, gather_vertex_coord, en_min_gather, en_or_gather
  use en_gather_impl, only: en_sum_scatter, en_or_scatter, en_min_scatter, en_max_scatter
  use mesh_impl, only: mesh, DEGENERATE_FACE
  use mesh_impl, only: is_face_ngbr, mesh_collate_vertex
  use cell_impl, only: cell, orthogonal_mesh
  use nn_gather_impl, only: vertex_ngbr_all, vertex_ngbr_all_orig, nn_gather_boundarydata
  use vertex_impl, only: vertex_data, vertex, vertex_collate
  use mesh_face_set_impl, only: mesh_face_set
  use mesh_impl, only: CELL_TET, CELL_PYRAMID, CELL_PRISM, CELL_HEX
  use mesh_impl, only: GAP_ELEMENT_1, GAP_ELEMENT_3, GAP_ELEMENT_5
  use cell_impl, only: linear_prop
  implicit none
  public

  logical, parameter :: mesh_has_cblockid_data = .true.
  integer, parameter :: ndim=3, nfc=6, nvc=8, nvf=4, nfv=3, nec=12

  ! Define the edges surrounding a hex cell
  integer :: Cell_Edge(2,12)
  data Cell_Edge/1,2, 2,3, 3,4, 4,1, 2,6, 3,7, 4,8, 1,5, 5,6, 6,7, 7,8, 8,5/

  ! Face-Vertex and Vertex-Face mappings.
  integer :: Face_Vrtx(nfc,nvf)
  data Face_Vrtx /3,1,4,2,3,8,  4,2,1,3,2,5,  8,6,5,7,1,6,  7,5,8,6,4,7/

  integer :: Vrtx_Face(nvc,ndim)
  data Vrtx_Face /2,2,1,1,2,2,1,1,  3,4,4,3,3,4,4,3,  5,5,5,5,6,6,6,6/

contains

  subroutine init_legacy_mesh_api

    use common_impl, only: init_common_impl
    use mesh_impl, only: init_mesh_impl
    use ee_gather_impl, only: init_ee_gather_impl
    use en_gather_impl, only: init_en_gather_impl
    use nn_gather_impl, only: init_nn_gather_impl
    use vertex_impl, only: init_vertex_impl
    use cell_impl, only: init_cell_impl
    use mesh_face_set_impl, only: init_mesh_face_set_impl

    call init_common_impl
    call init_mesh_impl
    call init_vertex_impl
    call init_ee_gather_impl
    call init_en_gather_impl
    call init_nn_gather_impl
    call init_cell_impl
    call init_mesh_face_set_impl

  end subroutine init_legacy_mesh_api

end module legacy_mesh_api

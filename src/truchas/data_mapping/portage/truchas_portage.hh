// -----------------------------------------------------------------------------
// This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
// -----------------------------------------------------------------------------
//
// This implements a custom interface to the Portage data remapping driver
// that is adapted to Truchas mesh data structures. This handles part-by-part
// remapping.
//
// Neil N. Carlson <nnc@lanl.gov> with help of Gary Dilts and Rao Garimella
// November 2019
//
// This works with tag v2.2.0 of Portage: https://github.com/laristra/portage
// NOTE: SERIAL IMPLEMENTATION ONLY


#ifndef TRUCHAS_PORTAGE_H_
#define TRUCHAS_PORTAGE_H_

#include "portage/support/portage.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_rNd.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "wonton/support/Point.h"
#include "wonton/mesh/AuxMeshTopology.h"
#include "portage/driver/coredriver.h"

#include <stdio.h>

extern "C" {
typedef struct {
  int num_cell, num_face, num_node;
  int *xcnode, *cnode;
  int *xcface, *cface;
  int *xfnode, *fnode;
  int *xncell, *ncell;
  double *coord;
  int *cfdir;
  int *blockid;
  //int *xnodec, *nodec;
} unstr_mesh_data;
}

//using namespace Portage;
using namespace Wonton;

namespace Truchas {

class StateWrapper {
  private:
  double *data_;
  int size_;

  public:
  StateWrapper(int size=0, double *data=NULL) : size_(size), data_(data) {};

  void set_data(int size, double *data) {
    size_ = size;
    data_ = data;
  };

  Entity_kind get_entity(std::string const name) const {
    if (data_) {
      return Entity_kind::CELL;
    } else {
      std::runtime_error("Field data not defined.");
    }
  }

  int get_data_size(const Entity_kind on_what, const std::string name) const {
    if (data_)
      return size_;
    else
      return 0;
  }

  void mesh_get_data(const Entity_kind on_what, const std::string name,
      double **data) {
    (*data) = data_;
  }

  void mesh_get_data(const Entity_kind on_what, const std::string name,
      double const **data) const {
    (*data) = data_;
  }

  template<class T>
  void mat_get_celldata(const std::string& var_name, int matid, T const **data) const {
  }

  template<class T>
  void mat_get_celldata(const std::string& var_name, int matid, T **data) const {
  }

  int cell_index_in_material(int meshcell, int matid) const {
  }

  void mat_get_cells(int matid, std::vector<int> *matcells) const {
  }

  Field_type field_type(Entity_kind on_what, std::string const& var_name) const {
    return Field_type::MESH_FIELD;
  }

  int num_materials() const {return 1;}
};


class MeshWrapper : public AuxMeshTopology<MeshWrapper> {
 private:

  int num_cell, num_face, num_node;
  int *xcnode, *cnode;
  int *xcface, *cface, *cfdir;
  int *blockid;
  int *xfnode, *fnode;
  int *xncell, *ncell;
  double *coord;

 public:

  explicit MeshWrapper(unstr_mesh_data *mesh) :
      AuxMeshTopology<MeshWrapper>(true, true, true) {

    // Copy the mesh data from the interface container
    num_cell  = mesh->num_cell;
    num_face  = mesh->num_face;
    num_node  = mesh->num_node;
    xcnode = mesh->xcnode;
    cnode  = mesh->cnode;
    xcface = mesh->xcface;
    cface  = mesh->cface;
    cfdir  = mesh->cfdir;
    blockid = mesh->blockid;
    xfnode = mesh->xfnode;
    fnode  = mesh->fnode;
    xncell = mesh->xncell;
    ncell  = mesh->ncell;
    coord  = mesh->coord;

    AuxMeshTopology<MeshWrapper>::build_aux_entities();
  }

  // Disabled copy constructor
//  MeshWrapper(MeshWrapper const & inmesh) = delete;

  // Disabled assignment operator
//  MeshWrapper & operator=(MeshWrapper const & inmesh) = delete;

  // Destructor
  ~MeshWrapper() {}

  // The following methods are required by AuxMeshTopology
  // NB: This is currently a serial-only implementation.

  int space_dimension() const { return 3; }

    // Number of OWNED entities
  int num_owned_cells() const { return num_cell; }
  int num_owned_faces() const { return num_face; }
  int num_owned_nodes() const { return num_node; }

  // Number of ghost data entities (SERIAL ONLY)
  int num_ghost_cells() const { return 0; }
  int num_ghost_faces() const { return 0; }
  int num_ghost_nodes() const { return 0; }

  // All cells are OWNED (SERIAL ONLY)
  Entity_type cell_get_type(int const cellid) const {
    //printf("Called into cell_get_type\n");
    return Entity_type::PARALLEL_OWNED;
  }

  // All nodes are OWNED (SERIAL ONLY)
  Entity_type node_get_type(int const nodeid) const {
    printf("Called into node_get_type\n");
    return Entity_type::PARALLEL_OWNED;
  }

  Element_type cell_get_element_type(int const cellid) const {
    printf("Called into cell_get_element_type\n");
    int n = xcnode[cellid+1] - xcnode[cellid]; //FIXME? 0-based cell IDs?
    switch (n) {
      case 4:
        return Element_type::TET;
      case 5:
        return Element_type::PYRAMID;
      case 6:
        return Element_type::PRISM;
      case 8:
        return Element_type::HEX;
      default:
        return Element_type::UNKNOWN_TOPOLOGY;
    }
  }

  // Get cell-to-face connectivity
  void cell_get_faces_and_dirs(int const cellid,
      std::vector<int> *cfaces, std::vector<int> *cfdirs) const {
    //printf("Called into cell_get_faces_and_dirs\n");
    cfaces->clear();
    cfdirs->clear();
    for (int i(xcface[cellid]); i < xcface[cellid+1]; i++) {
      cfaces->push_back(cface[i-1]-1);
      cfdirs->push_back(cfdir[i-1]);
    }
  }

  // Get cell-to-node connectivity
  void cell_get_nodes(int const cellid, std::vector<int> *nodes) const {
    //printf("Called into cell_get_nodes\n");
    nodes->clear();
    for (int i(xcnode[cellid]); i < xcnode[cellid+1]; i++) {
      nodes->push_back(cnode[i-1]-1); //FIXME? 0-based node IDs, cell IDs?
    }
  }

  // Get face-to-node connectivity
  void face_get_nodes(int const faceid, std::vector<int> *nodes) const {
    //printf("Called into face_get_nodes\n");
    nodes->clear();
    for (int i(xfnode[faceid]); i < xfnode[faceid+1]; i++) {
      nodes->push_back(fnode[i-1]-1); //FIXME? 0-based node IDs, face IDs?
    }
  }

  // Get the ID list of cells of a particular (parallel) type attached to a node.
  // SERIAL ONLY -- ignoring the ptype argument
  void node_get_cells(int const nodeid,
                      Entity_type const ptype,
                      std::vector<int> *nodecells) const {
    //printf("Called into node_get_cells\n");
    nodecells->clear();
    for (int i(xncell[nodeid]); i < xncell[nodeid+1]; i++) {
        nodecells->push_back(ncell[i-1]-1);
    }
  }

  // Get the global entity ID (SERIAL ONLY: GLOBAL == LOCAL)
  int get_global_id(int const id, Entity_kind const kind) const {
    printf("Called into get_global_id\n");
    return id;
  }

  template<int D>
  void node_get_coordinates(int const nodeid, Point<D>* pp) const
  {
    assert(D == space_dimension());
  }

};

template<>
void MeshWrapper::node_get_coordinates<3>(int const nodeid, Point<3> *pp) const{
    //printf("Called into node_get_coordinates\n");
  *pp = Point<3>(coord[3*nodeid], coord[3*nodeid+1], coord[3*nodeid+2]);
};

class PortageMapper {

 private:

  MeshWrapper mesh1;
  MeshWrapper mesh2;
  StateWrapper field1;
  StateWrapper field2;
  Portage::CoreDriver<3,Entity_kind::CELL,MeshWrapper,StateWrapper> *driver;

  using PartPair = Portage::PartPair<3, MeshWrapper, StateWrapper>;
  std::vector<PartPair> parts_manager;

  std::vector<std::vector<int>> part_cells1, part_cells2;

  Wonton::vector<std::vector<Portage::Weights_t>> weights;

  // parameters for interpolate_mesh_var
  double lower_bound = -1e99, upper_bound = 1e99;
  Portage::Limiter_type limiter = Portage::DEFAULT_LIMITER;
  Portage::Boundary_Limiter_type bnd_limiter = Portage::BND_NOLIMITER;
  Portage::Empty_fixup_type empty_fixup = Portage::LEAVE_EMPTY;
  double tol = 1e-12;
  int maxitr = 5; // really want 1 (i.e. none?) but give some rope to go bad

 public:

  PortageMapper(unstr_mesh_data *src, unstr_mesh_data *tgt) :
      mesh1(src), mesh2(tgt), field1(StateWrapper()), field2(StateWrapper()),
      driver(new Portage::CoreDriver<3,Entity_kind::CELL,MeshWrapper,StateWrapper>
           (mesh1, field1, mesh2, field2)) {

    // Set default tolerances
    Portage::NumericTolerances_t numtols = Portage::DEFAULT_NUMERIC_TOLERANCES<3>;

    // Tweak min allowable negative intersection volume
    numtols.minimal_intersection_volume = -1.0e-12;  // or -1e-12 or ....
    driver->set_num_tols(numtols);

    // The raw mesh-mesh intersection volumes.
    auto candidates = driver->search<Portage::SearchKDTree>();
    weights = driver->intersect_meshes<Portage::IntersectRnD>(candidates);

    // Generate the set of common block IDs; these are the mapable mesh parts.
    std::set<int> src_bid(src->blockid, src->blockid+src->num_cell); // source block IDs
    std::set<int> tgt_bid(tgt->blockid, tgt->blockid+tgt->num_cell); // target block IDs
    std::set<int> part_bid;
    std::set_intersection(tgt_bid.begin(),tgt_bid.end(),src_bid.begin(),src_bid.end(),
                    std::inserter(part_bid,part_bid.begin()));
    int npart = part_bid.size();

    // Form the list of cell IDs belonging to each part for each mesh. The
    // parts_manager member will hold references to these, so they must persist.
    part_cells1.resize(npart);
    part_cells2.resize(npart);
    int n = 0;
    for (auto &bid: part_bid) {
      for (int j = 0; j < src->num_cell; j++) {
        if (src->blockid[j] == bid) part_cells1[n].push_back(j);
      };
      for (int j = 0; j < tgt->num_cell; j++) {
        if (tgt->blockid[j] == bid) part_cells2[n].push_back(j);
      };
      n++;
    };

    // Initialize mesh part pairs for later part-by-part mapping.
    parts_manager.reserve(npart);
    Wonton::SerialExecutor_type serial_exec;
    // filter cells and populate lists
    for (int i = 0; i < npart; ++i) {
      parts_manager.emplace_back(mesh1, field1, mesh2, field2,
                                 part_cells1[i], part_cells2[i], &serial_exec);
      parts_manager[i].check_mismatch(weights);
    }

  }

  void map_field(int src_size, double *src, int dest_size, double *dest, int method){

    field1.set_data(src_size, src);
    field2.set_data(dest_size, dest);

    auto partial_fixup = Portage::DEFAULT_PARTIAL_FIXUP_TYPE;
    switch (method) {
      case 0: partial_fixup = Portage::LOCALLY_CONSERVATIVE; break;
      case 1: partial_fixup = Portage::CONSTANT; break;
      case 2: partial_fixup = Portage::GLOBALLY_CONSERVATIVE; break;
    }

    for (int i = 0; i < parts_manager.size(); i++) {
      int matid = 0;  // dummy

      // 1st order interpolation the way you had it
      driver->interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
          "foo", "foo", weights,  &(parts_manager[i]));

      /* If you want 2nd order interpolation (You can even supply your
       * own gradients if you have a special method to calculate them)

         auto gradients =
         driver->compute_source_gradient("foo",
                                         limiter, bnd_limiter,
                                         matid,
                                         &(parts_manager[i].source()));

         driver->interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
            "foo", "foo", weights,  &(parts_manager[i]), gradients);
      */

      if (parts_manager[i].has_mismatch())
        parts_manager[i].fix_mismatch("foo", "foo", lower_bound, upper_bound,
                                      tol, maxitr, partial_fixup, empty_fixup);

    }

  }

};

}

#endif // TRUCHAS_PORTAGE_H_

// -----------------------------------------------------------------------------
// This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
// -----------------------------------------------------------------------------
//
// This implements a Fortran-callable interface to a custom Portage data
// remapping driver that is adapted to Truchas mesh data structures. This
// handles part-by-part remapping.
//
// Neil N. Carlson <nnc@lanl.gov>
// November 2019

#include "truchas_portage.hh"

extern "C" {

struct CPortageMapper {
  Truchas::PortageMapper pm;
};

CPortageMapper* portage_mapper_new(unstr_mesh_data *src_mesh, unstr_mesh_data *tgt_mesh)
{
  return new CPortageMapper({Truchas::PortageMapper(src_mesh, tgt_mesh)});
}

void portage_mapper_delete(CPortageMapper *mapper) {
  delete mapper;
}

void portage_map_field(CPortageMapper *self, int src_size, double *src,
                                                int dest_size, double *dest,
                          int method) {
  self->pm.map_field(src_size, src, dest_size, dest, method);
};

}

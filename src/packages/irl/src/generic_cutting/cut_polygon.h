// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_CUT_POLYGON_H_
#define SRC_GENERIC_CUTTING_CUT_POLYGON_H_

#include <algorithm>

#include "src/data_structures/stack_vector.h"
#include "src/generic_cutting/recursive_simplex_cutting/lookup_tables.h"
#include "src/geometry/general/plane.h"
#include "src/geometry/polygons/polygon.h"
#include "src/geometry/polygons/tri.h"
#include "src/geometry/polyhedrons/hexahedron.h"
#include "src/geometry/polyhedrons/rectangular_cuboid.h"
#include "src/helpers/geometric_cutting_helpers.h"
#include "src/parameters/defined_types.h"
#include "src/planar_reconstruction/planar_separator.h"

namespace IRL {

template <class PolygonType, class ReconstructionType>
PolygonType cutPolygonByReconstruction(
    const PolygonType& a_polygon, const Plane* a_plane_that_created_polygon,
    const ReconstructionType& a_reconstruction_to_intersect);

template <class PolygonType, class ReconstructionType>
void cutPolygonByReconstructionInPlace(
    PolygonType* a_polygon, const Plane* a_plane_that_created_polygon,
    const ReconstructionType& a_reconstruction_to_intersect);

/// \brief Cuts a plane by `a_hexahedron`, returning a `PlanePolygon`
/// which is a polygon of the plane contained inside the rectangular_cuboid.
///
/// NOTE: This function does not set the centroid stored in
/// a DividedPolygon when PolygonType is DividedPolygon.
///
/// Template requirements for `HexahedronType`:
/// - Must be a `Hexahedron` or `RectangularCuboid` class object.
///
/// \param[in] a_hexahedron Rectangular cuboid that a_plane
/// will be applied to.
/// \param[in] a_plane Plane that will be cut by `a_hexahedron`.
template <class PolygonType, class HexahedronType>
PolygonType cutPlaneByHexahedron(const HexahedronType& a_hexahedron,
                                 const Plane& a_plane);

template <class PolygonType>
PolygonType cutPlaneByPolyhedron(const RectangularCuboid& a_bouding_box,
                                 const Plane& a_plane);

template <class PolygonType>
PolygonType cutPlaneByPolyhedron(const Hexahedron& a_bouding_box,
                                 const Plane& a_plane);

template <class PolygonType, class HexahedronType>
PolygonType cutPlaneByPolyhedron(const HexahedronType& a_bouding_box,
                                 const PlanarLocalizer& a_planar_localizer,
                                 const Plane& a_plane);

template <class PolygonType>
PolygonType cutPlaneByPolyhedron(const RectangularCuboid& a_bouding_box,
                                 const PlanarLocalizer& a_planar_localizer,
                                 const Plane& a_plane);

template <class PolygonType>
PolygonType cutPlaneByPolyhedron(const Hexahedron& a_bouding_box,
                                 const PlanarLocalizer& a_planar_localizer,
                                 const Plane& a_plane);

/// \brief Loops through the vertices of a polygon and computes a signed
/// distance from `a_cutting_plane` for each. Also sets all remaining
/// distances in array to distance value of first vertex to be cyclic.
///
/// This function computes the signed distance from `a_cutting_plane` to
/// the vertices in `a_polygon` using the `signedDistanceToPoint` method
/// of `a_cutting_plane`. `flip_cut` is needed in order to flip the sign
/// of the distance and cut for the correct phase when the reconstruction
/// is comprised of multiple planes. A supplied array of double[6] is
/// modified in the function.
///
/// \param[in] a_polygon Polygon needing signed distances
///  for its vertices.
/// \param[in] a_cutting_plane The plane that is being
/// intersected with the polygon.
/// \param[in] a_flip_cut Whether the sign should be flipped
///  (indicating cutting to obtain quantities for above
/// the plane instead of below).
/// \param[out] a_distance_to_vertices Array of double[6] to place
///  computed distances for the vertices. to the vertices.
template <class PolygonType>
inline void signedDistanceToPolygonVertices(
    const PolygonType& a_polygon, const Plane& a_cutting_plane,
    const double a_flip_cut, std::array<double, 6>* a_distance_to_vertices);

/// \brief Cuts a `PlanePolygon` by a given plane, returning a new
/// Polygon of type `PolygonType` which is a polygon of the plane after being
/// intersected with `a_cutting_plane`.
///
/// NOTE: This function does not set the centroid stored in
/// a DividedPolygon when PolygonType is DividedPolygon.
///
/// Template requirements for `PolygonType`:
/// - Must be either Polygon or DividedPolygon.
///
/// \param[in] a_polygon Plane that will be cut
/// by `a_cutting_plane`.
/// \param[in] a_cutting_plane Plane that will cut `a_plane_polygon`.
/// \param[in] a_flip_cut Indicates whether below planes (if 1.0) or
/// above planes (if -1.0) polygon should be returned.
template <class PolygonType>
PolygonType cutPolygonByPlane(const PolygonType& a_polygon,
                              const Plane& a_cutting_plane,
                              const double a_flip_cut);

template <class PolygonType>
void cutPolygonByPlaneInPlace(const PolygonType* a_polygon,
                              const Plane& a_cutting_plane,
                              const double a_flip_cut);

/// \brief Return ConvexPolygon that represents the plane in
/// `a_reconstruction[a_plane_to_make_polygon]` on `a_hexahedron.
///
/// NOTE: This function does not set the centroid stored in
/// a DividedPolygon when PolygonType is DividedPolygon.
///
/// Template requirements for `HexahedronType`:
/// - Must be a `Hexahedron` or `RectangularCuboid` class object.
///
/// Template requirements for `PolygonType`:
/// - Must be either Polygon or DividedPolygon.
///
/// \param[in] a_hexahedron Rectangular cuboid `a_reconstruction`
/// will be applied on.
/// \param[in] a_reconstruction Reconstruction to get plane from
/// \param[in] a_plane_to_make_polygon Plane index in `a_reconstruction`
/// that will be made into a polygon.
template <class PolygonType, class PolyhedronType>
PolygonType getPlanePolygonFromReconstruction(
    const PolyhedronType& a_polyhedron, const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon);

template <class PolygonType>
PolygonType getPlanePolygonFromReconstruction(
    const RectangularCuboid& a_polyhedron,
    const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon);

template <class PolygonType>
PolygonType getPlanePolygonFromReconstruction(
    const Hexahedron& a_polyhedron, const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon);

template <class PolygonType, class PolyhedronType>
PolygonType getPlanePolygonFromReconstruction(
    const PolyhedronType& a_polyhedron,
    const PlanarLocalizer& a_polyhedron_localizer,
    const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon);

template <class PolygonType>
PolygonType getPlanePolygonFromReconstruction(
    const RectangularCuboid& a_polyhedron,
    const PlanarLocalizer& a_polyhedron_localizer,
    const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon);

template <class PolygonType>
PolygonType getPlanePolygonFromReconstruction(
    const Hexahedron& a_polyhedron,
    const PlanarLocalizer& a_polyhedron_localizer,
    const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon);

/// \brief Cut a supplied rectangular cuboid by a given reconstruction
///  and return the surface area.
///
/// This function cuts `a_hexahedron` by `a_reconstruction` to obtain
/// the interfacial surface area in that cell. This makes heavy use
/// of large lookup tables held in `lookup_tables.h`.
///
/// NOTE: If `a_reconstruction` contains multiple planes, it will
/// also cut the planes by one another before calculating surface area.
///
/// Template requirements for `HexahedronType`:
/// - Must be a `Hexahedron` or `RectangularCuboid` class object.
///
/// Template requirements for `PolygonType`:
/// - Must be either Polygon or DividedPolygon.
///
/// \param[in] a_hexahedron Rectangular cuboid (such as a
/// computational cell) that the reconstruction will be used for.
/// \param[in] a_reconstruction Planar reconstruction that will be cut
/// by `a_hexahedron` and have its surface area calculated.
template <class HexahedronType>
double getReconstructionSurfaceArea(const HexahedronType& a_hexahedron,
                                    const PlanarSeparator& a_reconstruction);

}  // namespace IRL

#include "src/generic_cutting/cut_polygon.tpp"

#endif  // SRC_GENERIC_CUTTING_CUT_POLYGON_H_

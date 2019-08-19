// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_CUT_POLYGON_TPP_
#define SRC_GENERIC_CUTTING_CUT_POLYGON_TPP_

namespace IRL {

//******************************************************************* //
//     Inlined function definitions placed below this.
//******************************************************************* //
template <class PolygonType>
inline void signedDistanceToPolygonVertices(
    const PolygonType& a_polygon, const Plane& a_cutting_plane,
    const double a_flip_cut, std::array<double, 6>* a_distance_to_vertices) {
  const LookupIndex_t nvert =
      static_cast<LookupIndex_t>(a_polygon.getNumberOfVertices());
  assert(nvert <= 6);
  // Get distance for valid vertices
  for (LookupIndex_t n = 0; n < nvert; ++n) {
    (*a_distance_to_vertices)[n] =
        a_flip_cut * a_cutting_plane.signedDistanceToPoint(a_polygon[n]);
  }
  // Get distance for wrap and carry same distance to rest of entries
  // Only have lookup tables to cut polygons with <=6 vertices
  for (LookupIndex_t n = nvert; n < 6; ++n) {
    (*a_distance_to_vertices)[n] = (*a_distance_to_vertices)[0];
  }
}

//******************************************************************* //
//     Template function definitions placed below this.
//******************************************************************* //
template <class PolygonType, class ReconstructionType>
PolygonType cutPolygonByReconstruction(
    const PolygonType& a_polygon, const Plane* a_plane_that_created_polygon,
    const ReconstructionType& a_reconstruction_to_intersect) {
  PolygonType plane_polygon_to_return = a_polygon;
  cutPolygonByReconstructionInPlace(&plane_polygon_to_return,
                                    a_plane_that_created_polygon,
                                    a_reconstruction_to_intersect);
  plane_polygon_to_return.setPlaneOfExistence(a_plane_that_created_polygon);
  return plane_polygon_to_return;
}

template <class PolygonType, class ReconstructionType>
void cutPolygonByReconstructionInPlace(
    PolygonType* a_polygon, const Plane* a_plane_that_created_polygon,
    const ReconstructionType& a_reconstruction_to_intersect) {
  for (const auto& other_plane : a_reconstruction_to_intersect) {
    // Cut by all other planes in this reconstruction
    if (a_polygon->getNumberOfVertices() == 0) {
      return;
    }
    if (a_plane_that_created_polygon != &other_plane) {
      cutPolygonByPlaneInPlace<PolygonType>(
          a_polygon, other_plane, a_reconstruction_to_intersect.flip());
    }
  }
}

template <class PolygonType, class HexahedronType>
PolygonType cutPlaneByHexahedron(const HexahedronType& a_hexahedron,
                                 const Plane& a_plane) {
  using PtType = typename HexahedronType::pt_type;
  std::array<double, HexahedronType::number_of_vertices> vertex_distance;
  signedDistanceToVertices(a_hexahedron, a_plane, &vertex_distance);
  LookupIndex_t cutting_case = getGeometricCaseId(vertex_distance);

  PolygonType plane_polygon_to_return;
  plane_polygon_to_return.setNumberOfVertices(
      cut_plane_by_cuboid::number_of_plane_intersections_with_cuboid
          [cutting_case]);
  for (LookupIndex_t n = 0;
       n < cut_plane_by_cuboid::number_of_plane_intersections_with_cuboid
               [cutting_case];
       ++n) {
    LookupIndex_t v1 =
        cut_plane_by_cuboid::cuboid_cut_vertices[cutting_case][0][n];
    LookupIndex_t v2 =
        cut_plane_by_cuboid::cuboid_cut_vertices[cutting_case][1][n];
    plane_polygon_to_return[n] = a_hexahedron[v1].fromEdgeIntersection(
        a_hexahedron[v1], vertex_distance[v1], a_hexahedron[v2],
        vertex_distance[v2]);
  }
  plane_polygon_to_return.setPlaneOfExistence(a_plane);
  return plane_polygon_to_return;
}

template <class PolygonType>
PolygonType cutPlaneByPolyhedron(const RectangularCuboid& a_bounding_box,
                                 const PlanarLocalizer& a_planar_localizer,
                                 const Plane& a_plane) {
  auto polygon_to_return =
      cutPlaneByHexahedron<PolygonType>(a_bounding_box, a_plane);
  return polygon_to_return;
}

template <class PolygonType>
PolygonType cutPlaneByPolyhedron(const Hexahedron& a_bounding_box,
                                 const PlanarLocalizer& a_planar_localizer,
                                 const Plane& a_plane) {
  auto polygon_to_return =
      cutPlaneByHexahedron<PolygonType>(a_bounding_box, a_plane);
  return polygon_to_return;
}

template <class PolygonType, class PolyhedronType>
PolygonType cutPlaneByPolyhedron(const PolyhedronType& a_bounding_box,
                                 const PlanarLocalizer& a_planar_localizer,
                                 const Plane& a_plane) {
  auto polygon_to_return =
      cutPlaneByHexahedron<PolygonType>(a_bounding_box, a_plane);
  cutPolygonByReconstructionInPlace(&polygon_to_return, nullptr,
                                    a_planar_localizer);
  return polygon_to_return;
}

template <class PolygonType>
PolygonType cutPlaneByPolyhedron(const RectangularCuboid& a_bounding_box,
                                 const Plane& a_plane) {
  auto polygon_to_return =
      cutPlaneByHexahedron<PolygonType>(a_bounding_box, a_plane);
  return polygon_to_return;
}

template <class PolygonType>
PolygonType cutPlaneByPolyhedron(const Hexahedron& a_bounding_box,
                                 const Plane& a_plane) {
  auto polygon_to_return =
      cutPlaneByHexahedron<PolygonType>(a_bounding_box, a_plane);
  return polygon_to_return;
}

template <class PolygonType>
PolygonType cutPolygonByPlane(const PolygonType& a_polygon,
                              const Plane& a_cutting_plane,
                              const double a_flip_cut) {
  PolygonType polygon_to_return = a_polygon;
  cutPolygonByPlaneInPlace(&polygon_to_return, a_cutting_plane, a_flip_cut);
  polygon_to_return.setPlaneOfExistence(a_cutting_plane);
  return polygon_to_return;
}

template <class PolygonType>
double getIntersectionPts(
    const PolygonType& a_polygon, const Plane& a_cutting_plane,
    const double a_flip_cut, StackVector<Pt, 2>* a_intersection_pts,
    StackVector<UnsignedIndex_t, 2>* a_index_after_intersections) {
  double distance_to_vertex =
      a_flip_cut * a_cutting_plane.signedDistanceToPoint(a_polygon[0]);
  // Assumed here that polygon wraps around on itself, with a_polygon[0] ==
  // a_polygon[a_polygon.size()-1]
  for (UnsignedIndex_t n = 0; n < a_polygon.getNumberOfVertices() - 1; ++n) {
    double distance_to_next_vertex =
        a_flip_cut * a_cutting_plane.signedDistanceToPoint(a_polygon[n + 1]);
    if (distance_to_vertex * distance_to_next_vertex < 0.0) {
      a_intersection_pts->push_back(
          Pt::fromEdgeIntersection(a_polygon[n], distance_to_vertex,
                                   a_polygon[n + 1], distance_to_next_vertex));
      a_index_after_intersections->push_back(n + 1);
      if (a_index_after_intersections->size() == 2) {
        return distance_to_next_vertex;
      }
    }
    distance_to_vertex = distance_to_next_vertex;
  }
  return distance_to_vertex;
}

template <class PolygonType>
void cutPolygonByPlaneInPlace(PolygonType* a_polygon,
                              const Plane& a_cutting_plane,
                              const double a_flip_cut) {
  assert(a_polygon->getNumberOfVertices() >= 3);
  const UnsignedIndex_t starting_polygon_size =
      a_polygon->getNumberOfVertices();
  StackVector<Pt, 2> intersection_pts;
  StackVector<UnsignedIndex_t, 2> index_after_intersection;
  a_polygon->addVertex((*a_polygon)[0]);
  double distance_to_vertex_after_second_intersection =
      getIntersectionPts(*a_polygon, a_cutting_plane, a_flip_cut,
                         &intersection_pts, &index_after_intersection);
  if (intersection_pts.size() == 0) {
    // No intersection with plane, distance to any point will tell us if fully
    // under or over plane
    if (distance_to_vertex_after_second_intersection <= 0.0) {
      // Just removing the vertex we added at start
      a_polygon->removeLastVertex();
      return;
    } else {
      a_polygon->setNumberOfVertices(0);
      return;
    }
  } else {  // have intersections
    if (distance_to_vertex_after_second_intersection < 0.0) {
      // Second intersection went from above plane to below plane
      UnsignedIndex_t new_size =
          intersection_pts.size() + index_after_intersection[0] +
          (starting_polygon_size - index_after_intersection[1]);
      std::copy_backward(a_polygon->begin() + index_after_intersection[1],
                         a_polygon->end() - 1, a_polygon->begin() + new_size);
      // Need to copy_backward before adding in intersections
      // because copying of locations later overwritten might be done.
      (*a_polygon)[index_after_intersection[0]] = intersection_pts[0];
      (*a_polygon)[index_after_intersection[0] + 1] = intersection_pts[1];
      a_polygon->setNumberOfVertices(new_size);
      return;
    } else {
      // Second intersection went from below plane to above plane
      UnsignedIndex_t new_size = intersection_pts.size() +
                                 index_after_intersection[1] -
                                 index_after_intersection[0];
      (*a_polygon)[0] = intersection_pts[0];
      std::copy(a_polygon->begin() + index_after_intersection[0],
                a_polygon->begin() + index_after_intersection[1],
                a_polygon->begin() + 1);
      (*a_polygon)[new_size - 1] = intersection_pts[1];
      a_polygon->setNumberOfVertices(new_size);
      return;
    }
  }
}

template <class PolygonType, class PolyhedronType>
PolygonType getPlanePolygonFromReconstruction(
    const PolyhedronType& a_polyhedron, const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon) {
  // Get plane that exists inside cube
  const auto polyhedron_localizer = a_polyhedron.getLocalizer();
  return getPlanePolygonFromReconstruction<PolygonType>(a_polyhedron, polyhedron_localizer,
                                           a_reconstruction,
                                           a_plane_to_make_polygon);
}

template <class PolygonType>
PolygonType getPlanePolygonFromReconstruction(
    const RectangularCuboid& a_polyhedron,
    const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon) {
  // Get plane that exists inside cube
  auto plane_polygon_to_return =
      cutPlaneByPolyhedron<PolygonType>(a_polyhedron, a_plane_to_make_polygon);
  cutPolygonByReconstructionInPlace(&plane_polygon_to_return,
                                    &a_plane_to_make_polygon, a_reconstruction);
  return plane_polygon_to_return;
}

template <class PolygonType>
PolygonType getPlanePolygonFromReconstruction(
    const Hexahedron& a_polyhedron, const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon) {
  // Get plane that exists inside cube
  auto plane_polygon_to_return =
      cutPlaneByPolyhedron<PolygonType>(a_polyhedron, a_plane_to_make_polygon);
  cutPolygonByReconstructionInPlace(&plane_polygon_to_return,
                                    &a_plane_to_make_polygon, a_reconstruction);
  return plane_polygon_to_return;
}

template <class PolygonType, class PolyhedronType>
PolygonType getPlanePolygonFromReconstruction(
    const PolyhedronType& a_polyhedron,
    const PlanarLocalizer& a_polyhedron_localizer,
    const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon) {
  auto lower_pt = a_polyhedron.getLowerLimits();
  auto upper_pt = a_polyhedron.getUpperLimits();
  auto plane_polygon_to_return = cutPlaneByPolyhedron<PolygonType>(
      RectangularCuboid::fromBoundingPts(lower_pt, upper_pt),
      a_polyhedron_localizer, a_plane_to_make_polygon);
  cutPolygonByReconstructionInPlace(&plane_polygon_to_return,
                                    &a_plane_to_make_polygon, a_reconstruction);
  return plane_polygon_to_return;
}

template <class PolygonType>
PolygonType getPlanePolygonFromReconstruction(
    const RectangularCuboid& a_polyhedron,
    const PlanarLocalizer& a_polyhedron_localizer,
    const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon) {
  return getPlanePolygonFromReconstruction<PolygonType>(
      a_polyhedron, a_reconstruction, a_plane_to_make_polygon);
}

template <class PolygonType>
PolygonType getPlanePolygonFromReconstruction(
    const Hexahedron& a_polyhedron,
    const PlanarLocalizer& a_polyhedron_localizer,
    const PlanarSeparator& a_reconstruction,
    const Plane& a_plane_to_make_polygon) {
  return getPlanePolygonFromReconstruction<PolygonType>(
      a_polyhedron, a_reconstruction, a_plane_to_make_polygon);
}

template <class HexahedronType>
double getReconstructionSurfaceArea(const HexahedronType& a_hexahedron,
                                    const PlanarSeparator& a_reconstruction) {
  double surface_area = {0.0};
  for (const auto& plane : a_reconstruction) {
    auto plane_polygon = getPlanePolygonFromReconstruction<Polygon>(
        a_hexahedron, a_reconstruction, plane);
    // Calculate and accumulate surface area
    surface_area += plane_polygon.calculateConvexVolume();
  }
  return surface_area;
}

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_CUT_POLYGON_TPP_

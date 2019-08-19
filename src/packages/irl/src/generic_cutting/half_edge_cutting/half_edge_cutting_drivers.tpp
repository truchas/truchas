// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_DRIVERS_TPP_
#define SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_DRIVERS_TPP_

#include "src/generic_cutting/general/encountered_id_list.h"

namespace IRL {

// Foward declare getVolumeMoments to avoid circular dependency.
template <class ReturnType, class CuttingMethod, class SegmentedPolytopeType,
          class HalfEdgePolytopeType, class ReconstructionType>
__attribute__((hot)) inline ReturnType getVolumeMoments(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction);

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
enable_if_t<HasAReconstructionLink<ReconstructionType>::value>
getVolumeMomentsForPolytope(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_originating_reconstruction,
    const ReconstructionType& a_reconstruction, EncounteredIdList* a_id_list,
    ReturnType* a_moments_to_return);

namespace getVolumeMomentsForPolytopeDetails {
template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType, typename Enable = void>
struct getVolumeMomentsForPolytopeStaticStructWrapper;

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_reconstruction,
      ReturnType* a_moments_to_return);
};

// For ReconstructionLink cutting with informed origin.
template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType, typename Enable = void>
struct getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper;

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_originating_reconstruction,
      const ReconstructionType& a_reconstruction, EncounteredIdList* a_id_list,
      ReturnType* a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_originating_reconstruction,
      const ReconstructionType& a_reconstruction, EncounteredIdList* a_id_list,
      ReturnType* a_moments_to_return);
};

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
struct getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>> {
  static void getVolumeMomentsForPolytopeImplementation(
      SegmentedPolytopeType* a_polytope,
      HalfEdgePolytopeType* a_complete_polytope,
      const ReconstructionType& a_originating_reconstruction,
      const ReconstructionType& a_reconstruction, EncounteredIdList* a_id_list,
      ReturnType* a_moments_to_return);
};

}  // namespace getVolumeMomentsForPolytopeDetails

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType>
inline void localizeInternalToReconstruction(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction);

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
inline void splitAndShareThroughLinks(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction, EncounteredIdList* a_id_list,
    ReturnType* a_moments_to_return);

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void splitAndShareThroughLinks(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_originating_reconstruction,
    const ReconstructionType& a_reconstruction, EncounteredIdList* a_id_list,
    ReturnType* a_moments_to_return);

//******************************************************************* //
//     Function template definitions placed below this
//******************************************************************* //

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytope(SegmentedPolytopeType* a_polytope,
                                 HalfEdgePolytopeType* a_complete_polytope,
                                 const ReconstructionType& a_reconstruction,
                                 ReturnType* a_moments_to_return) {
  getVolumeMomentsForPolytopeDetails::
      getVolumeMomentsForPolytopeStaticStructWrapper<
          SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType,
          ReturnType>::
          getVolumeMomentsForPolytopeImplementation(
              a_polytope, a_complete_polytope, a_reconstruction,
              a_moments_to_return);
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
enable_if_t<HasAReconstructionLink<ReconstructionType>::value>
getVolumeMomentsForPolytope(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_originating_reconstruction,
    const ReconstructionType& a_reconstruction, EncounteredIdList* a_id_list,
    ReturnType* a_moments_to_return) {
  getVolumeMomentsForPolytopeDetails::
      getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
          SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType,
          ReturnType>::
          getVolumeMomentsForPolytopeImplementation(
              a_polytope, a_complete_polytope, a_originating_reconstruction,
              a_reconstruction, a_id_list, a_moments_to_return);
}

namespace getVolumeMomentsForPolytopeDetails {

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  localizeInternalToReconstruction(a_polytope, a_complete_polytope,
                                   a_reconstruction);
  if (a_polytope->getNumberOfFaces() > 0) {
    *a_moments_to_return += IRL::getVolumeMoments<ReturnType, HalfEdgeCutting>(
        a_polytope, a_complete_polytope,
        a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  EncounteredIdList id_list;
  splitAndShareThroughLinks(a_polytope, a_complete_polytope, a_reconstruction,
                            &id_list, a_moments_to_return);
  if (a_polytope->getNumberOfFaces() > 0) {
    *a_moments_to_return += IRL::getVolumeMoments<ReturnType, HalfEdgeCutting>(
        a_polytope, a_complete_polytope,
        a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<DoesNotHaveAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  localizeInternalToReconstruction(a_polytope, a_complete_polytope,
                                   a_reconstruction);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  if (a_polytope->getNumberOfFaces() > 0) {
    *a_moments_to_return +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, HalfEdgeCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  EncounteredIdList id_list;
  splitAndShareThroughLinks(a_polytope, a_complete_polytope, a_reconstruction,
                            &id_list, a_moments_to_return);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;

  if (a_polytope->getNumberOfFaces() > 0) {
    *a_moments_to_return +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, HalfEdgeCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_reconstruction,
        ReturnType* a_moments_to_return) {
  EncounteredIdList id_list;
  splitAndShareThroughLinks(a_polytope, a_complete_polytope, a_reconstruction,
                            &id_list, a_moments_to_return);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;

  if (a_polytope->getNumberOfFaces() > 0) {
    assert(a_reconstruction.isIdSet());
    (*a_moments_to_return)[a_reconstruction.getId()] +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, HalfEdgeCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                DoesNotHaveANestedType<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_originating_reconstruction,
        const ReconstructionType& a_reconstruction,
        EncounteredIdList* a_id_list, ReturnType* a_moments_to_return) {
  splitAndShareThroughLinks(a_polytope, a_complete_polytope,
                            a_originating_reconstruction, a_reconstruction,
                            a_id_list, a_moments_to_return);
  if (a_polytope->getNumberOfFaces() > 0) {
    *a_moments_to_return += IRL::getVolumeMoments<ReturnType, HalfEdgeCutting>(
        a_polytope, a_complete_polytope,
        a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasANestedType<ReturnType>::value &&
                DoesNotHaveACollection<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_originating_reconstruction,
        const ReconstructionType& a_reconstruction,
        EncounteredIdList* a_id_list, ReturnType* a_moments_to_return) {
  splitAndShareThroughLinks(a_polytope, a_complete_polytope,
                            a_originating_reconstruction, a_reconstruction,
                            a_id_list, a_moments_to_return);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  if (a_polytope->getNumberOfFaces() > 0) {
    *a_moments_to_return +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, HalfEdgeCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void getVolumeMomentsForPolytopeKnownOriginStaticStructWrapper<
    SegmentedPolytopeType, HalfEdgePolytopeType, ReconstructionType, ReturnType,
    enable_if_t<HasAReconstructionLink<ReconstructionType>::value &&
                HasACollection<ReturnType>::value>>::
    getVolumeMomentsForPolytopeImplementation(
        SegmentedPolytopeType* a_polytope,
        HalfEdgePolytopeType* a_complete_polytope,
        const ReconstructionType& a_originating_reconstruction,
        const ReconstructionType& a_reconstruction,
        EncounteredIdList* a_id_list, ReturnType* a_moments_to_return) {
  splitAndShareThroughLinks(a_polytope, a_complete_polytope,
                            a_originating_reconstruction, a_reconstruction,
                            a_id_list, a_moments_to_return);
  using WantedVolumeMomentsType = typename ReturnType::contained_type;
  if (a_polytope->getNumberOfFaces() > 0) {
    assert(a_reconstruction.isIdSet());
    (*a_moments_to_return)[a_reconstruction.getId()] +=
        IRL::getVolumeMoments<WantedVolumeMomentsType, HalfEdgeCutting>(
            a_polytope, a_complete_polytope,
            a_reconstruction.getNextReconstruction());
  }
}

}  // namespace getVolumeMomentsForPolytopeDetails

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType>
void localizeInternalToReconstruction(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_reconstruction) {
  const auto& cutting_reconstruction =
      a_reconstruction.getCurrentReconstruction();

  for (const auto& plane : cutting_reconstruction) {
    const auto cutting_plane = cutting_reconstruction.isFlipped()
                                   ? plane.generateFlippedPlane()
                                   : plane;
    if (a_polytope->getNumberOfFaces() == 0) {
      return;
    }
    truncateHalfEdgePolytope(a_polytope, a_complete_polytope, cutting_plane);
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void splitAndShareThroughLinks(SegmentedPolytopeType* a_polytope,
                               HalfEdgePolytopeType* a_complete_polytope,
                               const ReconstructionType& a_reconstruction,
                               EncounteredIdList* a_id_list,
                               ReturnType* a_moments_to_return) {
  const auto& cutting_reconstruction =
      a_reconstruction.getCurrentReconstruction();

  a_id_list->addId(a_reconstruction.getId());
  const auto starting_id_size = a_id_list->size();

  SegmentedPolytopeType clipped_polytope;
  for (UnsignedIndex_t plane_index = 0;
       plane_index < cutting_reconstruction.getNumberOfPlanes();
       ++plane_index) {

    if (a_polytope->getNumberOfFaces() == 0) {
	  return;
	}

    if (a_reconstruction.hasNeighbor(plane_index)) {
      if (a_id_list->isIdPresent(
              a_reconstruction.getNeighbor(plane_index).getId())) {
        continue;
      }
    }

    const auto cutting_plane =
        cutting_reconstruction.isFlipped()
            ? cutting_reconstruction[plane_index].generateFlippedPlane()
            : cutting_reconstruction[plane_index];

    if (a_reconstruction.hasNeighbor(plane_index)) {
      splitHalfEdgePolytope(a_polytope, &clipped_polytope, a_complete_polytope,
                            cutting_plane);
      getVolumeMomentsForPolytope(&clipped_polytope, a_complete_polytope,
                                  a_reconstruction,
                                  a_reconstruction.getNeighbor(plane_index),
                                  a_id_list, a_moments_to_return);
      //clipped_polytope.clear(a_complete_polytope);
      a_id_list->resize(starting_id_size);
    } else {
      truncateHalfEdgePolytope(a_polytope, a_complete_polytope, cutting_plane);
    }
  }
}

template <class SegmentedPolytopeType, class HalfEdgePolytopeType,
          class ReconstructionType, class ReturnType>
void splitAndShareThroughLinks(
    SegmentedPolytopeType* a_polytope,
    HalfEdgePolytopeType* a_complete_polytope,
    const ReconstructionType& a_originating_reconstruction,
    const ReconstructionType& a_reconstruction, EncounteredIdList* a_id_list,
    ReturnType* a_moments_to_return) {
  const auto& cutting_reconstruction =
      a_reconstruction.getCurrentReconstruction();

  a_id_list->addId(a_reconstruction.getId());
  const auto starting_id_size = a_id_list->size();

  SegmentedPolytopeType clipped_polytope;
  for (UnsignedIndex_t plane_index = 0;
       plane_index < cutting_reconstruction.getNumberOfPlanes();
       ++plane_index) {
    if (a_polytope->getNumberOfFaces() == 0) {
      return;
    }

    if (a_reconstruction.hasNeighbor(plane_index)) {
      if (a_id_list->isIdPresent(
              a_reconstruction.getNeighbor(plane_index).getId())) {
        continue;
      }
    }
    const auto cutting_plane =
        cutting_reconstruction.isFlipped()
            ? cutting_reconstruction[plane_index].generateFlippedPlane()
            : cutting_reconstruction[plane_index];

    if (a_reconstruction.hasNeighbor(plane_index)) {
      splitHalfEdgePolytope(a_polytope, &clipped_polytope, a_complete_polytope,
                            cutting_plane);
      getVolumeMomentsForPolytope(&clipped_polytope, a_complete_polytope,
                                  a_reconstruction,
                                  a_reconstruction.getNeighbor(plane_index),
                                  a_id_list, a_moments_to_return);
      //clipped_polytope.clear(a_complete_polytope);      
      a_id_list->resize(starting_id_size);
    } else {
      truncateHalfEdgePolytope(a_polytope, a_complete_polytope, cutting_plane);
    }
  }
}

}  // namespace IRL

#endif  // SRC_GENERIC_CUTTING_HALF_EDGE_CUTTING_HALF_EDGE_CUTTING_DRIVERS_TPP_

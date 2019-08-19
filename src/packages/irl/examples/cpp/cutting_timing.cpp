// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "examples/cpp/cutting_timing.h"

#include <chrono>
#include <random>
#include <utility>

#include "src/generic_cutting/generic_cutting.h"
#include "src/geometry/polyhedrons/tet.h"

int main(void) {
  static constexpr IRL::UnsignedIndex_t ncycles = 100000; //1000000;
  std::cout << "Timing comparison : " << std::endl;

  auto time_by_cutting_into_tets = cutWithRecursiveSimplex(ncycles);
  std::cout << "Recursive cutting : " << time_by_cutting_into_tets << std::endl;

  auto time_by_cutting_polyhedron = cutWithHalfEdge(ncycles);
  std::cout << "Half-Edge cutting : " << time_by_cutting_polyhedron
            << std::endl;

  auto time_by_cutting_simplex = cutWithSimplex(ncycles);
  std::cout << "Simplex cutting : " << time_by_cutting_simplex << std::endl;
}

static double randomDoubleFromNormalDist(void){
  auto u1 = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
  auto u2 = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
  return std::sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);
}
IRL::Tet generateRandomTet(const double minimum_volume) {
  IRL::Tet tet_to_return;
  double tet_vol = 0.0;
  do {
    for (auto& a_vertex : tet_to_return) {
      IRL::Normal tet_vertex_stand_in;
      for (auto& component : tet_vertex_stand_in) {
        component = randomDoubleFromNormalDist();
      }
      tet_vertex_stand_in.normalize();
      a_vertex = IRL::Normal::toPt(tet_vertex_stand_in);
    }
    tet_vol = tet_to_return.calculateVolume();
    if (tet_vol < -minimum_volume) {
      tet_vol = -tet_vol;
      std::swap(tet_to_return[2], tet_to_return[3]);
    }
  } while (tet_vol < minimum_volume);
  return tet_to_return;
}

double cutWithRecursiveSimplex(
    const IRL::UnsignedIndex_t number_of_iterations) {
  std::chrono::duration<double> duration(0.0);
  auto start = std::chrono::system_clock::now();
  double old_minimum_volume = IRL::global_constants::MINIMUM_VOLUME_TO_TRACK;
  constexpr double minimum_volume = 1.0e-8;
  IRL::setMinimumVolumeToTrack(minimum_volume);
  for (IRL::UnsignedIndex_t cycle = 0; cycle < number_of_iterations; ++cycle) {
    IRL::Tet tet_for_volume = generateRandomTet(minimum_volume);
    IRL::PlanarLocalizer tet_localizer =
        generateRandomTet(minimum_volume).getLocalizer();
    auto volume_moments =
        IRL::getNormalizedVolumeMoments<IRL::VolumeMoments,
                                        IRL::RecursiveSimplexCutting>(
            tet_for_volume, tet_localizer);
  }
  IRL::setMinimumVolumeToTrack(old_minimum_volume);
  auto end = std::chrono::system_clock::now();
  duration = end-start;
  return duration.count();
}

double cutWithSimplex(const IRL::UnsignedIndex_t number_of_iterations) {
  std::chrono::duration<double> duration(0.0);
  auto start = std::chrono::system_clock::now();
  double old_minimum_volume = IRL::global_constants::MINIMUM_VOLUME_TO_TRACK;
  constexpr double minimum_volume = 1.0e-8;
  IRL::setMinimumVolumeToTrack(minimum_volume);
  for (IRL::UnsignedIndex_t cycle = 0; cycle < number_of_iterations; ++cycle) {
    IRL::Tet tet_for_volume = generateRandomTet(minimum_volume);
    IRL::PlanarLocalizer tet_localizer =
        generateRandomTet(minimum_volume).getLocalizer();
    auto volume_moments = IRL::getNormalizedVolumeMoments<IRL::VolumeMoments,
                                                          IRL::SimplexCutting>(
        tet_for_volume, tet_localizer);
  }
  IRL::setMinimumVolumeToTrack(old_minimum_volume);
  auto end = std::chrono::system_clock::now();
  duration = end-start;
  return duration.count();
}

double cutWithHalfEdge(const IRL::UnsignedIndex_t number_of_iterations) {
  std::chrono::duration<double> duration(0.0);
  auto start = std::chrono::system_clock::now();
  double old_minimum_volume = IRL::global_constants::MINIMUM_VOLUME_TO_TRACK;
  constexpr double minimum_volume = 1.0e-8;
  IRL::setMinimumVolumeToTrack(minimum_volume);
  for (IRL::UnsignedIndex_t cycle = 0; cycle < number_of_iterations; ++cycle) {
    IRL::Tet tet_for_volume = generateRandomTet(minimum_volume);
    IRL::PlanarLocalizer tet_localizer =
        generateRandomTet(minimum_volume).getLocalizer();
    auto volume_moments = IRL::getNormalizedVolumeMoments<IRL::VolumeMoments,
                                                          IRL::HalfEdgeCutting>(
        tet_for_volume, tet_localizer);
  }
  IRL::setMinimumVolumeToTrack(old_minimum_volume);
  auto end = std::chrono::system_clock::now();
  duration = end-start;
  return duration.count();
}

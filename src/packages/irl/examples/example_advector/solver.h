// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_EXAMPLE_ADVECTOR_SOLVER_H_
#define EXAMPLES_EXAMPLE_ADVECTOR_SOLVER_H_

#include <chrono>
#include <iostream>

#include "examples/example_advector/basic_mesh.h"
#include "examples/example_advector/data.h"
#include "examples/example_advector/vof_advection.h"
#include "src/planar_reconstruction/localized_separator_link.h"
#include "src/planar_reconstruction/planar_localizer.h"
#include "src/planar_reconstruction/planar_separator.h"

/// \brief Handles running and advancing the solution according to provided
/// static functions in structs.
template <class SimulationType, class VOFAdvector, class InterfaceReconstructor>
int runSimulation(const double a_dt, const double a_end_time,
                  const int a_visualization_frequency);

// \brief Convert and store the mesh cells into localizers.
void initializeLocalizers(Data<IRL::PlanarLocalizer>* a_localizers);

/// \brief Initialize the linked localized separators used during advection.
void initializeLocalizedSeparators(
    const Data<IRL::PlanarLocalizer>& a_cell_localizers,
    const Data<IRL::PlanarSeparator>& a_interface,
    Data<IRL::LocalizedSeparatorLink>* a_linked_localized_separators);

/// \brief Set phase quantities according to the given
void setPhaseQuantities(const Data<IRL::PlanarSeparator>& a_interface,
                        Data<double>* a_liquid_volume_fraction,
                        Data<IRL::Pt>* a_liquid_centroid,
                        Data<IRL::Pt>* a_gas_centroid);

/// \brief Write out the header for the diagnostics.
void writeDiagnosticsHeader(void);

/// \brief Write out mesh to file to be plotted during visualization.
void writeOutMesh(const BasicMesh& a_mesh);

/// \brief Write out diagnostics to the screen.
void writeOutDiagnostics(const int a_iteration, const double a_dt,
                         const double a_simulation_time,
                         const Data<double>& a_U, const Data<double>& a_V,
                         const Data<double>& a_W,
                         const Data<double>& a_liquid_volume_fraction,
                         const Data<IRL::PlanarSeparator>& a_interface,
                         std::chrono::duration<double> a_VOF_duration,
                         std::chrono::duration<double> a_recon_duration);

/// \brief Write out visualization files for python.
void writeOutVisualization(const int a_iteration,
                           const int a_visualization_frequency,
                           const double a_simulation_time,
                           const Data<double>& a_liquid_volume_fraction);

//******************************************************************* //
//     Template function definitions placed below this.
//******************************************************************* //
template <class SimulationType, class VOFAdvector, class InterfaceReconstructor>
int runSimulation(const double a_dt, const double a_end_time,
                  const int a_visualization_frequency) {
  // Set mesh
  BasicMesh cc_mesh = SimulationType::setMesh();

  // Allocate local data
  Data<double> velU(&cc_mesh);
  Data<double> velV(&cc_mesh);
  Data<double> velW(&cc_mesh);
  Data<double> liquid_volume_fraction(&cc_mesh);
  Data<IRL::Pt> liquid_centroid(&cc_mesh);
  Data<IRL::Pt> gas_centroid(&cc_mesh);

  // Allocate
  Data<IRL::PlanarLocalizer> cell_localizers(&cc_mesh);
  initializeLocalizers(&cell_localizers);
  Data<IRL::PlanarSeparator> interface(&cc_mesh);
  Data<IRL::LocalizedSeparatorLink> link_localized_separators(&cc_mesh);
  initializeLocalizedSeparators(cell_localizers, interface,
                                &link_localized_separators);
  connectMesh(cc_mesh, &link_localized_separators);
  // Set constants in IRL
  IRL::setMinimumVolumeToTrack(10.0 * DBL_EPSILON * cc_mesh.dx() *
                               cc_mesh.dy() * cc_mesh.dz());
  IRL::setVolumeFractionBounds(1.0e-8);
  IRL::setVolumeFractionTolerance(1.0e-13);

  // Initialize data
  SimulationType::initialize(&velU, &velV, &velW, &interface);
  setPhaseQuantities(interface, &liquid_volume_fraction, &liquid_centroid,
                     &gas_centroid);

  double simulation_time = 0.0;
  int iteration = 0;
  writeDiagnosticsHeader();
  writeOutMesh(cc_mesh);
  std::chrono::duration<double> advect_VOF_time(0.0);
  std::chrono::duration<double> recon_time(0.0);
  writeOutDiagnostics(iteration, a_dt, simulation_time, velU, velV, velW,
                      liquid_volume_fraction, interface, advect_VOF_time,
                      recon_time);
  writeOutVisualization(iteration, a_visualization_frequency, simulation_time,
                        liquid_volume_fraction);
  while (simulation_time < a_end_time) {
    SimulationType::setVelocity(simulation_time, &velU, &velV, &velW);

    auto start = std::chrono::system_clock::now();
    VOFAdvector::advectVOF(a_dt, velU, velV, velW, &link_localized_separators,
                           &liquid_volume_fraction, &liquid_centroid,
                           &gas_centroid);
    auto advect_end = std::chrono::system_clock::now();
    advect_VOF_time = advect_end - start;
    InterfaceReconstructor::getReconstruction(
        liquid_volume_fraction, liquid_centroid, gas_centroid,
        link_localized_separators, a_dt, velU, velV, velW, &interface);
    auto recon_end = std::chrono::system_clock::now();
    recon_time = recon_end - advect_end;

    ++iteration;
    simulation_time += a_dt;
    writeOutDiagnostics(iteration, a_dt, simulation_time, velU, velV, velW,
                        liquid_volume_fraction, interface, advect_VOF_time,
                        recon_time);

    if (iteration % a_visualization_frequency == 0) {
      writeOutVisualization(iteration, a_visualization_frequency,
                            simulation_time, liquid_volume_fraction);
    }
  }

  return 0;
}

#endif  // EXAMPLES_EXAMPLE_ADVECTOR_SOLVER_H_

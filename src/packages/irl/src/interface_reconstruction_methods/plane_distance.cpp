// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/interface_reconstruction_methods/plane_distance.h"

#include "src/helpers/helper.h"
#include "src/parameters/constants.h"

namespace IRL {

double findDistanceOnePlane(const RectangularCuboid& a_rectangular_cuboid,
                            const double a_volume_fraction,
                            const Normal& a_normal) {
  double mm[3] = {a_normal[0] * a_rectangular_cuboid.calculateSideLength(0),
                  a_normal[1] * a_rectangular_cuboid.calculateSideLength(1),
                  a_normal[2] * a_rectangular_cuboid.calculateSideLength(2)};
  double norm = std::fabs(mm[0]) + std::fabs(mm[1]) + std::fabs(mm[2]);
  double factor = {0.0};
  for (int d = 0; d < 3; ++d) {
    mm[d] /= norm;
    if (mm[d] < 0.0) {
      factor += mm[d];
    }
  }
  double VOFo =
      a_volume_fraction > 0.5 ? 1.0 - a_volume_fraction : a_volume_fraction;
  double alpha = getAlpha(mm, VOFo);
  alpha = a_volume_fraction > 0.5 ? 1.0 - alpha : alpha;
  alpha += (factor - 0.5 * mm[0] - 0.5 * mm[1] - 0.5 * mm[2]);
  return alpha * norm + a_normal * a_rectangular_cuboid.calculateCentroid();
}

double getAlpha(const double* a_mm, const double a_VOFo) {
  // Copy and take absolute value of mm values
  double m[3] = {std::fabs(a_mm[0]), std::fabs(a_mm[1]), std::fabs(a_mm[2])};
  sort3Ascending(m);
  double m12 = {m[0] + m[1]};
  double V1 = {m[0] * m[0] / safelyEpsilon(6.0 * m[1] * m[2])};
  double V2 = {V1 + 0.5 * (m[1] - m[0]) / m[2]};
  double V3 = m12 <= m[2] ? 0.5 * m12 / m[2]
                          : (m[2] * m[2] * (3.0 * m12 - m[2]) +
                             m[0] * m[0] * (m[0] - 3.0 * m[2]) +
                             m[1] * m[1] * (m[1] - 3.0 * m[2])) /
                                safelyEpsilon(6.0 * m[0] * m[1] * m[2]);
  if (a_VOFo >= 0.0 && a_VOFo < V1) {
    return std::pow((6.0 * m[0] * m[1] * m[2] * a_VOFo), 1.0 / 3.0);
  } else if (a_VOFo >= V1 && a_VOFo < V2) {
    return 0.5 * (m[0] + sqrt(m[0] * m[0] + 8.0 * m[1] * m[2] * (a_VOFo - V1)));
  } else if (a_VOFo >= V2 && a_VOFo < V3) {
    double a0 = -(m[0] * m[0] * m[0] + m[1] * m[1] * m[1] -
                  6.0 * m[0] * m[1] * m[2] * a_VOFo);
    double a1 = 3.0 * (m[0] * m[0] + m[1] * m[1]);
    double a2 = -3.0 * m12;
    double p0 = -(a1 / 3.0 - a2 * a2 / 9.0);
    double q0 = (a1 * a2 - 3.0 * a0) / 6.0 - a2 * a2 * a2 / 27.0;
    double theta =
        std::acos(clipBetween(-1.0, q0 / std::sqrt(p0 * p0 * p0), 1.0)) / 3.0;
    return std::sqrt(p0) *
               (std::sqrt(3.0) * std::sin(theta) - std::cos(theta)) -
           a2 / 3.0;
  } else {
    if (m[2] >= m12) {
      return m[2] * a_VOFo + 0.5 * m12;
    } else {
      double a0 =
          -0.5 * (m[0] * m[0] * m[0] + m[1] * m[1] * m[1] + m[2] * m[2] * m[2] -
                  6.0 * m[0] * m[1] * m[2] * a_VOFo);
      double a1 = 1.5 * (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);
      double a2 = -1.5;
      double p0 = -(a1 / 3.0 - a2 * a2 / 9.0);
      double q0 = (a1 * a2 - 3.0 * a0) / 6.0 - a2 * a2 * a2 / 27.0;
      double theta =
          std::acos(clipBetween(-1.0, q0 / std::sqrt(p0 * p0 * p0), 1.0)) / 3.0;
      return std::sqrt(p0) *
                 (std::sqrt(3.0) * std::sin(theta) - std::cos(theta)) -
             a2 / 3.0;
    }
  }
}

}  // namespace IRL

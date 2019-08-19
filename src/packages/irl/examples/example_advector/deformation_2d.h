// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef EXAMPLES_EXAMPLE_ADVECTOR_DEFORMATION_2D_H_
#define EXAMPLES_EXAMPLE_ADVECTOR_DEFORMATION_2D_H_

#include "src/planar_reconstruction/planar_separator.h"
#include "examples/example_advector/basic_mesh.h"
#include "examples/example_advector/data.h"

struct Deformation2D {
  static BasicMesh setMesh(void);

  static void initialize(Data<double>* a_U, Data<double>* a_V,
                         Data<double>* a_W,
                         Data<IRL::PlanarSeparator>* a_separators);

  static void setVelocity(const double a_time, Data<double>* a_U,
                          Data<double>* a_V, Data<double>* a_W);
};

int main(void);

#endif  // EXAMPLES_EXAMPLE_ADVECTOR_DEFORMATION_2D_H_

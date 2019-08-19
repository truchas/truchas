// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "src/geometry/polyhedrons/tet.h"
#include "src/parameters/defined_types.h"

int main(void);

IRL::Tet generateRandomTet(const double minimum_volume);

double cutWithRecursiveSimplex(const IRL::UnsignedIndex_t number_of_iterations);

double cutWithSimplex(const IRL::UnsignedIndex_t number_of_iterations);

double cutWithHalfEdge(const IRL::UnsignedIndex_t number_of_iterations);

# Truchas Solid Mechanics Package

## Introduction

The solid mechanics package in this directory implements the thermoelastic small-displacement algorithm described by Baily and Cross [1]. Additional boundary conditions are introduced, allowing for sliding contact between different components.

This is a node-centered finite-volume algorithm, with a different set of access patterns compared to other Truchas packages. Nodes are surrounded by face-like integration points, which lie throughout the interior of neighboring cells. Each "integration point" is the center of a face created by joining edge center, adjacent face centers, and cell center.

This algorithm solves for displacement at nodes by integrating surface integrals of stress over all surrounding integration points. Stress is a local calculation given material Lame constants and strain. Strain is related to the displacement gradient, which must be calculated at all integration points.

The calculation of the displacement gradient at integration points within a given cell depends on the displacement at the nodes of that cell. This part of the algorithm is FEM-like, amounting to a multiplication of a matrix of shape-function gradients by the list of node-displacements. The shape function gradients are computed at initialization in `integration_geometry_type.F90` by inverting the Jacobian to a reference element where the shape functions are straighforwardly differentiated.

- [1] C. Bailey and M. Cross. A finite volume procedure to solve elastic solid mechanics problems in three dimensions on an unstructured mesh. International Journal for Numerical Methods in Engineering, 38:1757–1776, 1995.

## Tour and Logic Description

The most important files here are:

- `integration_geometry_type.F90` and `integration_cell_type.F90`: These modules handle the discretization and Jacobian calculation. One could think of the `integration_geometry` type
as analogous to an `unstr_mesh`, but holding connectivity data for the discretization.
- `sm_model_type.F90`: This is the top-level module for calculating the residual and forces.
- `sm_bc_manager_type.F90`: This module holds a type for managing boundary conditions, including
initializing, computing, and applying.
- `sm_bc_cXdY_type.F90`: These are a collection of files, each one applying a specific kind of boundary condition as described in the Physics & Algorithms manual, appendix K. The notation denotes the number of contact conditions and number of Dirichlet-displacement conditions at a given node. E.g., c1d2 handles nodes with 1 contact and 2 displacement conditions.

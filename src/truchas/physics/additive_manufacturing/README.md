## Toolheads

The code in this directory implements the `toolhead` object and its associated
models and facilities. The object is intended to describe a toolhead and
its attached energy and mass sources that are moving with respect to the
computational domain. The object couples a `toolpath` object that describes
the movement of the toolhead with reference configuration source functions
to yield moving sources. It serves then, in part, as a factory for allocating
standard abstract function objects (e.g., `vector_func`) that implement
the moving sources, and which can be used by the Truchas physics models.
It also manages updating the toolpath between path segments and the associated
updates to the sources during simulations.

Currently only a laser irradiation energy source is implemented; work on
a mass source for DED AM exists on a development branch.

A brief description of each of the modules follows. This should provide a rough
overview of how the code and its functionality is organized.

---

* `toolhead_type`
  This module defines the `toolhead` derived type.

* `laser_irrad_class`

  Laser irradiation source functions are vector-valued functions directed along
  the axis of the laser beam. Its flux density in planes orthogonal to the axis
  is a scalar. This module defines the abstract base class `laser_irrad` that
  establishes the interface for concrete implementations of this flux density
  function. The implementations are with respect to a reference configuration
  in which the z-axis is the axis of the beam and the origin is the focal
  point of the beam (if any). There are the following extensions of the class:

  * `gauss_laser_irrad_type`

    A flux density that is Gaussian in distance from the z-axis and invariant
    in the z coordinate.

  * `beam_laser_irrad_type`

    A flux density given by the Gaussian beam equation. See the reference
    manual for the detailed description.

* `laser_irrad_factory`

  This module provides the parameter list-driven subroutine `alloc_laser_irrad`
  that allocates a `laser_irrad` class object according the the parameter list.

* `toolhead_table_type`

  This module defines the `toolhead_table` derived type. As the name suggests,
  objects of this type hold a table of instantiated `toolhead` objects that are
  referenced by name. Since a `toolhead` is not necessarily specific to any one
  physics model, the table serves as the owner of the objects. Physics models
  will obtain a `toolhead` reference from the table in order to instantiate
  source functions (which they will own) but make little, if any, direct use
  of a `toolhead` object.

* `toolhead_namelist`

  This module provides the familiar input procedure that reads a collection
  of TOOLHEAD namelists from an input file. It converts the input to parameter
  lists, instantiates the specified `toolhead` objects, and adds them to the
  `toolhead_table` object held by the `toolhead_driver` module.

* `toolhead_driver`

  This module holds the Truchas `toolhead_table` object as a singleton, and
  provides several procedures which operate on that instance. Conceptually
  this module is part of the overarching Truchas multiphysics driver, and
  functions as an adapter between its present non-OO disign and the OO design
  of the toolhead component.

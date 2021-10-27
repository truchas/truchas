## Toolpaths

The code in this directory implements the `toolpath` object and associated
facilities. The object is designed to describe a path through space like that
taken by a machine tool in the course of a manufacturing process. The path
is specified using a simple command language that is compatible with common
CNC machine languages (e.g., G-code). The original motivating use case was the
motion of the laser/powder head in a directed energy deposition AM machine,
but `toolpath` objects have found use in non-AM applications as well.

A brief description of each of the modules follows. This should provide a rough
overview of how the code and its functionality is organized.

---

* `toolpath_type`

  This module defines the `toolpath` derived type. Internally the path is
  represented by a sequence of simple path segments. Those are objects of
  class `xyz_motion`.

* `xyz_motion_class.F90`

  This defines the abstract class `xyz_motion` that establishes the interface
  the path segments will implement, such as position given time, and the end
  points of the path segment. Its concrete extensions are:

  * `dwell_xyz_motion_type`

    This module defines the derived type `dwell_xyz_motion` that implements a
    constant path (no motion) over a time interval.`

  * `linear_xyz_motion_type`

    This module defines the derived type `linear_xyz_motion` that implements
    linear displacement from a point with a given speed and possible initial
    acceleration and final deceleration.

* `toolpath_factory`

  This module provides the parameter list driven subroutine `alloc_toolpath`
  that instantiates a `toolpath` object. It also implements the parsing of the
  JSON-format toolpath command file or string, which constitutes the major
  portion of the instantiation process. Any changes to the toolpath command
  language will be implemented here.

* `toolpath_factory_type`

  This module defines the derived type `toolpath_factory` that is a high-level
  wrapper around the `alloc_toolpath` subroutine. An object of this type
  holds a parameter list that consists of a collection of named sublists,
  each of which is a parameter list for a toolpath. Client code can use
  a factory object to instantiate a particular toolpath by referencing its
  name. This is useful in situations where input and instantiation based on
  that input occur at different times during program execution, as is
  currently the case for Truchas.

* `toolpath_namelist`

  This module provides the familiar input procedure that reads a collection
  of TOOLPATH namelists from an input file, converting the input to parameter
  lists and initializing a `toolpath_factory` object for later use by clients.

* `toolpath_event_type`

  This module defines the derived type `toolpath_event`, which extends the
  abstract base class `event_action` defined in the module `sim_event_queue`.
  Users of a `toolpath` object are responsible for registering events of this
  type at the toolpath segment endpoint times, which has the effect of
  advancing the toolpath to the next path segment.

* `toolpath_driver`

  This module holds the Truchas `toolpath_factory` object as a singleton, and
  provides several procedures which access that instance. Conceptually this
  module is part of the overarching Truchas multiphysics driver, and functions
  as an adapter between its present non-OO design and the OO design of the
  toolpath component.

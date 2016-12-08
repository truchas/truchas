The code in this directory implements a TOOLPATH object that is designed to
describe a path through space like that taken by a machine tool in the course
of a manufacturing process.  The path is specified using a simple command
lannguage that is adapted to common CNC machine languages.  The motivating
use case is the motion of the laser/powder head in an AM machine.  A brief
description of each of the source files follows.  This should provide a rough
overview of how the code and its functionality is organized.

Neil Carlson <nnc@lanl.gov>
December 2016

================================================================================

toolpath_type.F90
  This defines the fundamental TOOLPATH type.  Internally the path is
  represented as a sequence of simple path segments.  Those are objects
  of class XYZ_MOTION.  This also provides some basic tools that are
  used to instantiate a TOOLPATH object.

xyz_motion_class.F90
  This defines the abstract type XYZ_MOTION that establishes the interface
  the path segments will implement, such as position given time, and the
  end points of the path segment.

dwell_xyz_motion_type.F90
  This defines the type DWELL_XYZ_MOTION of class XYZ_MOTION that implements
  a constant path (no motion) over a time interval.

linear_xyz_motion_type.F90
  This defines the type LINEAR_XYZ_MOTION of class XYZ_MOTION that implements
  linear displacement from a point with speed and possible initial acceleration
  and final deceleration.

toolpath_factory.F90
  This provides a parameter list driven subroutine that instantantiates a
  TOOLPATH object.

toolpath_namelist.F90
  This provides standard Truchas input procedure which reads a collection
  of TOOLPATH namelists and instantiates the specified TOOLPATH objects
  using the factory method.

toolpath_table.F90
  This module holds the named toolpath objects instantiated by the input
  procedure.  Client code can then query the table to get a reference to
  one of the toolpaths.

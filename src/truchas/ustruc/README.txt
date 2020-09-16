MICROSTRUCTURE ANALYSIS KERNEL

Neil Carlson <nnc@lanl.gov>
August 2014

This code is the initial implementation of a framework for the prediction
of microstructure during solidification.  It is aimed at experimentally-
informed analytic point models that use local information like the thermal
gradient, solidification front velocity, and time in the mushy zone to
predict the type of microstructure growth (planar, dendritic, etc.) and
its characteristics (primary and secondary arm spacing, e.g.) that form
at a location.  This work is funded by LDRD 20140639ER "Solute and
Microstructure Prediction during Processing".

Code organization.  The modeling framework is implemented using an object-
oriented design.  The primary class is USTRUC_MODEL which encapsulates the
modeling kernel.  An object of that class is held by the USTRUC_DRIVER
module, which provides adapter-like procedures for driving the model from
from the top-most Truchas levels (which are the antithesis of OO design.)
All interaction with a USTRUC_MODEL object is via its type bound procedures;
it accesses no external data, consistent with OO principles.  That class is
responsible for generating derived state data requiring mesh-based
computations, before handing off the actual point-wise modeling to an object
of abstract class USTRUC_COMP.  For maximum flexibility and extensibility,
the goal is to decompose the modeling into individual analysis components
that could be combined dynamically at run time.  This design was implemented
using the decorator programming pattern.  In this pattern, the concrete
implementation USTRUC_CORE of USTRUC_COMP provides the core functionality.
Other optional analysis components are realized as concrete implementations
of the abstract class USTRUC_PLUGIN which is itself an extension of the
USTRUC_COMP class.  The USTRUC_COMP object that carries out the point-wise
microstructure modeling is a core object wrapped by zero or more plugin
class objects, like a set of nested russian dolls.

TODO

This is a work in progess, and some remains to be done.  A number of the
issues that need to be addressed are described in code comments; look for
the the tags "FIXME" or "TODO".  Other larger issues are listed below
roughly in order of importance.

Redistribution of work.  The microstructure modeling is limited to a fixed
subset of mesh cells where the subject material may be present.  Currently
each process works on the cells in its subdomain that belong to this subset.
For real applications, this will lead to a very unbalanced work load, and a
redistribution of that work should be considered at some point.

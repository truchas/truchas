# Materials

This directory contains the source code for the software component that
manages materials and their properties and attributes for the Truchas solvers.

Here is an overview of the component with links to pages describing the
individual pieces in detail.

* The `material` class and its companion `phase` derived type form the
  foundation of the "materials" code component. These types store the
  properties (e.g., thermal conductivity function) and attributes (e.g.,
  "fluid") of both simple single-phase materials and multiphase materials.
  The companion `phase_change` class stores additional information about
  the transformation connecting two phases of a multiphase material. These
  are described in detail by the [materials.md](./doc/materials.md) and
  [phase-change.md](./doc/phase-change.md) pages.

* The `material_factory` module provides procedures for instantiating
  `material` class objects, using inputs supplied by a parameter list.
  The procedures and the format of the parameter list input are described
  in detail by the [material_factory.md](doc/material_factory.md) page.

* The `material_database` derived type is a simple container for a collection
  of named `material` objects. It serves as a potentially large in-memory
  database of materials from which Truchas will select specific materials
  to include in a simulation. This is described by the
  [material_database.md](./doc/material_database.md) page. The
  `material_factory` module provides a procedure for populating a material
  database.

* The `material_model` derived type is a container that stores the specific
  materials involved in a Truchas simulation. The materials and constituent
  phases are numbered and can be referenced in an array-like fashion that
  allows the use of corresponding volume fraction arrays describing mixed
  material/phase mesh cells. This derived type also introduces the notion of
  a virtual "void" material that is important to many Truchas applications.
  This derived type and its methods are described in detail by the
  [material_model.md](./doc/material_model.md) page.
  Truchas currently uses a single instance of this type, which provides
  the only access to the materials and phases available to the solvers.

* The computation of material property values is the one aspect of materials
  that is essential to the solvers. The `matl_prop` class provides the
  interface for the evaluation of a material property. It is akin to the
  `scalar_func` class used throughout Truchas, but is a higher-level structure
  associated with the `material` class which is capable of handling properties
  of multiphase materials and their temperature-dependent mixture of phases.
  It is described in detail by the [matl_prop.md](./doc/matl_prop.md) page.
  More generally, the `avg_matl_prop` and `avg_phase_prop` derived types
  provide the interface for the evaluation of volume fraction averaged
  properties on mixed material and phase cells. These are described in detail
  by the [average-property-functions.md](doc/average-property-functions.md)
  page. These latter two derived types are the ones used most often by the
  solvers.

The remaining code in this directory should be regarded as being part of the
top-level Truchas driver which has an old-fashion, non-OO, implementation that
eschews passing data through procedure arguments in favor of using module
variables accessed through `use` association.

* The [`material_model_driver`](./material_model_driver.F90) module hosts the
  single instances of the `material_model` derived type currently used. It
  also hosts a single instance of the `material_database` that is initialized
  through the Truchas input file. It provides a single procedure that invokes
  the initialization of both module variables.

* The [`material_namelist`](./material_namelist.F90) module provides the
  input procedure that reads the `MATERIAL`, `PHASE`, and `PHASE_CHANGE`
  namelist data, and translates it into a parameter list that will later be
  passed to a `material_factory` module procedure that populates the
  `material_database` instance.

* The `material_utilities` module provides an ad hoc collection of high-level
  convenience procedures, mostly used during Truchas initialization, that
  operate on the single `material_model` instance. Its procedures are
  documented by the [material_utilities.md](./doc/material_utilities.md) page.

# Truchas Material Volume Fraction Design

## Background

A Truchas simulation generally involves multiple materials (a dozen would not
be unheard of) and any given cell may contain multiple materials. This is
described by the volume fraction of each material (summing to 1). Further each
material has an associated set of properties (viscosity, thermal conductivity,
etc.) and those properties need to be evaluated on a cell using those volume
fractions and property values. The latter is the key bridge between the
distribution of materials (volume fractions) and the physics model parameters.
This design question is about how most efficiently store that volume fraction
data (which may evolve from one time step to the next) and perhaps how to
associate it with the specific materials.

## Important considerations

* Materials are either fluid and able to move between cells, or immobile and
  fixed within cells. The former imply volume fractions that evolve in time,
  and the latter volume fractions that are constant.

* Materials are conserved

* A material either consists of a single phase (so the same as the material
  itself) or multiple phases. The phase composition is determined by physics
  models; e.g., it may be a function of temperature. Phases are not conserved.

  * some phases may be fluid and others immobile.
  * each phase may have its own set of property values, however the material
    model is able to evaluate a net property for the material as a whole, so
    that is outside of the scope of volume fractions.
  * Note that mainline truchas stores phase volume fractions, but I think that
    is not the right way to handle things. Individual physics that require only
    specific phases (like fluid for flow) should derive that information.

* Initially each mesh region is assigned a single material, either using cell
  sets so that those cells hold exactly 1 material, or using geometric
  primitives like a disk, where a cell cut by the geometric boundary will
  contain two or more materials.

* VOID is a special virtual material representing "empty space", typically
  used by the flow model where it is considered "fluid".

* Some physics models may be applied to only a subset of the domain, for
  example flow, where it is known precisely which materials may be involved.

* In parallel, an individual mesh partition may only contain a subset of the
  simulation materials, allowing for different structures between partitions.

* I would generally expect the volume fraction data structure(s) to belong with
  the simulation object and be provided to the individual physics models/solvers.

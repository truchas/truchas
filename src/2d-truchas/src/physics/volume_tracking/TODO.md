# Remaining work

These are the remaining meaningful follow-up items for the 2D volume-tracking
implementation. Changes to the mainline 3D volume tracker are out of scope
unless they are required to support 2D.

- Handle nonconvergence of the plane locator through a proper status/error
  return instead of silently ignoring the root solver result.
- Validate that the material-priority list is a permutation of the material
  indices.
- Decide whether the remaining geometric tracker controls currently hard-coded
  in `t2d_geometric_volume_tracker%init`—the interface-location iteration
  limit and nested-dissection setting—should be exposed as input parameters.
- Either implement material-specific inflow handling for
  `t2d_simple_volume_tracker` or explicitly document that the simple algorithm
  does not support it.
- Implement nested-dissection reconstruction for cells containing three or
  more materials. Onion-skin reconstruction is the current implementation.
- Revisit the conservation-repair failure path when no available material can
  accept a reassigned flux.

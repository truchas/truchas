meshes<M>-<N>.jou are Cubit input files to generate both the HT and EM meshes.

* M=1 is the original COMSOL geometry with a point contact between the
  sphere and plate. This has a cusp plate/sphere surface intersection that
  results in badly shaped tetrahedra -- something that should generally
  be avoided.

* M=2 is a modified geometry that removes a very thin spherical cap so
  that the sphere contacts the plate in a disk. The plate/sphere surfaces
  intersect in about a 15 degree angle, and is meshed without badly shaped
  tetrahedra.

* N is the number of points per wavelength. The linear mesh resolution is
  set in each region accordingly; e.g., the edge lengths in the sphere are
  about 8 times smaller than in the free-space region.

geom<M>.jou create the geometries. These are input by the meshes*.jou files
and are not used directly.

Note: All meshes use the same high-resolution mesh in the small waveguide
section with the intent of producing the same energy input into the oven
from the port-feed boundary condition.

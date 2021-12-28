**IMPORTANT NOTE**
----
The `grid_mapping_module.F90` source file used in Truchas differs from the
original version contained in the `grid_mapping_lanl.tar.gz` distribution
found in this directory. The latter has a defect in the `get_adjhex_relation`
subroutine that sometimes results in bogus neighbor information (leading
to later failures) for degenerate hex cells. Thus the original source code
should only be used with tet and/or **non**-degenerate hex meshes.

See https://gitlab.com/truchas/truchas/-/issues/57 for more details.

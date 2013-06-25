!!
!! EXODUS
!!
!! Neil N. Carlson <nnc@lanl.gov> 22 Sep 2004
!! Last revised 31 Oct 2005
!!
!! This module provides procedures for reading and writing Exodus II mesh
!! files [1], and derived data types that encapsulate the mesh data.  The
!! package uses the Fortran 90 interface to the netCDF library, defined in
!! the module "netcdf", to read and write the data.  That module is part
!! of the netCDF distribution.
!!
!! EXODUS is the only module that an application code needs to use.  See
!! the accompanying documentation for detailed instructions on its use.
!! The following is a brief list of the facilities provided:
!!
!! o Derived types EXODUS_MESH, ELEM_BLK, NODE_SET, SIDE_SET.
!! o Subroutine DESTROY for deallocating storage associated with these types.
!! o Function DEFINED for testing the definition state of these types.
!! o Defined relational operators == and /= for these types.
!! o Subroutine READ_EXODUS_MESH for reading a mesh from an Exodus II file.
!! o Subroutine WRITE_EXODUS_MESH for writing a mesh to an Exodus II file.
!! o Function EXO_ERR_STR for retrieving error strings when things go bad.
!! o Subroutine PRINT_EXODUS_MESH for dumping the data in a EXODUS_MESH variable.
!!
!! [1] L.A.Schoof and V.R.Yarberry, "Exodus II: A Finite Element Data Model",
!!     Sandia report SAND92-2137.  This can be obtained at
!!     http://endo.sandia.gov/SEACAS/Documentation/exodusII.pdf
!!
!! REVISION HISTORY
!! 10/31/2005 -- added exodus_mesh_utilities module
!!

module exodus

  use exodus_mesh_type
  use exodus_mesh_reader
  use exodus_mesh_writer
  use exodus_errors, only: exo_err_str
  use exodus_truchas_hack
  use exodus_mesh_utilities
  public

  !! These are from EXODUS_MESH_TYPE and aren't relevant to application code.
  private :: MAX_STR_LENGTH, MAX_LINE_LENGTH

end module exodus

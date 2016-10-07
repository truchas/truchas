// HDF5 Identifier handling 
//%import "hdf5.h"

// Can not get the above import to pick up the typedef for herr_t
// Without this typedef SWIG doesn't 'understand' what herr_t is
// Will need to change this if HDF5 changes the typedef
// herr_t defined in H5IPublic.h
%typedef int herr_t;


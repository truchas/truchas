// HDF5 Identifier handling 
//%import "hdf5.h"

// Can not get the above import to pick up the typedef for hid_t
// Without this typedef SWIG doesn't 'understand' what hid_t is
// Will need to change this if HDF5 changes the typedef
// hid_t defined in H5IPublic.h
%typedef int hid_t;


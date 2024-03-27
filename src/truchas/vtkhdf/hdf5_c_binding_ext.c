#include "hdf5.h"
hid_t H5P_DATASET_CREATE_value() { return H5P_DATASET_CREATE; }
hid_t H5T_NATIVE_INTEGER_value() { return H5Tcopy(H5T_NATIVE_INT); }
hid_t H5T_NATIVE_DOUBLE_value() { return H5Tcopy(H5T_NATIVE_DOUBLE); }
hid_t H5T_STD_U8LE_value() { return H5Tcopy(H5T_STD_U8LE); }
hid_t H5T_NATIVE_CHARACTER_value() {
  hid_t type_id;
  if ((type_id = H5Tcopy(H5T_FORTRAN_S1)) < 0) return type_id;
  if (H5Tset_size(type_id, 1) < 0) return type_id;
  if (H5Tset_strpad(type_id, H5T_STR_SPACEPAD) < 0) return type_id;
  return type_id;
}

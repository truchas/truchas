/********************************************************************/
/*
   Create a dataset with a string in it.
*/
/********************************************************************/

#include "hdf5.h"
#define FILE "ex.h5"

int main() {

   hid_t       file_id, dataset_id, dataspace_id;  /* identifiers */
   herr_t      status;
   char        buf[]={"test left call mesh " };
   hid_t       atype;
   size_t      size;
   hsize_t     dims[2] = {2, 2};

   file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   printf ("H5Fcreate returns: %i\n", file_id);
   dataspace_id = H5Screate_simple (2, dims, NULL);
   printf ("H5Screate_simple returns: %i\n", dataspace_id);

   atype = H5Tcopy (H5T_C_S1);
   printf ("H5Tcopy returns: %i\n", atype);
   size = 5;
   status = H5Tset_size (atype, size);
   printf ("H5Tset_size returns: %i\n", status);

   dataset_id = H5Dcreate(file_id, "strdata", atype, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
   printf ("H5Dcreate returns: %i\n", dataset_id);

   status = H5Dwrite (dataset_id, atype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
   printf ("H5Dwrite returns: %i\n", status);

   status = H5Dclose(dataset_id);
   printf ("H5Dclose returns: %i\n", status);
   status = H5Sclose(dataspace_id);
   printf ("H5Sclose returns: %i\n", status);
   status = H5Fclose(file_id);
   printf ("H5Fclose returns: %i\n", status);

   return 0;
}


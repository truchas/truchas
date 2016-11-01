/* UPGRADE VIEW FACTOR FILE (upgradevff.c)
 *
 * This little utility adds the "icount" array variable to a view factor file.
 * This array gives the number of nonzero elements in each row of the view
 * factor matrix.  It replaces the original "ia" array that is one of the
 * arrays in a standard CSR representation, which gives the starting index in
 * the view factor array of each row.  The advantage of the "icount" array is
 * that it can use 32-bit integers, whereas "ia" really needs to be 64-bit in
 * order to move to larger problems. The conversion between "ia" and "icount"
 * is trivial.
 *
 * Neil N. Carlson <nnc@lanl.gov>
 * October 2016
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include "netcdf.h"

void handle_status(int status, const char *string) {
  if (status != NC_NOERR) {
    printf("Error: %s: %s\n", string, nc_strerror(status));
    exit(1);
  }
}

int main(int argc, char **argv)
{
  int status, ncid, varid, dimid;
  size_t len, len2;

  /* check for the expected number of arguments */
  if (argc != 2) {
    printf("Usage: %s VIEW_FACTOR_FILE\n", argv[0]);
    exit(1);
  }

  /* check the file exists and has rw permissions */
  if (access(argv[1], F_OK|R_OK|W_OK)) {
    printf("Error: file %s does not exist or does not have rw permissions\n", argv[1]);
    exit(1);
  }

  /* open the netcdf file in rw mode */
  status = nc_open(argv[1], NC_WRITE, &ncid);
  handle_status(status, "error opening the netCDF file");

  /* check if the file is already upgraded */
  status = nc_inq_varid(ncid, "icount", &varid);
  if (status == NC_NOERR) {
    printf("view factor file is already upgraded\n");
    exit(0);
  }

  /* get the id and length of the ia variable */
  status = nc_inq_varid(ncid, "ia", &varid);
  handle_status(status, "file is missing view factor data");
  status = nc_inq_vardimid(ncid, varid, &dimid);
  status = nc_inq_dimlen(ncid, dimid, &len);
  handle_status(status, "error getting view factor data size");

  /* read the ia array */
  int* ia = (int*) malloc(len*sizeof(int));
  status = nc_get_var_int(ncid, varid, ia);
  handle_status(status, "error reading variable ia");

  /* initialize the icount array */
  int* icount = (int*) malloc((len-1)*sizeof(int));
  for (int j=0; j < len-1; j++) {
    icount[j] = ia[j+1] - ia[j];
  }

  /* get the num_faces dimension */
  status = nc_inq_dimid(ncid, "num_faces", &dimid);
  handle_status(status, "file is missing num_faces dimension");
  status = nc_inq_dimlen(ncid, dimid, &len2);
  handle_status(status, "error reading num_faces dimension");
  if (len2 != len-1) {
    printf("length of ia is inconsistent with num_faces\n");
    exit(1);
  }

  /* define and write the icount variable */
  status = nc_redef(ncid);
  handle_status(status, "error putting file into definition mode");
  status = nc_def_var(ncid, "icount", NC_INT, 1, &dimid, &varid);
  handle_status(status, "error defining icount variable");
  status = nc_enddef(ncid);
  handle_status(status, "error putting file into data mode");
  status = nc_put_var(ncid, varid, icount);
  handle_status(status, "error writing icount variable");

  /* close file and report success */
  status = nc_close(ncid);
  handle_status(status, "error closing file");
  printf("upgraded view factor file\n");
  exit(0);
}

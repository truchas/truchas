/****************************************************************************
 *   								                                        *
 *    vlslab.c  - Create dataset using variable length string datatype and  *
 *                writing by slab                                           *
 *   								                                        *
 ****************************************************************************/

#include <hdf5.h>

#define FILE   "vlslab.h5"

/* 1-D dataset with starting size of 4 */
#define SPACE1_NAME  "Space1"
#define SPACE1_RANK	1
#define SPACE1_DIM1	4
#define SLABS       3


int main (void)
{
    const char *wdata [SPACE1_DIM1]= {
        "Four score and seven years ago our forefathers brought forth on this continent a new nation,",
        "conceived in liberty and dedicated to the proposition that all men are created equal.",
        "Now we are engaged in a great civil war,",
        "testing whether that nation or any nation so conceived and so dedicated can long endure."
        };   
    const char *wdata1 [SPACE1_DIM1]= {
        "Five score and seven years ago our forefathers brought forth on this continent a new nation,",
        "conceived in liberty and dedicated to the proposition that all men are created equal.",
        "Now we are engaged in a great civil war,",
        "testing whether that nation or any nation so conceived and so dedicated can long endure."
        };   
    const char *wdata2 [SPACE1_DIM1]= {
        "Nine score and seven years ago our forefathers brought forth on this continent a new nation,",
        "conceived in liberty and dedicated to the proposition that all men are created equal.",
        "Now we are engaged in a great civil war,",
        "testing whether that nation or any nation so conceived and so dedicated can long endure."
        };   

    char        *rdata [SPACE1_DIM1*SLABS]; /* Information read in */
    hid_t		fid1;		      
    hid_t		dataset;	    
    hid_t		sid1, space, memspace;
    hid_t		tid1, cparms, xfer_pid;  
    hsize_t		dims1[] = {SPACE1_DIM1};
    hsize_t		maxdims1[1] = {SPACE1_DIM1*SLABS};
    hsize_t		chkdims1[1] = {2};
    int         i;          /* counting variable */
    herr_t		ret;		/* Generic return value  */
    hssize_t    offset[1] = {0};
    hssize_t    offset1[1] = {SPACE1_DIM1};
    hsize_t     count[1] = {SPACE1_DIM1};


    /* Create file */
    fid1 = H5Fcreate (FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    printf ("H5Fcreate (fid1): %i\n", fid1);

    /* Create dataspace for datasets */
    sid1 = H5Screate_simple (SPACE1_RANK, maxdims1, NULL);
    printf ("H5Screate_simple (sid1): %i\n", sid1);

    /* Create a datatype to refer to */
    tid1 = H5Tcopy (H5T_C_S1);
    printf ("H5Tcopy (tid1): %i\n", tid1);

    ret = H5Tset_size (tid1,H5T_VARIABLE);
    printf ("H5Tset_size: %i\n", ret);

    cparms = H5Pcreate (H5P_DATASET_CREATE);
    printf ("H5Pcreate returns: %i\n", cparms);

    ret = H5Pset_chunk ( cparms, 1, chkdims1);
    printf ("H5Pset_chunk returns: %i\n", ret);

    /* Create a dataset */
    dataset = H5Dcreate (fid1, SPACE1_NAME, tid1, sid1, H5P_DEFAULT,cparms,H5P_DEFAULT);

    memspace = H5Screate_simple (SPACE1_RANK, dims1, NULL);
    printf ("H5Screate_simple: %i\n", memspace);

    for (i=0; i<SLABS; i++)
    {
       space = H5Screate_simple (SPACE1_RANK, maxdims1, NULL);
       printf ("H5Screate_simple (sid1): %i\n", space);

       offset[0] = offset1[0]*i;
       ret = H5Sselect_hyperslab (space, H5S_SELECT_SET, offset, NULL,
                count, NULL);

       printf ("H5Sselect_hyperslab: %i\n", ret);

       if (i==0)
         ret = H5Dwrite (dataset, tid1, memspace, space, H5P_DEFAULT, wdata);
       else if (i==1)
         ret = H5Dwrite (dataset, tid1, memspace, space, H5P_DEFAULT, wdata1);
       else if (i==2)
         ret = H5Dwrite (dataset, tid1, memspace, space, H5P_DEFAULT, wdata2);

       printf ("H5Dwrite: %i\n", ret);

       ret = H5Sclose (space);
       printf ("H5Sclose: %i\n", ret);

    }

    xfer_pid = H5Pcreate (H5P_DATASET_XFER);
    printf ("H5Pcreate: %i\n", xfer_pid);

    ret= H5Dread (dataset, tid1, H5S_ALL, H5S_ALL, xfer_pid, rdata);
    printf ("H5Dread returns: %i\n", ret);

    printf ("DATA:\n ");

    for (i=0; i< SPACE1_DIM1*SLABS; i++)
      printf ("%s\n ",rdata[i]);
    printf ("\n");

    space = H5Dget_space (dataset);
    printf ("H5Dget_space: %i\n", space);
    
    ret = H5Dvlen_reclaim (tid1, space, xfer_pid, rdata);
    printf ("H5Dvlen_reclaim: %i\n", ret);

    /* Close Everything */
    ret = H5Dclose (dataset);
    printf ("Close everything: %i ", ret); 
    ret = H5Tclose (tid1);
    printf ("%i ", ret);
    ret = H5Sclose (sid1);
    printf ("%i ", ret);
    ret = H5Sclose (space);
    printf ("%i ", ret);
    ret = H5Sclose (memspace);
    printf ("%i ", ret);
    ret = H5Pclose (cparms);
    printf ("%i ", ret);
    ret = H5Pclose (xfer_pid);
    printf ("%i ", ret);
    ret = H5Fclose (fid1);
    printf ("%i\n", ret);

    return 0;

}


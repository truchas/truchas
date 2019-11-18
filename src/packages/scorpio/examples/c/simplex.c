/*******************************************************************************
* Copyright Notice
*  + 2010-2012 North Carolina State University
*  + 2010-2012 Pacific Northwest National Laboratory
* 
* This file is part of SCORPIO.
* 
* SCORPIO is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* any later version.
* 
* SCORPIO is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
* 
* You should have received a copy of the GNU Lesser General Public License
* along with SCORPIO.  If not, see <http://www.gnu.org/licenses/>.
* 
*******************************************************************************/
/** 
 * @file simplex.c
 * @brief Simple test driver for parallel I/O
 * @author Sarat Sreepathi
 * @version 1.0
 * @date 2010-12-02
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "scorpio.h"

/* Number of dimensions of dataset*/
#define NDIMS 2

int main(int argc, char **argv)
{
	int errcode;
	int fhandle;
	char filename[80] = "sample.h5";
	int i,j;
	int localcount;
	int ndims;
	int globaldims[NDIMS];
	int localdims[NDIMS];
	double **matrix;
	double *vector;
	int groupid;
	int verifyError;
	int nprocs;

	verifyError = 0;

	iogroup_conf_t myIOconfig;
	iogroup_t myIOgroup;

	errcode = MPI_Init(&argc, &argv);
	if (errcode != MPI_SUCCESS) { fprintf(stderr, "Error initializing MPI\n"); exit(2); }

	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	myIOconfig.numIOgroups = 1;
	myIOconfig.commIncoming = MPI_COMM_WORLD;

	scorpio_IOgroup_init(&myIOconfig, &myIOgroup);

	/* Configure global dataset dimensions */
	ndims = NDIMS;

	// Sample non-uniform case, even ranked nodes have 3 rows and odd ranked nodes have 2 rows each
	/* Even number of processes */
	if ( nprocs%2 == 0)
		globaldims[0] = nprocs/2 *  (3 + 2) ;
	/* Odd number of processes */
	else
		/* This also handles the case of 1 process */
		globaldims[0] = (nprocs+1)/2 *  (3 + 2) - 2 ;

	/* second dimension is arbitrary */
	globaldims[1] = 24;
	// ... and so on for the rest of the dimensions

	/* Distribute dataset uniformly using row-wise partitioning */
	/* Data from each node is aggregated at the IO nodes */
	/* Then written collectively by IO nodes to the file */
	/* For uniform row-wise partitioning, localdims[0] = globaldims[0]/myIOgroup.globalsize;  */

	// Sample non-uniform case, even ranked nodes have 3 rows and odd ranked nodes have 2 rows each
	if (myIOgroup.globalrank%2 == 0) 
		localdims[0] = 3;
	else
		localdims[0] = 2;

	for (i=1; i < ndims; i++)
		localdims[i] = globaldims[i];

	localcount = 1;
	for (i=0; i < ndims; i++)
		localcount *= localdims[i];

	vector = (double *) malloc( localdims[0] * localdims[1] * sizeof(double));
	matrix = (double **) malloc ( localdims[0] * sizeof(double *));

	for (i=0; i < localdims[0] ; i++)
		matrix[i] = (double *) malloc( localdims[1] * sizeof(double) );	

	/* Initialize data using arbitrary values */
	for (i=0; i < localdims[0] ; i++)
	{
		for (j=0; j < localdims[1] ; j++)
		{
			matrix[i][j] = (i+1)*(myIOgroup.globalrank+1) + j;
			vector[i*localdims[1]+j] = matrix[i][j];
		}
	}

	/* --------------------------------------------------------------------------- */
	/* Parallel I/O Calls  */
	/* --------------------------------------------------------------------------- */

	fhandle = scorpio_open_file(filename, &myIOgroup, SCORPIO_FILE_CREATE);
	if (fhandle == -1)
	{
		fprintf(stderr, "Opening file for writing failed.\n");
		exit(-1);
	}

	/* groupid = scorpio_create_dataset_group("mygroup", fhandle, &myIOgroup); */

	scorpio_write_dataset( vector, SCORPIO_DOUBLE, ndims, globaldims, localdims, fhandle, "testdata" ,  &myIOgroup, SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE);

	/* scorpio_close_dataset_group(groupid, fhandle, &myIOgroup); */
	scorpio_close_file( fhandle, &myIOgroup);

	/* Reset the vector before we read back into it for verification */
	memset(vector, 0, localcount * sizeof(double));

	/* Open file for reading now */
	fhandle = scorpio_open_file(filename, &myIOgroup, SCORPIO_FILE_READONLY);
	if (fhandle == -1)
	{
		fprintf(stderr, "Opening file for reading failed.\n");
		exit(-1);
	}

	/*
	int flag;
	flag = scorpio_dataset_exists("testdata", fhandle, &myIOgroup);
	printf("\n Dataset exists: %d \n", flag);
	flag = scorpio_dataset_exists("testdata2", fhandle, &myIOgroup);
	printf("\n Dataset exists: %d \n", flag);
	*/

	scorpio_get_dataset_ndims(&ndims, fhandle, "testdata", &myIOgroup);
	scorpio_get_dataset_dims(globaldims, fhandle, "testdata", &myIOgroup);

	int tmpsize;
	scorpio_get_dataset_size(&tmpsize, fhandle, "testdata", &myIOgroup);

	/*
	printf("\n Dataset size: %d", tmpsize);
	printf("\n Number of dimensions: %d", ndims);
	printf("\n Dims: %d %d \n", globaldims[0], globaldims[1]);
	fflush(stdout);
	*/

	/* When using groups, remember to update dataset name by prepending "/groupname/" for reading */

	scorpio_read_dataset(vector, SCORPIO_DOUBLE, ndims, globaldims, localdims, fhandle, "testdata" , &myIOgroup, SCORPIO_NONUNIFORM_CONTIGUOUS_READ);

	scorpio_close_file( fhandle, &myIOgroup);

	/* Verify */
	for (i=0; i < localdims[0] ; i++)
	{
		for (j=0; j < localdims[1] ; j++)
		{
			if ( vector[i*localdims[1]+j] != (i+1)*(myIOgroup.globalrank+1) + j )
			{
				fprintf(stderr, "Validation failed\n");	;
				verifyError = 1;
			}
		}
	}

	/* Print vector to stdout */
	/*
	if ( myIOgroup.localrank == 0 )
	{
		for (i=0; i < localdims[0] ; i++)
		{
			printf("Rank %d: Row %d :",myIOgroup.globalrank, i );
			for (j=0; j < localdims[1] ; j++)
			{
				printf("%.2f ", vector[i*localdims[1]+j]);
			}
			printf("\n");
		}
	}
	*/

	/* Cleanup */
	scorpio_IOgroup_cleanup(&myIOgroup);

	for (i=0; i < localdims[0] ; i++)
		free(matrix[i]);
	free(matrix);

	free(vector);

	MPI_Finalize();

	if (verifyError > 0)
	{
		fprintf(stderr, "Test failed.\n");
		return -1;
	}
	else
	{
		fprintf(stderr, "Test success.\n");
		return 0;
	}
}



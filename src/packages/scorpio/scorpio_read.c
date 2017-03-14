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
#include "scorpio.h"

int scorpio_read_dataset1( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int fhandle, char *dset_name,  iogroup_t *myIOgroup);
int scorpio_read_dataset2( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int fhandle, char *dset_name,  iogroup_t *myIOgroup);
int scorpio_read_entire_dataset( void *vector, datatype_t mytype, int ndims, int *globaldims, int fhandle, char *dset_name,  iogroup_t *myIOgroup);

int scorpio_read_dataset( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int fhandle, char *dset_name,  iogroup_t *myIOgroup, iopattern_t mypattern)
{
	if ( mypattern == SCORPIO_UNIFORM_CONTIGUOUS_READ)
		return scorpio_read_dataset1(vector, mytype, ndims, globaldims, localdims, fhandle, dset_name, myIOgroup);
	else if ( mypattern == SCORPIO_NONUNIFORM_CONTIGUOUS_READ)
		return scorpio_read_dataset2(vector, mytype, ndims, globaldims, localdims, fhandle, dset_name, myIOgroup);
	else if ( mypattern == SCORPIO_EVERYONE_ENTIRE_DATASET_READ)
		return scorpio_read_entire_dataset(vector, mytype, ndims, globaldims, fhandle, dset_name, myIOgroup);

	/* You shouldn't get here */
	return -1;
}

int scorpio_get_dataset_size( int *mydataset_size, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{
	hssize_t dataset_size;

	if ( myIOgroup->localrank != 0)
	{
		MPI_Bcast(&dataset_size, 1, MPI_INT,  0, myIOgroup->localcomm);
	}
	/* only do this for IO nodes */
	else
	{	 
		iofile_t *currfile;
		currfile = myIOgroup->file[fhandle];

		herr_t ret;         	/* Generic return value */

		// Open dataset collectively 
		currfile->dataset = H5Dopen2(currfile->fid, dset_name, H5P_DEFAULT); 
		assert(currfile->dataset != FAILURE);

		/* create a file dataspace independently */
		currfile->file_dataspace = H5Dget_space (currfile->dataset);
		assert(currfile->file_dataspace != FAILURE);

		/* H5Sget_simple_extent_dims(currfile->file_dataspace, dims, maxdims) */
		dataset_size = H5Sget_simple_extent_npoints(currfile->file_dataspace);

		/* release all temporary handles. */
		H5Sclose(currfile->file_dataspace);

		/* All reads completed.  Close dataset collectively */
		ret=H5Dclose(currfile->dataset);
		assert(ret != FAILURE);

		MPI_Bcast(&dataset_size, 1, MPI_INT,  0, myIOgroup->localcomm);

	}

	*mydataset_size = dataset_size;

	return 0;

}

int scorpio_get_dataset_dims( int *dataset_dims, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{
	int dataset_ndims;
	hsize_t *tmpdims;
	int i;

	if ( myIOgroup->localrank != 0)
	{
		MPI_Bcast(&dataset_ndims, 1, MPI_INT,  0, myIOgroup->localcomm);
		MPI_Bcast(dataset_dims, dataset_ndims, MPI_INT,  0, myIOgroup->localcomm);
		/* MPI_Bcast(dataset_dims, dataset_ndims * sizeof(hsize_t), MPI_BYTE,  0, myIOgroup->localcomm); */
	}
	/* only do this for IO nodes */
	else
	{	 
		iofile_t *currfile;
		currfile = myIOgroup->file[fhandle];

		herr_t ret;         	/* Generic return value */

		// Open dataset collectively 
		currfile->dataset = H5Dopen2(currfile->fid, dset_name, H5P_DEFAULT); 
		assert(currfile->dataset != FAILURE);

		/* create a file dataspace independently */
		currfile->file_dataspace = H5Dget_space (currfile->dataset);
		assert(currfile->file_dataspace != FAILURE);

		dataset_ndims = H5Sget_simple_extent_ndims(currfile->file_dataspace);
		tmpdims = (hsize_t *) calloc(dataset_ndims, sizeof(hsize_t));

		dataset_ndims = H5Sget_simple_extent_dims(currfile->file_dataspace, tmpdims, NULL);

		for ( i=0; i<dataset_ndims; i++)
		{
			/* cast back to integer for now */
			dataset_dims[i] = (int) tmpdims[i];
			PRINT_MSG((SCORPIO_INFO, "dataset dims[i]: %d ", i, dataset_dims[i]));
		}

		/* printf("\nTEST tmp dims[0]: %d dims[1]: %d\n", tmpdims[0],tmpdims[1]); */
		/* printf("\nTEST dataset dims[0]: %d dims[1]: %d", dataset_dims[0], dataset_dims[1]); */
		free(tmpdims);

		/* release all temporary handles. */
		H5Sclose(currfile->file_dataspace);

		/* All reads completed.  Close dataset collectively */
		ret=H5Dclose(currfile->dataset);
		assert(ret != FAILURE);

		MPI_Bcast(&dataset_ndims, 1, MPI_INT,  0, myIOgroup->localcomm);
		MPI_Bcast(dataset_dims, dataset_ndims, MPI_INT,  0, myIOgroup->localcomm);
		/* MPI_Bcast(dataset_dims, dataset_ndims * sizeof(hsize_t), MPI_BYTE,  0, myIOgroup->localcomm); */
	}

	return 0;
}

int scorpio_get_dataset_ndims( int *pndims, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{
	int dataset_ndims; 

	if ( myIOgroup->localrank != 0)
	{
		MPI_Bcast(&dataset_ndims, 1, MPI_INT,  0, myIOgroup->localcomm);
	}
	/* only do this for IO nodes */
	else
	{	 
		iofile_t *currfile;
		currfile = myIOgroup->file[fhandle];

		herr_t ret;         	/* Generic return value */

		// Open dataset collectively 
		currfile->dataset = H5Dopen2(currfile->fid, dset_name, H5P_DEFAULT); 
		assert(currfile->dataset != FAILURE);

		/* create a file dataspace independently */
		currfile->file_dataspace = H5Dget_space (currfile->dataset);
		assert(currfile->file_dataspace != FAILURE);

		dataset_ndims = H5Sget_simple_extent_ndims(currfile->file_dataspace);

		/* release all temporary handles. */
		H5Sclose(currfile->file_dataspace);

		/* All reads completed.  Close dataset collectively */
		ret=H5Dclose(currfile->dataset);
		assert(ret != FAILURE);

		MPI_Bcast(&dataset_ndims, 1, MPI_INT,  0, myIOgroup->localcomm);
	}

	*pndims = dataset_ndims;

	return 0;

}

// Everyone reads the entire dataset -- generalized datatype
int scorpio_read_entire_dataset( void *vector, datatype_t mytype, int ndims, int *globaldims, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{

	void *buffer = NULL;
	int localcount;
	int i;

	initialize_datatype(mytype, myIOgroup);

	localcount = 1;
	// Essentially localcount = globaldims[0] * globaldims[1] ... (as we are reading the entire dataset)
	for(i=0; i<ndims; i++)
		localcount *= globaldims[i];

	if ( myIOgroup->localrank != 0)
	{
		MPI_Scatter(buffer, localcount, myIOgroup->mpi_type, vector, localcount, myIOgroup->mpi_type, 0, myIOgroup->localcomm );
	}
	/* only do this for IO nodes */
	else
	{	 
		int i;
		iofile_t *currfile;
		currfile = myIOgroup->file[fhandle];

		buffer = calloc(myIOgroup->localsize*localcount, myIOgroup->datatype_size);
		if ( buffer == NULL)
		{
			perror("Read_entire_dataset-calloc failure");
			PRINT_MSG(( SCORPIO_ERROR, "Failed trying to allocate %d bytes for the buffer.", myIOgroup->localsize*localcount*myIOgroup->datatype_size));
			MPI_Abort(myIOgroup->globalcomm, 1111);
		}

		/* dataspace dim sizes */
		hsize_t *dimsf; 
		hsize_t *start,*count, *stride;			/* for hyperslab setting */
		herr_t ret;         	/* Generic return value */

		MPI_Info info = MPI_INFO_NULL;

		dimsf = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		count = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		start = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		stride = (hsize_t *) calloc(ndims, sizeof(hsize_t) );

		/* Set up dimensions of the slab this process accesses.  */
		/* Dataset: each process takes a block of rows. */

		// Set the dimension values from incoming argument
		for( i = 0; i < ndims; i++)
			dimsf[i] = globaldims[i];

		// Contiguous in all dimensions, so stride is 1 for all dimensions
		for( i = 0; i < ndims; i++)
			stride[i] = 1;

		// Everyone reads the entire dataset
		// For all dimensions, count[i] is equal to globaldims which is copied to dimsf[i]
		for( i = 0; i < ndims; i++)
			count[i] = dimsf[i];

		// Everyone reads the entire dataset
		// So for all dimensions, offset is zero as everyone reads all elements in that dimension.
		for( i = 0; i < ndims; i++)
			start[i] = 0;

		// Open dataset collectively 
		currfile->dataset = H5Dopen2(currfile->fid, dset_name, H5P_DEFAULT); 
		assert(currfile->dataset != FAILURE);

		PRINT_MSG(( SCORPIO_INFO, "Read - IO Node: %d - start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu", myIOgroup->iogroupRank, (unsigned long)start[0], (unsigned long)start[1], (unsigned long)count[0], (unsigned long)count[1], (unsigned long)(count[0]*count[1]) ));

		/* create a file dataspace independently */
		currfile->file_dataspace = H5Dget_space (currfile->dataset);
		assert(currfile->file_dataspace != FAILURE);
		ret=H5Sselect_hyperslab(currfile->file_dataspace, H5S_SELECT_SET, start, stride, count, NULL);
		assert(ret != FAILURE);

		/* create a memory dataspace independently */
		currfile->mem_dataspace = H5Screate_simple (ndims, count, NULL);
		assert (currfile->mem_dataspace != FAILURE);

		/* set up the collective transfer properties list */
		currfile->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
		assert(currfile->xfer_plist != FAILURE);
#if !defined(SERIAL_HDF5)
		ret=H5Pset_dxpl_mpio(currfile->xfer_plist, H5FD_MPIO_COLLECTIVE);
		assert(ret != FAILURE);
#endif		

		/* write data collectively */
		ret = H5Dread(currfile->dataset, myIOgroup->hdf_type, currfile->mem_dataspace, currfile->file_dataspace, currfile->xfer_plist, buffer);
		assert(ret != FAILURE);

		/* release all temporary handles. */
		H5Sclose(currfile->file_dataspace);
		H5Sclose(currfile->mem_dataspace);
		H5Pclose(currfile->xfer_plist);

		/* All reads completed.  Close dataset collectively */
		ret=H5Dclose(currfile->dataset);
		assert(ret != FAILURE);

		MPI_Scatter(buffer, localcount, myIOgroup->mpi_type, vector, localcount, myIOgroup->mpi_type, 0, myIOgroup->localcomm );

		free(buffer);
		free(dimsf);
		free(start);
		free(count);
		free(stride);
	}

	/* MPI_Barrier(myIOgroup->globalcomm); */
	return 0;
}

// Everyone reads the same portion of a dataset 
int scorpio_read_same_sub_dataset( void *vector, datatype_t mytype, int ndims, int *globaldims, int *mystart, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{

	void *buffer = NULL;
	int localcount;
	int i;

	initialize_datatype(mytype, myIOgroup);

	localcount = 1;
	// Essentially localcount = globaldims[0] * globaldims[1] ... (as we are reading the entire dataset)
	for(i=0; i<ndims; i++)
		localcount *= globaldims[i];

	if ( myIOgroup->localrank != 0)
	{
		MPI_Scatter(buffer, localcount, myIOgroup->mpi_type, vector, localcount, myIOgroup->mpi_type, 0, myIOgroup->localcomm );
	}
	/* only do this for IO nodes */
	else
	{	 
		int i;
		iofile_t *currfile;
		currfile = myIOgroup->file[fhandle];

		buffer = calloc(myIOgroup->localsize*localcount, myIOgroup->datatype_size);
		if ( buffer == NULL)
		{
			perror("Read_entire_dataset-calloc failure");
			PRINT_MSG(( SCORPIO_ERROR, "Failed trying to allocate %d bytes for the buffer.", myIOgroup->localsize*localcount*myIOgroup->datatype_size));
			MPI_Abort(myIOgroup->globalcomm, 1111);
		}

		/* dataspace dim sizes */
		hsize_t *dimsf; 
		hsize_t *start,*count, *stride;			/* for hyperslab setting */
		herr_t ret;         	/* Generic return value */

		MPI_Info info = MPI_INFO_NULL;

		dimsf = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		count = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		start = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		stride = (hsize_t *) calloc(ndims, sizeof(hsize_t) );

		/* Set up dimensions of the slab this process accesses.  */
		/* Dataset: each process takes a block of rows. */

		// Set the dimension values from incoming argument
		for( i = 0; i < ndims; i++)
			dimsf[i] = globaldims[i];

		// Contiguous in all dimensions, so stride is 1 for all dimensions
		for( i = 0; i < ndims; i++)
			stride[i] = 1;

		// Everyone reads the entire dataset
		// For all dimensions, count[i] is equal to globaldims which is copied to dimsf[i]
		for( i = 0; i < ndims; i++)
			count[i] = dimsf[i];

		// So for all dimensions, offset is read from the incoming argument
		for( i = 0; i < ndims; i++)
			start[i] = mystart[i];

		// Open dataset collectively 
		currfile->dataset = H5Dopen2(currfile->fid, dset_name, H5P_DEFAULT); 
		assert(currfile->dataset != FAILURE);

		PRINT_MSG(( SCORPIO_INFO, "Read - IO Node: %d - start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu", myIOgroup->iogroupRank, (unsigned long)start[0], (unsigned long)start[1], (unsigned long)count[0], (unsigned long)count[1], (unsigned long)(count[0]*count[1]) ));

		/* create a file dataspace independently */
		currfile->file_dataspace = H5Dget_space (currfile->dataset);
		assert(currfile->file_dataspace != FAILURE);
		ret=H5Sselect_hyperslab(currfile->file_dataspace, H5S_SELECT_SET, start, stride, count, NULL);
		assert(ret != FAILURE);

		/* create a memory dataspace independently */
		currfile->mem_dataspace = H5Screate_simple (ndims, count, NULL);
		assert (currfile->mem_dataspace != FAILURE);

		/* set up the collective transfer properties list */
		currfile->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
		assert(currfile->xfer_plist != FAILURE);
#if !defined(SERIAL_HDF5)
		ret=H5Pset_dxpl_mpio(currfile->xfer_plist, H5FD_MPIO_COLLECTIVE);
		assert(ret != FAILURE);
#endif		

		/* write data collectively */
		ret = H5Dread(currfile->dataset, myIOgroup->hdf_type, currfile->mem_dataspace, currfile->file_dataspace, currfile->xfer_plist, buffer);
		assert(ret != FAILURE);

		/* release all temporary handles. */
		H5Sclose(currfile->file_dataspace);
		H5Sclose(currfile->mem_dataspace);
		H5Pclose(currfile->xfer_plist);

		/* All reads completed.  Close dataset collectively */
		ret=H5Dclose(currfile->dataset);
		assert(ret != FAILURE);

		MPI_Scatter(buffer, localcount, myIOgroup->mpi_type, vector, localcount, myIOgroup->mpi_type, 0, myIOgroup->localcomm );

		free(buffer);
		free(dimsf);
		free(start);
		free(count);
		free(stride);
	}

	/* MPI_Barrier(myIOgroup->globalcomm); */
	return 0;
}

// Everyone reads the entire dataset
int scorpio_read_dataset_all( void *vector, datatype_t mytype, int ndims, int *globaldims, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{

	void *buffer = NULL;
	int localcount;
	int i;

	initialize_datatype(mytype, myIOgroup);

	localcount = 1;
	// Essentially localcount = globaldims[0] * globaldims[1] ... (as we are reading the entire dataset)
	for(i=0; i<ndims; i++)
		localcount *= globaldims[i];

	if ( myIOgroup->localrank != 0)
	{
		MPI_Scatter(buffer, localcount, myIOgroup->mpi_type, vector, localcount, myIOgroup->mpi_type, 0, myIOgroup->localcomm );
	}
	/* only do this for IO nodes */
	else
	{	 
		int i;
		iofile_t *currfile;
		currfile = myIOgroup->file[fhandle];

		buffer =  calloc(myIOgroup->localsize*localcount, myIOgroup->datatype_size);

		/* dataspace dim sizes */
		hsize_t *dimsf; 
		hsize_t *start,*count, *stride;			/* for hyperslab setting */
		herr_t ret;         	/* Generic return value */

		MPI_Info info = MPI_INFO_NULL;

		dimsf = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		count = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		start = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		stride = (hsize_t *) calloc(ndims, sizeof(hsize_t) );

		/* Set up dimensions of the slab this process accesses.  */
		/* Dataset: each process takes a block of rows. */

		// Set the dimension values from incoming argument
		for( i = 0; i < ndims; i++)
			dimsf[i] = globaldims[i];

		// Contiguous in all dimensions, so stride is 1 for all dimensions
		for( i = 0; i < ndims; i++)
			stride[i] = 1;

		// Everyone reads the entire dataset
		// For all dimensions, count[i] is equal to globaldims which is copied to dimsf[i]
		for( i = 0; i < ndims; i++)
			count[i] = dimsf[i];

		// Everyone reads the entire dataset
		// So for all dimensions, offset is zero as everyone reads all elements in that dimension.
		for( i = 0; i < ndims; i++)
			start[i] = 0;

		// Open dataset collectively 
		currfile->dataset = H5Dopen2(currfile->fid, dset_name, H5P_DEFAULT); 
		assert(currfile->dataset != FAILURE);

		PRINT_MSG(( SCORPIO_INFO, "Read - IO Node: %d - start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu", myIOgroup->iogroupRank, (unsigned long)start[0], (unsigned long)start[1], (unsigned long)count[0], (unsigned long)count[1], (unsigned long)(count[0]*count[1]) ));

		/* create a file dataspace independently */
		currfile->file_dataspace = H5Dget_space (currfile->dataset);
		assert(currfile->file_dataspace != FAILURE);
		ret=H5Sselect_hyperslab(currfile->file_dataspace, H5S_SELECT_SET, start, stride, count, NULL);
		assert(ret != FAILURE);

		/* create a memory dataspace independently */
		currfile->mem_dataspace = H5Screate_simple (ndims, count, NULL);
		assert (currfile->mem_dataspace != FAILURE);

		/* set up the collective transfer properties list */
		currfile->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
		assert(currfile->xfer_plist != FAILURE);
#if !defined(SERIAL_HDF5)
		ret=H5Pset_dxpl_mpio(currfile->xfer_plist, H5FD_MPIO_COLLECTIVE);
		assert(ret != FAILURE);
#endif		

		/* write data collectively */
		ret = H5Dread(currfile->dataset, myIOgroup->hdf_type, currfile->mem_dataspace, currfile->file_dataspace, currfile->xfer_plist, buffer);
		assert(ret != FAILURE);

		/* release all temporary handles. */
		H5Sclose(currfile->file_dataspace);
		H5Sclose(currfile->mem_dataspace);
		H5Pclose(currfile->xfer_plist);

		/* All reads completed.  Close dataset collectively */
		ret=H5Dclose(currfile->dataset);
		assert(ret != FAILURE);

		MPI_Scatter(buffer, localcount, myIOgroup->mpi_type, vector, localcount, myIOgroup->mpi_type, 0, myIOgroup->localcomm );

		free(buffer);
		free(dimsf);
		free(start);
		free(count);
		free(stride);
	}

	/* MPI_Barrier(myIOgroup->globalcomm); */
	return 0;
}

int scorpio_read_dataset1( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{

	void *buffer = NULL;
	int localcount;
	int i;

	initialize_datatype(mytype, myIOgroup);

	localcount = 1;
	// Essentially localcount = localdims[0] * localdims[1] ...
	for(i=0; i<ndims; i++)
		localcount *= localdims[i];

	if ( myIOgroup->localrank != 0)
	{
		MPI_Scatter(buffer, localcount, myIOgroup->mpi_type, vector, localcount, myIOgroup->mpi_type, 0, myIOgroup->localcomm );
	}
	/* only do this for IO nodes */
	else
	{	 
		int i;
		iofile_t *currfile;
		currfile = myIOgroup->file[fhandle];

		buffer =  calloc(myIOgroup->localsize*localcount, myIOgroup->datatype_size);

		/* dataspace dim sizes */
		hsize_t *dimsf; 
		hsize_t *start,*count, *stride;			/* for hyperslab setting */
		herr_t ret;         	/* Generic return value */

		MPI_Info info = MPI_INFO_NULL;

		dimsf = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		count = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		start = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		stride = (hsize_t *) calloc(ndims, sizeof(hsize_t) );

		/* Set up dimensions of the slab this process accesses.  */
		/* Dataset: each process takes a block of rows. */

		// Set the dimension values from incoming argument
		for( i = 0; i < ndims; i++)
			dimsf[i] = globaldims[i];

		// Contiguous in all dimensions, so stride is 1 for all dimensions
		for( i = 0; i < ndims; i++)
			stride[i] = 1;

		// Count in first dimension is equal to the number of rows at local ionode
		/* count[0] = dimsf[0]/myIOgroup->numIOgroups; */
		/* Note: Sarat changed Sept 12, 2011: Support for different group size in last IO group */
		count[0] = dimsf[0]/myIOgroup->globalsize * myIOgroup->localsize;

		// For other dimensions, count[i] is equal to globaldims which is copied to dimsf[i]
		for( i = 1; i < ndims; i++)
			count[i] = dimsf[i];

		// Calculate offset at IO node by summing up prior IO node counts
		/* start[0] = myIOgroup->iogroupRank * count[0];  */
		/* Note: Sarat changed Sept 12, 2011: Support for different group size in last IO group */
		/* last group has to be handled specially */
		if (myIOgroup->iogroupRank == (myIOgroup->numIOgroups-1) )
			/* start[0] = myIOgroup->iogroupRank * ( dimsf[0]/(myIOgroup->numIOgroups-1) ); */
			start[0] = myIOgroup->iogroupRank * ( dimsf[0]/myIOgroup->globalsize * myIOgroup->preferredGroupSize );
		else
			start[0] = myIOgroup->iogroupRank * count[0]; 


		// For other dimensions, offset is zero as everyone writes all elements in that dimension.
		for( i = 1; i < ndims; i++)
			start[i] = 0;

		// Open dataset collectively 
		currfile->dataset = H5Dopen2(currfile->fid, dset_name, H5P_DEFAULT); 
		assert(currfile->dataset != FAILURE);

		PRINT_MSG(( SCORPIO_INFO, "Read - IO Node: %d - start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu", myIOgroup->iogroupRank, (unsigned long)start[0], (unsigned long)start[1], (unsigned long)count[0], (unsigned long)count[1], (unsigned long)(count[0]*count[1]) ));

		/* create a file dataspace independently */
		currfile->file_dataspace = H5Dget_space (currfile->dataset);
		assert(currfile->file_dataspace != FAILURE);
		ret=H5Sselect_hyperslab(currfile->file_dataspace, H5S_SELECT_SET, start, stride, count, NULL);
		assert(ret != FAILURE);

		/* create a memory dataspace independently */
		currfile->mem_dataspace = H5Screate_simple (ndims, count, NULL);
		assert (currfile->mem_dataspace != FAILURE);

		/* set up the collective transfer properties list */
		currfile->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
		assert(currfile->xfer_plist != FAILURE);
#if !defined(SERIAL_HDF5)
		ret=H5Pset_dxpl_mpio(currfile->xfer_plist, H5FD_MPIO_COLLECTIVE);
		assert(ret != FAILURE);
#endif		

		/* write data collectively */
		ret = H5Dread(currfile->dataset, myIOgroup->hdf_type, currfile->mem_dataspace, currfile->file_dataspace, currfile->xfer_plist, buffer);
		assert(ret != FAILURE);

		/* release all temporary handles. */
		H5Sclose(currfile->file_dataspace);
		H5Sclose(currfile->mem_dataspace);
		H5Pclose(currfile->xfer_plist);

		/* All reads completed.  Close dataset collectively */
		ret=H5Dclose(currfile->dataset);
		assert(ret != FAILURE);

		MPI_Scatter(buffer, localcount, myIOgroup->mpi_type, vector, localcount, myIOgroup->mpi_type, 0, myIOgroup->localcomm );

		free(buffer);
		free(dimsf);
		free(start);
		free(count);
		free(stride);
	}

	/* MPI_Barrier(myIOgroup->globalcomm); */
	return 0;
}

// SCORPIO_NONUNIFORM_CONTIGUOUS_READ
int scorpio_read_dataset2( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{

	void *buffer = NULL;
	int *slavecount = NULL;
	int localcount;
	int i,j;
	int *displs = NULL;
	int *slave_nrows = NULL;

	initialize_datatype(mytype, myIOgroup);

	localcount = 1;
	// Essentially localcount = localdims[0] * localdims[1] ...
	for(i=0; i<ndims; i++)
		localcount *= localdims[i];


	/* only do this for IO nodes */
	if ( myIOgroup->localrank != 0)
	{
		MPI_Gather(&localcount, 1, MPI_INT, slavecount, 1, MPI_INT, 0, myIOgroup->localcomm );
		MPI_Gather(&localdims[0], 1, MPI_INT, slave_nrows, 1, MPI_INT, 0, myIOgroup->localcomm );
		PRINT_MSG(( SCORPIO_VERBOSE, "Read Before Scatterv"));
		MPI_Scatterv(buffer, slavecount, displs, myIOgroup->mpi_type, vector, localcount, myIOgroup->mpi_type, 0, myIOgroup->localcomm );
		PRINT_MSG(( SCORPIO_VERBOSE, "Read After Scatterv"));
	}
	else
	{	 
		iofile_t *currfile;

		int local_ionode_count, local_ionode_nrows;
		int *global_ionode_count, *global_ionode_nrows;

		// slavecount contains number of elements at each slave
		slavecount = (int *) calloc(myIOgroup->localsize, sizeof(int) );
		// slave_nrows contains number of rows at each slave
		slave_nrows = (int *) calloc(myIOgroup->localsize, sizeof(int) );
		displs = (int *) calloc(myIOgroup->localsize, sizeof(int) );
		global_ionode_count = (int *) calloc(myIOgroup->numIOgroups, sizeof(int));
		global_ionode_nrows = (int *) calloc(myIOgroup->numIOgroups, sizeof(int));

		// Gather all localcount from slave nodes in the local group.
		MPI_Gather(&localcount, 1, MPI_INT, slavecount, 1, MPI_INT, 0, myIOgroup->localcomm );
		MPI_Gather(&localdims[0], 1, MPI_INT, slave_nrows, 1, MPI_INT, 0, myIOgroup->localcomm );

		// Calculate displacements for MPI_Scatterv
		// Each displs is equal to sum of localcounts of prior members.
		// NOTE: displs[i] is set to zero using calloc
		for(i=0; i < myIOgroup->localsize; i++)
			for ( j=0; j<i; j++)
				displs[i] += slavecount[j];

		// Calculate total number of elements at the IO node level (i.e., sum of all slave node counts)
		local_ionode_count = 0;
		local_ionode_nrows = 0;
		for(i=0; i< myIOgroup->localsize; i++)
		{
			local_ionode_count += slavecount[i];
			local_ionode_nrows += slave_nrows[i];
		}

		// Allocate buffer to hold all elements at IO node.
		buffer =  calloc(local_ionode_count, myIOgroup->datatype_size );

		PRINT_MSG(( SCORPIO_INFO, "Read Rank:%d  localcount: %d slavecount[0]: %d displs[0] : %d ", myIOgroup->globalrank, localcount, slavecount[0], displs[0])); 

		PRINT_MSG(( SCORPIO_VERBOSE, "IOnode Before Allgather"));
		// Gather all ionode_counts from all IO nodes 
		MPI_Allgather(&local_ionode_count, 1, MPI_INT, global_ionode_count, 1, MPI_INT, myIOgroup->iocomm );
		PRINT_MSG(( SCORPIO_VERBOSE, "IOnode After Allgather"));

		PRINT_MSG(( SCORPIO_VERBOSE, "IOnode Before Allgather2"));
		// Gather all ionode_nrows from all IO nodes 
		MPI_Allgather(&local_ionode_nrows, 1, MPI_INT, global_ionode_nrows, 1, MPI_INT, myIOgroup->iocomm );
		PRINT_MSG(( SCORPIO_VERBOSE, "IOnode After Allgather2"));

		currfile = myIOgroup->file[fhandle];

		/* dataspace dim sizes */
		hsize_t *dimsf; 
		hsize_t *start,*count, *stride;         /* for hyperslab setting */
		herr_t ret;             /* Generic return value */

		MPI_Info info = MPI_INFO_NULL;

		dimsf = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		count = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		start = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		stride = (hsize_t *) calloc(ndims, sizeof(hsize_t) );

		/* Set up dimensions of the slab this process accesses.  */
		/* Dataset: each process takes a block of rows. */

		// Set the dimension values from incoming argument
		for( i = 0; i < ndims; i++)
			dimsf[i] = globaldims[i];

		// Contiguous in all dimensions, so stride is 1 for all dimensions
		for( i = 0; i < ndims; i++)
			stride[i] = 1;

		// Calculate offset at IO node by summing up prior IO node counts
		for(i=0; i < myIOgroup->iogroupRank; i++)
			start[0] += global_ionode_nrows[i];

		// For other dimensions, offset is zero as everyone writes all elements in that dimension.
		for( i = 1; i < ndims; i++)
			start[i] = 0;

		// Count in first dimension is equal to the number of rows at local ionode
		count[0] = local_ionode_nrows;

		// For other dimensions, count[i] is equal to globaldims which is copied to dimsf[i]
		for( i = 1; i < ndims; i++)
			count[i] = dimsf[i];

		PRINT_MSG((SCORPIO_INFO, "Rank:%d  iogroup_Rank: %d local_ionode_count : %d start[0] %d start[1] : %d", myIOgroup->globalrank, myIOgroup->iogroupRank, local_ionode_nrows, start[0], start[1])); 

		PRINT_MSG(( SCORPIO_INFO, "Write - IO Node: %d - start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu", myIOgroup->iogroupRank, \
					(unsigned long)start[0], (unsigned long)start[1], \
					(unsigned long)count[0], (unsigned long)count[1], \
					(unsigned long)(count[0]*count[1]) ));

		// Open dataset collectively 
		currfile->dataset = H5Dopen2(currfile->fid, dset_name, H5P_DEFAULT); 
		assert(currfile->dataset != FAILURE);

		/* create a file dataspace independently */
		currfile->file_dataspace = H5Dget_space (currfile->dataset);
		assert(currfile->file_dataspace != FAILURE);
		ret=H5Sselect_hyperslab(currfile->file_dataspace, H5S_SELECT_SET, start, stride, count, NULL);
		assert(ret != FAILURE);

		/* create a memory dataspace independently */
		currfile->mem_dataspace = H5Screate_simple (ndims, count, NULL);
		assert (currfile->mem_dataspace != FAILURE);

		/* set up the collective transfer properties list */
		currfile->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
		assert(currfile->xfer_plist != FAILURE);
#if !defined(SERIAL_HDF5)
		ret=H5Pset_dxpl_mpio(currfile->xfer_plist, H5FD_MPIO_COLLECTIVE);
		assert(ret != FAILURE);
#endif		

		/* write data collectively */
		ret = H5Dread(currfile->dataset, myIOgroup->hdf_type, currfile->mem_dataspace, currfile->file_dataspace, currfile->xfer_plist, buffer);
		assert(ret != FAILURE);

		/* release all temporary handles. */
		H5Sclose(currfile->file_dataspace);
		H5Sclose(currfile->mem_dataspace);
		H5Pclose(currfile->xfer_plist);

		/* All reads completed.  Close dataset collectively */
		ret=H5Dclose(currfile->dataset);
		assert(ret != FAILURE);

		/* Scatter data to all slave nodes.  */
		PRINT_MSG(( SCORPIO_VERBOSE, "IOnode Read Before Scatterv"));
		MPI_Scatterv(buffer, slavecount, displs, myIOgroup->mpi_type, vector, localcount, myIOgroup->mpi_type, 0, myIOgroup->localcomm );
		PRINT_MSG(( SCORPIO_VERBOSE, "IOnode Read After Scatterv"));

		free(buffer);
		free(dimsf);
		free(start);
		free(count);
		free(stride);

		free(slavecount);
		free(slave_nrows);
		free(displs);
		free(global_ionode_count);
		free(global_ionode_nrows);

	}

	/* MPI_Barrier(myIOgroup->globalcomm); */
	return 0;
}

int scorpio_read_str_array(char ***pstrarray, int *pnstrings, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{
	int nstrings;
	int i,j;
	void* buffer;
	char* pbuffer;
	char **strarray;
	int maxstr_len;

	if ( myIOgroup->globalrank == 0)
	{
		/* only do this on root process - rank 0 */
		iofile_t *currfile;
		herr_t ret;

		hsize_t ndims;
		/* We are creating a string array where each row contains one string */
		/* int dims[1]; */
		ndims=1;

		currfile = myIOgroup->file[fhandle];

		/* scorpio_get_dataset_dims(dims, fhandle, dset_name, myIOgroup); */
		currfile->dataset = H5Dopen2(currfile->fid, dset_name, H5P_DEFAULT); 

		currfile->file_dataspace = H5Dget_space (currfile->dataset);

		/* tmpdims = (hsize_t *) calloc(ndims, sizeof(hsize_t)); */
		/* H5Sget_simple_extent_dims(currfile->file_dataspace, dims, NULL); */
		/* nstrings = dims[0]; */
		nstrings = H5Sget_simple_extent_npoints(currfile->file_dataspace);

		/* myIOgroup->hdf_type = H5Tcopy(H5T_C_S1); */
		myIOgroup->hdf_type = H5Dget_type(currfile->dataset);
		maxstr_len = H5Tget_size(myIOgroup->hdf_type);
		/* maxstr_len = 38; */

		/* printf("\n Maximum string length: %d ", maxstr_len); */
		PRINT_MSG((SCORPIO_VERBOSE, "Maximum string length: %d ", maxstr_len));

		currfile = myIOgroup->file[fhandle];
		buffer = (char *) malloc( nstrings * maxstr_len * sizeof(char));

		/* printf("Reading dataset : %s\n", dset_name); fflush(stdout); */
		ret = H5Dread(currfile->dataset, myIOgroup->hdf_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
		/* printf("I read: %s", buffer); fflush(stdout); */

		/* * All writes completed.  Close dataset */
		ret=H5Dclose(currfile->dataset);
		H5Tclose(myIOgroup->hdf_type);

	}

	MPI_Bcast(&nstrings, 1, MPI_INT,  0, myIOgroup->globalcomm);
	MPI_Bcast(&maxstr_len, 1, MPI_INT,  0, myIOgroup->globalcomm);

	/* Root already has this allocated */
	if ( myIOgroup->globalrank != 0)
		buffer = (char *) malloc( nstrings * maxstr_len * sizeof(char));

	MPI_Bcast(buffer, (nstrings * maxstr_len * sizeof(char)), MPI_CHAR,  0, myIOgroup->globalcomm);

	pbuffer = buffer;
	strarray = (char **) malloc( nstrings * sizeof(char*));
	/* send back to calling function  */
	*pstrarray = strarray;
	*pnstrings = nstrings;

	for (i=0; i < nstrings; i++)
	{
		strarray[i] = (char *) malloc( maxstr_len * sizeof(char));
		strncpy(strarray[i], pbuffer, maxstr_len);
		pbuffer += maxstr_len;

		/* We don't need to add terminating \0 because that is included in the string during write */
		/* If this is a fresh dataset that was not written by our library,  */
		   /* then the following line acts as a safeguard */
		strarray[i][maxstr_len-1] = '\0';

		/* printf("\tRead: %s$", strarray[i]); fflush(stdout); */
		/* trim space on the right side */
		/*
		   j = maxstr_len-1;
		   while(strarray[i][j] == ' ')
		   {
			   j--;
			   continue;
		   }
		   strarray[i][j+1] = '\0';
		   */

		/* would not handle space inside strings */
		/* sscanf(pbuffer, "%s", strarray[i]); */

	}

	/* Processing buffer is no longer needed */
	free(buffer);

 
	return 0;
}

/* primary object name can be a dataset or group or committed datatype name */
int scorpio_read_simple_attr(char *attr_name, void **pattr_data, datatype_t mytype, int fhandle, char *primary_obj_name,  iogroup_t *myIOgroup)
{
	int ndims;
	int *dummydims;
	/* printf("\n Point 1 \n"); fflush(stdout); */
	/* This setups the following call to read a scalar attribute */
	return scorpio_read_attr(attr_name, pattr_data, mytype, &ndims, &dummydims, fhandle, primary_obj_name, myIOgroup);
}

int scorpio_read_attr(char *attr_name, void **pattr_data, datatype_t mytype, int *pndims, int **padims, int fhandle, char *primary_obj_name,  iogroup_t *myIOgroup)
{
	int i;
	void *attr_data;
	int ndims;
	int *adims = NULL;
	size_t size, strsize;

	initialize_datatype(mytype, myIOgroup);
	if ( myIOgroup->globalrank == 0)
	{
		/* only do this on root process - rank 0 */
		iofile_t *currfile;

		herr_t status;
		hid_t aspace_id, attr_id, loc_id;
		hid_t space_type;

		hsize_t *dims = NULL;

		currfile = myIOgroup->file[fhandle];

		loc_id = H5Oopen(currfile->fid, primary_obj_name, H5P_DEFAULT);

		attr_id = H5Aopen(loc_id, attr_name, H5P_DEFAULT);
		myIOgroup->hdf_type = H5Aget_type(attr_id);
		aspace_id = H5Aget_space(attr_id);

		space_type = H5Sget_simple_extent_type(aspace_id);

		if ( space_type == H5S_SCALAR)
		{
			ndims = 0;
			adims = NULL;
		}
		else if ( space_type == H5S_SIMPLE) 
		{
			ndims = H5Sget_simple_extent_ndims(aspace_id);
			dims = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
			/* adims is outgoing to calling function  */
			adims = (int *) calloc(ndims, sizeof(int) );

			H5Sget_simple_extent_dims(aspace_id, dims, NULL);

			/* Copy to outgoing dims */
			for (i=0;i<ndims;i++) 
				adims[i] = dims[i];

			free(dims);
			/* The other items allocated will be freed up by the caller once they are done */
		}

		size = H5Aget_storage_size(attr_id);
		attr_data = malloc(size);

		/* printf("\n Point 6 size %d \n",size); fflush(stdout); */

		if (H5T_STRING == H5Tget_class(myIOgroup->hdf_type))
			strsize = H5Tget_size(myIOgroup->hdf_type);

		status = H5Aread(attr_id, myIOgroup->hdf_type, attr_data);

		status = H5Sclose (aspace_id);
		status = H5Aclose(attr_id);
		status = H5Oclose(loc_id);

		MPI_Bcast(&ndims, 1, MPI_INT, 0, myIOgroup->globalcomm);
		MPI_Bcast(adims, ndims, MPI_INT, 0, myIOgroup->globalcomm);
		MPI_Bcast(&size, sizeof(size), MPI_BYTE, 0, myIOgroup->globalcomm);
		MPI_Bcast(attr_data, size, MPI_BYTE, 0, myIOgroup->globalcomm);
	}
	else
	{
		MPI_Bcast(&ndims, 1, MPI_INT, 0, myIOgroup->globalcomm);
		adims = (int *) calloc(ndims, sizeof(int) );

		MPI_Bcast(adims, ndims, MPI_INT, 0, myIOgroup->globalcomm);
		MPI_Bcast(&size, sizeof(size), MPI_BYTE, 0, myIOgroup->globalcomm);
		attr_data = malloc(size);
		MPI_Bcast(attr_data, size, MPI_BYTE, 0, myIOgroup->globalcomm);
	}

	*pattr_data = attr_data;
	*pndims = ndims;
	*padims = adims;

	return 0;
}

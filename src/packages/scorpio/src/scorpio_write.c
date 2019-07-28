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

int scorpio_write_dataset1( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int fhandle, char *dset_name,  iogroup_t *myIOgroup);
int scorpio_write_dataset2( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int fhandle, char *dset_name,  iogroup_t *myIOgroup);

int scorpio_write_dataset( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int fhandle, char *dset_name,  iogroup_t *myIOgroup, iopattern_t mypattern)
{
	if ( mypattern == SCORPIO_UNIFORM_CONTIGUOUS_WRITE)
		return scorpio_write_dataset1(vector, mytype, ndims, globaldims, localdims, fhandle, dset_name, myIOgroup);
	else if ( mypattern == SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE)
		return scorpio_write_dataset2(vector, mytype, ndims, globaldims, localdims, fhandle, dset_name, myIOgroup);

	/* You shouldn't get here */
	return -1;
}

int scorpio_close_dataset_group( int64_t groupid, int fhandle, iogroup_t *myIOgroup)
{
	if ( myIOgroup->localrank == 0)
	{
		iofile_t *currfile;
		currfile = myIOgroup->file[fhandle];

		currfile->groupid = -1;
		return H5Gclose((hid_t)groupid);
	}

	return 0;
}

int64_t scorpio_create_dataset_group( char *group_name, int fhandle, iogroup_t *myIOgroup)
{
	if ( myIOgroup->localrank == 0)
	{
		iofile_t *currfile;
		currfile = myIOgroup->file[fhandle];

		H5Eset_auto2(H5E_DEFAULT, NULL, stderr);

		/* Refer: http://www.hdfgroup.org/HDF5/doc/RM/RM_H5L.html#Link-Exists
		1.8.0 	Function introduced in this release.
		H5Lexists checks the existence of only the final element in a relative or absolute path; it does not check any other path elements. 
		The function will therefore fail when both of the following conditions exist:
			name is not local to the group specified by loc_id or, if loc_id is something other than a group identifier, name is not local to the root group.
			Any element of the relative path or absolute path in name, except the target link, does not exist. 
		*/

		/* Check if group already exists, if it does, open it */
		currfile->groupid = H5Gopen2(currfile->fid, group_name, H5P_DEFAULT);

		if (currfile->groupid == FAILURE)
		{
			/* set up the link properties list so that intermediate groups are created while creating a dataset */
			currfile->link_plist = H5Pcreate (H5P_LINK_CREATE);
			assert(currfile->link_plist != FAILURE);
			H5Pset_create_intermediate_group(currfile->link_plist, 1);

			/* create a dataset collectively */
			currfile->groupid = H5Gcreate2(currfile->fid, group_name, currfile->link_plist, H5P_DEFAULT, H5P_DEFAULT);
			/* release all temporary handles. */
			H5Pclose(currfile->link_plist);

		}
		assert(currfile->groupid != FAILURE);
		
		H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t) H5Eprint2, stderr);

		return (int64_t)currfile->groupid;

	}

	return 0;
}

int scorpio_create_link(char *target, int64_t link_loc_id, char *link_name, int fhandle, iogroup_t *myIOgroup)
{
	if ( myIOgroup->localrank == 0)
	{
		iofile_t *currfile;
		currfile = myIOgroup->file[fhandle];

		H5Eset_auto2(H5E_DEFAULT, NULL, stderr);

                herr_t ret = H5Lcreate_soft(target, (hid_t)link_loc_id, link_name, H5P_DEFAULT, H5P_DEFAULT);
                assert(ret != FAILURE);

		H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t) H5Eprint2, stderr);
	}

	return 0;
}

int scorpio_write_dataset1( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{
	/* Some compilers generate warning messages about variables like buffer being used before initializing */
	/* but this happens only on slave nodes where such variables are inconsequential */
	/* To avoid such warning messages, just initialize those pointers to NULL */
	void *buffer = NULL;
	int localcount;
	int i;

	PRINT_MSG(( SCORPIO_INFO, "Begin - Write uniform pattern" ));
	initialize_datatype(mytype, myIOgroup);

	localcount = 1;
	// Essentially localcount = localdims[0] * localdims[1] ...
	for(i=0; i<ndims; i++)
		localcount *= localdims[i];

	if ( myIOgroup->localrank != 0)
	{
		PRINT_MSG(( SCORPIO_VERBOSE, "Before Gather"));
		MPI_Gather(vector, localcount, myIOgroup->mpi_type, buffer, localcount, myIOgroup->mpi_type, 0, myIOgroup->localcomm );
		PRINT_MSG(( SCORPIO_VERBOSE, "After Gather"));

	}
	else
		/* only do this for IO nodes */
	{	 
		int i;
		iofile_t *currfile;

		buffer =  calloc(myIOgroup->localsize * localcount, myIOgroup->datatype_size );

		PRINT_MSG(( SCORPIO_INFO, "Rank:%d  localdims: %d %d, localsize: %d ", myIOgroup->globalrank, localdims[0], localdims[1], myIOgroup->localsize)); 

		/* Gather data from all slave nodes in the local group */
		PRINT_MSG(( SCORPIO_VERBOSE, "Before Gather"));
		MPI_Gather(vector, localcount, myIOgroup->mpi_type, buffer, localcount, myIOgroup->mpi_type, 0, myIOgroup->localcomm );
		PRINT_MSG(( SCORPIO_VERBOSE, "After Gather"));

		currfile = myIOgroup->file[fhandle];

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

		/* setup dimensionality object */
		currfile->sid = H5Screate_simple (ndims , dimsf, NULL);
		assert (currfile->sid != FAILURE);

		/* set up the link properties list so that intermediate groups are created while creating a dataset */
		currfile->link_plist = H5Pcreate (H5P_LINK_CREATE);
		assert(currfile->link_plist != FAILURE);
		H5Pset_create_intermediate_group(currfile->link_plist, 1);

		/* create a dataset collectively */
		currfile->dataset = H5Dcreate2(currfile->fid, dset_name, myIOgroup->hdf_type, currfile->sid, currfile->link_plist, H5P_DEFAULT, H5P_DEFAULT);
		assert(currfile->dataset != FAILURE);

		PRINT_MSG(( SCORPIO_INFO, "Write - IO Node: %d - start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu", myIOgroup->iogroupRank, \
					(unsigned long)start[0], (unsigned long)start[1], \
					(unsigned long)count[0], (unsigned long)count[1], \
					(unsigned long)(count[0]*count[1]) ));

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
		ret = H5Dwrite(currfile->dataset, myIOgroup->hdf_type, currfile->mem_dataspace, currfile->file_dataspace, currfile->xfer_plist, buffer);
		assert(ret != FAILURE);

		/* release all temporary handles. */
		H5Sclose(currfile->file_dataspace);
		H5Sclose(currfile->mem_dataspace);
		H5Pclose(currfile->xfer_plist);
		H5Pclose(currfile->link_plist);

		/* All writes completed.  Close datasets collectively */
		ret=H5Dclose(currfile->dataset);
		assert(ret != FAILURE);

		/* release all IDs created */
		H5Sclose(currfile->sid);

		free(buffer);

		free(start);
		free(stride);
		free(count);
		free(dimsf);
	}

	/* Optional: In case we wish to block all slaves till the write is complete. */
	/* MPI_Barrier(myIOgroup->globalcomm); */

	PRINT_MSG(( SCORPIO_INFO, "End - Write uniform pattern" ));

	return 0;
}

// SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE 
int scorpio_write_dataset2( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{
	void *buffer = NULL;
	int *slavecount = NULL;
	int localcount;
	int i,j;
	int *displs = NULL;
	int *slave_nrows = NULL;

	PRINT_MSG(( SCORPIO_INFO, "Begin - Write nonuniform pattern" ));
	initialize_datatype(mytype, myIOgroup);

	localcount = 1;
	// Essentially localcount = localdims[0] * localdims[1] ...
	for(i=0; i<ndims; i++)
		localcount *= localdims[i];

	if ( myIOgroup->localrank != 0)
	{
		// Note to self: remember to use & for scalar variables
		MPI_Gather(&localcount, 1, MPI_INT, slavecount, 1, MPI_INT, 0, myIOgroup->localcomm );
		MPI_Gather(&localdims[0], 1, MPI_INT, slave_nrows, 1, MPI_INT, 0, myIOgroup->localcomm );
		PRINT_MSG(( SCORPIO_VERBOSE, "Before Gatherv"));
		MPI_Gatherv(vector, localcount, myIOgroup->mpi_type, buffer, slavecount, displs, myIOgroup->mpi_type, 0, myIOgroup->localcomm );
		PRINT_MSG(( SCORPIO_VERBOSE, "After Gatherv"));
	}
	/* only do this for IO nodes */
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

		// Calculate displacements for MPI_Gatherv
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

		/* PRINT_MSG(( SCORPIO_INFO, " Rank:%d  localcount: %d slavecount[0]: %d displs[0] : %d ", myIOgroup->globalrank, localcount, slavecount[0], displs[0])); */
		/* Gather data from all slave nodes in the local group */
		PRINT_MSG(( SCORPIO_VERBOSE, "IOnode Before Gatherv"));
		MPI_Gatherv(vector, localcount, myIOgroup->mpi_type, buffer, slavecount, displs, myIOgroup->mpi_type, 0, myIOgroup->localcomm );
		PRINT_MSG(( SCORPIO_VERBOSE, "IOnode After Gatherv"));

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

		PRINT_MSG(( SCORPIO_INFO, " Rank:%d  iogroup_Rank: %d local_ionode_count : %d start[0] %d start[1] : %d", myIOgroup->globalrank, myIOgroup->iogroupRank, local_ionode_nrows, start[0], start[1])); 

		PRINT_MSG(( SCORPIO_INFO, "Write - IO Node: %d - start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu", myIOgroup->iogroupRank, \
					(unsigned long)start[0], (unsigned long)start[1], \
					(unsigned long)count[0], (unsigned long)count[1], \
					(unsigned long)(count[0]*count[1]) ));

		/* setup dimensionality object */
		currfile->sid = H5Screate_simple (ndims , dimsf, NULL);
		assert (currfile->sid != FAILURE);

		/* set up the link properties list so that intermediate groups are created while creating a dataset */
		currfile->link_plist = H5Pcreate (H5P_LINK_CREATE);
		assert(currfile->link_plist != FAILURE);
		H5Pset_create_intermediate_group(currfile->link_plist, 1);

		/* create a dataset collectively */
		if ( currfile->groupid != -1 )
			currfile->dataset = H5Dcreate2(currfile->groupid, dset_name, myIOgroup->hdf_type, currfile->sid, currfile->link_plist, H5P_DEFAULT, H5P_DEFAULT);
		else
			currfile->dataset = H5Dcreate2(currfile->fid, dset_name, myIOgroup->hdf_type, currfile->sid, currfile->link_plist, H5P_DEFAULT, H5P_DEFAULT);
		/* currfile->dataset = H5Dcreate2(currfile->fid, dset_name, myIOgroup->hdf_type, currfile->sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); */
		assert(currfile->dataset != FAILURE);

		/* create a file dataspace independently */
		currfile->file_dataspace = H5Dget_space (currfile->dataset);
		assert(currfile->file_dataspace != FAILURE);
		ret=H5Sselect_hyperslab(currfile->file_dataspace, H5S_SELECT_SET, start, stride, count, NULL);
		assert(ret != FAILURE);

			H5Eset_auto2(H5E_DEFAULT, NULL, stderr);
		/* create a memory dataspace independently */
		currfile->mem_dataspace = H5Screate_simple (ndims, count, NULL);
		/* assert (currfile->mem_dataspace != FAILURE); */
		if ( currfile->mem_dataspace == FAILURE)
		{
			currfile->mem_dataspace = H5Screate(H5S_NULL);
			if ( currfile->mem_dataspace == FAILURE)
			{
				PRINT_MSG(( SCORPIO_VERBOSE, "Mem_dataspace still not set properly. So setting H5S_ALL"));
				currfile->mem_dataspace = H5S_ALL;
			}
		}

			H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t) H5Eprint2, stderr);

		/* set up the collective transfer properties list */
		currfile->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
		assert(currfile->xfer_plist != FAILURE);

#if !defined(SERIAL_HDF5)
		ret=H5Pset_dxpl_mpio(currfile->xfer_plist, H5FD_MPIO_COLLECTIVE);
		assert(ret != FAILURE);
#endif		

		/* write data collectively */
		if (globaldims[0] > 0)
		  ret = H5Dwrite(currfile->dataset, myIOgroup->hdf_type, currfile->mem_dataspace, currfile->file_dataspace, currfile->xfer_plist, buffer);
		assert(ret != FAILURE);
		PRINT_MSG(( SCORPIO_VERBOSE, "IOnode after write"));

		/* release all temporary handles. */
		H5Sclose(currfile->file_dataspace);
		if (currfile->mem_dataspace != H5S_ALL )
			H5Sclose(currfile->mem_dataspace);
		H5Pclose(currfile->link_plist);
		H5Pclose(currfile->xfer_plist);

		/* * All writes completed.  Close datasets collectively */
		ret=H5Dclose(currfile->dataset);
		assert(ret != FAILURE);

		/* release all IDs created */
		H5Sclose(currfile->sid);

		free(buffer);

		free(start);
		free(stride);
		free(count);
		free(dimsf);

		free(slavecount);
		free(slave_nrows);
		free(displs);
		free(global_ionode_count);
		free(global_ionode_nrows);
	}

	/* Optional: In case we wish to block all slaves till the write is complete. */
	/* MPI_Barrier(myIOgroup->globalcomm); */

	PRINT_MSG(( SCORPIO_INFO, "End - Write nonuniform pattern" ));

	return 0;
}

int scorpio_write_dataset_block( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int *localstarts, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{
	void *buffer = NULL;
	int offset;
	int *slavecount = NULL;
	int localcount;
	int i,j;
	int *displs = NULL;
	int **slavedims = NULL, **slavestarts = NULL;
	int *tmpslavedims = NULL, *tmpslavestarts = NULL;
	int *slave_nrows = NULL;
	int *slave_ncols = NULL;

	PRINT_MSG(( SCORPIO_INFO, "Begin - Write block" ));
	initialize_datatype(mytype, myIOgroup);

	localcount = 1;
	// Essentially localcount = localdims[0] * localdims[1] ...
	for(i=0; i<ndims; i++)
		localcount *= localdims[i];

	if ( myIOgroup->localrank != 0)
	{
		// Note to self: remember to use & for scalar variables
		MPI_Gather(&localcount, 1, MPI_INT, slavecount, 1, MPI_INT, 0, myIOgroup->localcomm );
		MPI_Gather(localdims, ndims, MPI_INT, tmpslavedims, ndims, MPI_INT, 0, myIOgroup->localcomm );
		MPI_Gather(localstarts, ndims, MPI_INT, tmpslavestarts, ndims, MPI_INT, 0, myIOgroup->localcomm );

		/* MPI_Gather(&localdims[1], 1, MPI_INT, slave_ncols, 1, MPI_INT, 0, myIOgroup->localcomm ); */
		PRINT_MSG(( SCORPIO_VERBOSE, "Before Gatherv"));
		MPI_Gatherv(vector, localcount, myIOgroup->mpi_type, buffer, slavecount, displs, myIOgroup->mpi_type, 0, myIOgroup->localcomm );
		PRINT_MSG(( SCORPIO_VERBOSE, "After Gatherv"));
	}
	/* only do this for IO nodes */
	else
	{	 
		iofile_t *currfile;
		int local_ionode_count;
		int *global_ionode_count;

		// slavecount contains number of elements at each slave
		slavecount = (int *) calloc(myIOgroup->localsize, sizeof(int) );
		displs = (int *) calloc(myIOgroup->localsize, sizeof(int) );
		global_ionode_count = (int *) calloc(myIOgroup->numIOgroups, sizeof(int));

		tmpslavedims = (int *) calloc(myIOgroup->localsize * ndims, sizeof(int) );
		slavedims = (int **) calloc(myIOgroup->localsize, sizeof(int *) );
		tmpslavestarts = (int *) calloc(myIOgroup->localsize * ndims, sizeof(int) );
		slavestarts = (int **) calloc(myIOgroup->localsize, sizeof(int *) );

		for ( i=0; i<myIOgroup->localsize; i++)
		{
			slavedims[i] = (int *) calloc(ndims, sizeof(int));
			slavestarts[i] = (int *) calloc(ndims, sizeof(int));
		}

		// Gather all localcount from slave nodes in the local group.
		MPI_Gather(&localcount, 1, MPI_INT, slavecount, 1, MPI_INT, 0, myIOgroup->localcomm );
		MPI_Gather(localdims, ndims, MPI_INT, tmpslavedims, ndims, MPI_INT, 0, myIOgroup->localcomm );
		MPI_Gather(localstarts, ndims, MPI_INT, tmpslavestarts, ndims, MPI_INT, 0, myIOgroup->localcomm );

		/* slavedims[i][j] contains localdims in j'th dimension for slave i. */
		for ( i=0; i<myIOgroup->localsize; i++)
			for ( j=0; j<ndims; j++)
			{
				slavedims[i][j] = tmpslavedims[i*ndims+j];
				slavestarts[i][j] = tmpslavestarts[i*ndims+j];
			}

		// Calculate displacements for MPI_Gatherv
		// Each displs is equal to sum of localcounts of prior members.
		// NOTE: displs[i] is set to zero using calloc
		for(i=0; i < myIOgroup->localsize; i++)
			for ( j=0; j<i; j++)
				displs[i] += slavecount[j];

		// Calculate total number of elements at the IO node level (i.e., sum of all slave node counts)
		local_ionode_count = 0;
		for(i=0; i< myIOgroup->localsize; i++)
		{
			local_ionode_count += slavecount[i];
		}

		// Allocate buffer to hold all elements at IO node.
		buffer =  calloc(local_ionode_count, myIOgroup->datatype_size );
		if ( buffer == NULL)
		{
			PRINT_MSG(( SCORPIO_ERROR, "Failed trying to allocate %d bytes for the buffer.", local_ionode_count*myIOgroup->datatype_size));
			MPI_Abort(myIOgroup->globalcomm, 1111);
		}


		/* Gather data from all slave nodes in the local group */
		PRINT_MSG(( SCORPIO_VERBOSE, "IOnode Before Gatherv"));
		MPI_Gatherv(vector, localcount, myIOgroup->mpi_type, buffer, slavecount, displs, myIOgroup->mpi_type, 0, myIOgroup->localcomm );
		PRINT_MSG(( SCORPIO_VERBOSE, "IOnode After Gatherv"));

		// Gather all ionode_counts from all IO nodes 
		/* MPI_Allgather(&local_ionode_count, 1, MPI_INT, global_ionode_count, 1, MPI_INT, myIOgroup->iocomm ); */

		currfile = myIOgroup->file[fhandle];

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
		/* Dataset: each process takes a block. */

		// Set the dimension values from incoming argument
		for( i = 0; i < ndims; i++)
			dimsf[i] = globaldims[i];

		// Contiguous in all dimensions, so stride is 1 for all dimensions
		for( i = 0; i < ndims; i++)
			stride[i] = 1;


		PRINT_MSG(( SCORPIO_INFO, " Rank:%d  iogroup_Rank: %d start[0] %d start[1] : %d", myIOgroup->globalrank, myIOgroup->iogroupRank, start[0], start[1])); 

		PRINT_MSG(( SCORPIO_INFO, "Write - IO Node: %d - start[]=(%lu,%lu), count[]=(%lu,%lu), total datapoints=%lu", myIOgroup->iogroupRank, \
					(unsigned long)start[0], (unsigned long)start[1], \
					(unsigned long)count[0], (unsigned long)count[1], \
					(unsigned long)(count[0]*count[1]) ));


		/* check dimensions */
		for( i=0; i<ndims; i++)
			if ( dimsf[i] == 0 ) 
			{
				PRINT_MSG(( SCORPIO_INFO, " ndims: %d ", ndims));
				PRINT_MSG(( SCORPIO_INFO, " This shouldn't happen. dimsf[%d]: %d %d %d ", i, dimsf[i]));
				fprintf(stderr, "libscorpio.a: Bad argument - dimsf[%d]: %ld \n", i, dimsf[i] );
				/* dimsf[i] = 1; */
			}

		/* setup dimensionality object */
		currfile->sid = H5Screate_simple (ndims , dimsf, NULL);

		/* PRINT_MSG(( SCORPIO_INFO, "sid: %d ", currfile->sid)); */

		/* set up the link properties list so that intermediate groups are created while creating a dataset */
		currfile->link_plist = H5Pcreate (H5P_LINK_CREATE);
		assert(currfile->link_plist != FAILURE);
		H5Pset_create_intermediate_group(currfile->link_plist, 1);

		/* Supress errors while we perform our checks, reset this once we create/open the dataset */
		H5Eset_auto2(H5E_DEFAULT, NULL, stderr);

		// Open dataset collectively 
		// If this dataset is part of a existing group, use that groupid instead of file id.
		if ( currfile->groupid != -1 )
			currfile->dataset = H5Dopen2(currfile->groupid, dset_name, H5P_DEFAULT); 
		else
			currfile->dataset = H5Dopen2(currfile->fid, dset_name, H5P_DEFAULT); 

		if ( currfile->dataset == FAILURE) 
		{
			/* if dataset doesn't exist, create it */

			/* create a dataset collectively */
			if ( currfile->groupid != -1 )
				currfile->dataset = H5Dcreate2(currfile->groupid, dset_name, myIOgroup->hdf_type, currfile->sid, currfile->link_plist, H5P_DEFAULT, H5P_DEFAULT);
			else
				currfile->dataset = H5Dcreate2(currfile->fid, dset_name, myIOgroup->hdf_type, currfile->sid, currfile->link_plist, H5P_DEFAULT, H5P_DEFAULT);

			assert(currfile->dataset != FAILURE);
		}

		H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t) H5Eprint2, stderr);

		/* create a file dataspace independently */
		currfile->file_dataspace = H5Dget_space (currfile->dataset);
		assert(currfile->file_dataspace != FAILURE);

		/* set up the collective transfer properties list */
		currfile->xfer_plist = H5Pcreate (H5P_DATASET_XFER);
		assert(currfile->xfer_plist != FAILURE);

		/* Currently using collective mode */
#if !defined(SERIAL_HDF5)
		ret=H5Pset_dxpl_mpio(currfile->xfer_plist, H5FD_MPIO_COLLECTIVE);
		assert(ret != FAILURE);
#endif		

		/* Write each block from each slave in the local group consecutively */
		offset = 0;
		for (i=0; i<myIOgroup->localsize; i++)
		{

			for( j = 0; j < ndims; j++)
			{
				count[j] = slavedims[i][j]; 
				if (slavedims[i][j] == 0) 
				{
					PRINT_MSG(( SCORPIO_INFO, " Bad argument - slavedims[%d][%d] = %d ", i, j, slavedims[i][j]));
					/* fprintf(stderr, "libscorpio.a: Dataset dimension is zero - slave[%d] dims[%d] = 0 \n", i, j); */
				}
			}

			// Set offsets for each slave appropriately
			for( j = 0; j < ndims; j++)
				start[j] = slavestarts[i][j];

			ret=H5Sselect_hyperslab(currfile->file_dataspace, H5S_SELECT_SET, start, stride, count, NULL);
			if ( ret == FAILURE)
				ret = H5Sselect_all(currfile->file_dataspace);
			assert(ret != FAILURE);

			/* create a memory dataspace independently */
			currfile->mem_dataspace = H5Screate_simple (ndims, count, NULL);

			if ( currfile->mem_dataspace == FAILURE)
			{
				PRINT_MSG(( SCORPIO_INFO, "Mem_dataspace still not set properly. So setting it to NULL "));
				/* fprintf(stderr, "libscorpio.a: Mem_dataspace : setting to NULL \n"); */
				currfile->mem_dataspace = H5Screate(H5S_NULL);
				if ( currfile->mem_dataspace == FAILURE)
				{
					PRINT_MSG(( SCORPIO_INFO, "Mem_dataspace still not set properly. So setting H5S_ALL"));
					currfile->mem_dataspace = H5S_ALL;
				}
			}

			/* write one block of data collectively */
			/* the offset is measured in bytes - see line 700 below */
			ret = H5Dwrite(currfile->dataset, myIOgroup->hdf_type, currfile->mem_dataspace, currfile->file_dataspace, currfile->xfer_plist, ( (char*) buffer + offset));
			assert(ret != FAILURE);
			PRINT_MSG(( SCORPIO_VERBOSE, "IOnode after write"));

			/* release all temporary handles. */
			if ( currfile->mem_dataspace != H5S_ALL)
				H5Sclose(currfile->mem_dataspace);

			offset += slavecount[i]*myIOgroup->datatype_size;

		}

		/* * All writes completed.  Close datasets collectively */
		ret=H5Dclose(currfile->dataset);
		assert(ret != FAILURE);

		/* release all IDs created */
		if ( currfile->sid != H5S_ALL)
			H5Sclose(currfile->sid);
		H5Pclose(currfile->xfer_plist);
		H5Pclose(currfile->link_plist);

		H5Sclose(currfile->file_dataspace);

		free(buffer);

		free(start);
		free(stride);
		free(count);
		free(dimsf);

		for ( i=0; i<myIOgroup->localsize; i++)
		{
			free(slavedims[i]);
			free(slavestarts[i]);
		}

		free(slavecount);
		free(slavestarts);
		free(slavedims);
		free(tmpslavedims);
		free(tmpslavestarts);
		free(displs);
		free(global_ionode_count);
	}

	/* Optional: In case we wish to block all slaves till the write is complete. */
	/* MPI_Barrier(myIOgroup->globalcomm); */

	PRINT_MSG(( SCORPIO_INFO, "End - Write block" ));
	return 0;
}

/* primary object name can be a dataset or group or committed datatype name */
int scorpio_write_simple_attr(char *attr_name, void *attr_data, datatype_t mytype, int fhandle, char *primary_obj_name,  iogroup_t *myIOgroup)
{
	int ndims;
	int dummydims[1];
	ndims = 0;
	dummydims[0] = 1;
	/* This setups the following call to write a scalar attribute */
	return scorpio_write_attr(attr_name, attr_data, mytype, ndims, dummydims, fhandle, primary_obj_name, myIOgroup);
}

int scorpio_write_attr(char *attr_name, void *attr_data, datatype_t mytype, int ndims, int *adims, int fhandle, char *primary_obj_name,  iogroup_t *myIOgroup)
{
	int i;

	initialize_datatype(mytype, myIOgroup);
	if ( myIOgroup->localrank == 0)
	{
		/* only do this on root process - rank 0 */
		iofile_t *currfile;

		herr_t status;
		hid_t aspace_id, attr_id, loc_id;
		/* size of strings */
		size_t size;

		hsize_t *dims;

		dims = (hsize_t *) calloc(ndims, sizeof(hsize_t) );
		/* Copy incoming dims */
		for (i=0;i<ndims;i++) 
			dims[i] = adims[i];

		currfile = myIOgroup->file[fhandle];

		loc_id = H5Oopen(currfile->fid, primary_obj_name, H5P_DEFAULT);

		if ( mytype == SCORPIO_STRING)
		{
			myIOgroup->hdf_type = H5Tcopy (H5T_C_S1);
			size = strlen(attr_data);
			status = H5Tset_size(myIOgroup->hdf_type, size);
		}

		if (ndims == 0)
			aspace_id = H5Screate(H5S_SCALAR);
		else
			aspace_id = H5Screate_simple(ndims, dims, NULL); 

		attr_id = H5Acreate (loc_id, attr_name, myIOgroup->hdf_type, aspace_id, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Awrite(attr_id, myIOgroup->hdf_type, attr_data);

		status = H5Aclose(attr_id);
		status = H5Sclose (aspace_id);
		status = H5Oclose(loc_id);
		if ( mytype == SCORPIO_STRING)
			status = H5Tclose(myIOgroup->hdf_type);

		free(dims);
	}

	return 0;
}

int scorpio_write_str_array(char **strarray, int nstrings, int fhandle, char *dset_name,  iogroup_t *myIOgroup)
{
	if ( myIOgroup->globalrank == 0)
	{
		/* only do this on root process - rank 0 */
		iofile_t *currfile;
		void* buffer;
		char* pbuffer;
		herr_t ret;
		int maxstr_len, curr_len,i;

		/* We are creating a string array where each row contains one string */
		hsize_t dims[1];
		int ndims;
		ndims=1;

		dims[0] = nstrings;

		maxstr_len = INT_MIN;
		/* calculate maximum string length */
		for (i=0; i < nstrings; i++)
		{
			curr_len = strlen(strarray[i]);
			if (curr_len > maxstr_len)
				maxstr_len = curr_len;
		}
		/* Make room for null character \0 */
		maxstr_len++;
		PRINT_MSG((SCORPIO_VERBOSE, "Maximum string length: %d ", maxstr_len));

		currfile = myIOgroup->file[fhandle];
		buffer = (char *) malloc( nstrings * maxstr_len * sizeof(char));
		pbuffer = buffer;

		for (i=0; i < nstrings; i++)
		{
			sprintf(pbuffer, "%-*s", maxstr_len, strarray[i]);
			pbuffer[strlen(strarray[i])] = '\0';
			pbuffer+=maxstr_len;
			/* sprintf(pbuffer, "%s ",strarray[i]); */
		}


		/* variable length strings would require version 1.8.8+ */
		myIOgroup->hdf_type = H5Tcopy(H5T_C_S1);
		H5Tset_size(myIOgroup->hdf_type, maxstr_len);
		/* H5Tset_size(myIOgroup->hdf_type, H5T_VARIABLE); */

		/* setup dimensionality object */
		currfile->sid = H5Screate_simple (ndims , dims, NULL);

		/* set up the link properties list so that intermediate groups are created while creating a dataset */
		currfile->link_plist = H5Pcreate (H5P_LINK_CREATE);
		H5Pset_create_intermediate_group(currfile->link_plist, 1);

		/* create a dataset collectively */
		currfile->dataset = H5Dcreate2(currfile->fid, dset_name, myIOgroup->hdf_type, currfile->sid, currfile->link_plist, H5P_DEFAULT, H5P_DEFAULT);

		ret = H5Dwrite(currfile->dataset, myIOgroup->hdf_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);

		H5Pclose(currfile->link_plist);

		/* * All writes completed.  Close datasets collectively */
		ret=H5Dclose(currfile->dataset);

		/* release all IDs created */
		H5Sclose(currfile->sid);
		H5Tclose(myIOgroup->hdf_type);
		free(buffer);

	}
 
	return 0;
}


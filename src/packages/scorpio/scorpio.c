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

void initialize_datatype(datatype_t mytype, iogroup_t *myIOgroup)
{
	/* NOTE: HDF5 Native datatypes are NOT THE SAME on different platforms.  */
	/* They are similar to C type names, and are aliased to the appropriate HDF5 standard pre-defined datatype for a given platform.  */
	/* If we have cross-platform datasets, we can add them with their defined datatype accordingly. */

	/* The datatype_size variable is used for allocating temporary buffer at I/O nodes during read/write operations.  */
	/* They only need to be aware of the size of the datatype as they don't perform any other operations on the data except copy. */
	switch( mytype )
	{
		case SCORPIO_INTEGER: 
			myIOgroup->datatype_size = sizeof(int);
			myIOgroup->mpi_type = MPI_INT;
			myIOgroup->hdf_type = H5T_NATIVE_INT;
			break;
		case SCORPIO_DOUBLE:
			myIOgroup->datatype_size = sizeof(double);
			myIOgroup->mpi_type = MPI_DOUBLE;
			myIOgroup->hdf_type = H5T_NATIVE_DOUBLE;
			break;
		case SCORPIO_FLOAT: 
			myIOgroup->datatype_size = sizeof(float);
			myIOgroup->mpi_type = MPI_FLOAT;
			myIOgroup->hdf_type = H5T_NATIVE_FLOAT;
			break;
		case SCORPIO_LONG: 
			myIOgroup->datatype_size = sizeof(long);
			myIOgroup->mpi_type = MPI_LONG;
			myIOgroup->hdf_type = H5T_NATIVE_LONG;
			break;
		case SCORPIO_CHAR: 
			myIOgroup->datatype_size = sizeof(char);
			myIOgroup->mpi_type = MPI_CHAR;
			myIOgroup->hdf_type = H5T_NATIVE_CHAR;
			break;
		case SCORPIO_BYTE: 
			myIOgroup->datatype_size = 1; /* 1 byte */
			myIOgroup->mpi_type = MPI_BYTE;
			/* uninterpreted bytes */
			myIOgroup->hdf_type = H5T_NATIVE_OPAQUE;
			break;
	}
}

/**
 * @brief Parallel IO framework initialization
 * Splits the given MPI communicator into requisite number of groups.
 * A configuration structure is used to represent incoming configuration to avoid function signature changes over time.
 * @param pf_conf_ptr pointer to the parallel framework configuration
 */ 
int scorpio_IOgroup_init(iogroup_conf_t *myIOgroupConf /* in */, iogroup_t *myIOgroup)
{
	int i, ierr;
	MPI_Group globalMPIgroup, ioMPIgroup; 
	/* List of IO nodes */
	int *ioNodes;

	/* Temp variables for debugging */
	int ioNodeRank, numIOnodes;

	ierr = MPI_Comm_dup(myIOgroupConf->commIncoming, &myIOgroup->globalcomm );
	if ( ierr != MPI_SUCCESS)
	{
		PRINT_MSG(( SCORPIO_ERROR, "MPI_Comm_dup failed."));
		exit(101);
	}


	MPI_Comm_rank(myIOgroup->globalcomm, &myIOgroup->globalrank);
	MPI_Comm_size(myIOgroup->globalcomm, &myIOgroup->globalsize);

	/* initialize debugging log (essentially passing in globalrank info) */
	scorpio_debuglog_init(myIOgroup->globalrank);

	if ( myIOgroupConf->numIOgroups > 0 )
	{
		if (myIOgroup->globalrank == 0 )
		{	
			/* fprintf(stderr, "SCORPIO_Info: You are currently using numIOgroups. Alternately you can use preferredGroupSize.\n"); */

			/* if ( myIOgroupConf->preferredGroupSize == 0) */
				/* fprintf(stderr, "SCORPIO_Info: Preferred group size is not set\n"); */

			if ( myIOgroup->globalsize % myIOgroupConf->numIOgroups != 0)
			{
				fprintf(stderr, "SCORPIO_Info: nprocs is not exactly divisible by numIOgroups\n");
				fprintf(stderr, "SCORPIO_Info: numIOgroups will be one more than requested, the last group containing the remaining processes. \n");
			}
			fprintf(stderr, "SCORPIO_Info: Preferred group size is set to nprocs/numIOgroups (%d/%d) \n", myIOgroup->globalsize, myIOgroupConf->numIOgroups);

		}

		myIOgroup->preferredGroupSize = myIOgroup->globalsize / myIOgroupConf->numIOgroups;

		/* Manually set preferredGroupSize in configuration struct too for subsequent calculations */
		myIOgroupConf->preferredGroupSize = myIOgroup->preferredGroupSize;

	}
	else
	{
		/* Just copy value from incoming configuration object (set by user program) */
		myIOgroup->preferredGroupSize = myIOgroupConf->preferredGroupSize;
	}

	if ( myIOgroupConf->preferredGroupSize <= 0 || myIOgroupConf->preferredGroupSize > myIOgroup->globalsize)
	{
		if (myIOgroup->globalrank == 0 )
		{
			fprintf(stderr, "SCORPIO_Error: Preferred group size should be > 0 and < number of processors. Given size = %d\n",myIOgroupConf->preferredGroupSize);
			fprintf(stderr, "SCORPIO_Info: Setting Preferred group size to total number of procs to continue execution \n");
		}

		/* Set group size to total number of processes */
		/* Manually set preferredGroupSize in configuration struct too for subsequent calculations */
		myIOgroupConf->preferredGroupSize = myIOgroup->preferredGroupSize = myIOgroup->globalsize ;

		/* Sarat:  Avoiding aborting the application by choosing sane defaults in case of user error */
		/* MPI_Abort(myIOgroup->globalcomm, 1111); */
	}

	/* Initially we are working under the assumption that globalsize is exact multiple of numIOgroups */
	if (myIOgroup->globalsize % myIOgroupConf->preferredGroupSize == 0 )
	{
		myIOgroup->numIOgroups = myIOgroup->globalsize/myIOgroupConf->preferredGroupSize;

		/* myIOgroup->iogroupSize = myIOgroup->globalsize / myIOgroup->numIOgroups; */
		myIOgroup->iogroupSize = myIOgroupConf->preferredGroupSize;
		myIOgroup->iogroupRank = myIOgroup->globalrank / myIOgroup->preferredGroupSize;

	}
	else
	{
		fprintf(stderr, "SCORPIO_Info: nprocs is not exactly divisible by preferredGroupSize\n");

		/* Last group has fewer processes than the rest  */
		myIOgroup->numIOgroups = myIOgroup->globalsize/myIOgroupConf->preferredGroupSize + 1;

		myIOgroup->iogroupRank = myIOgroup->globalrank / myIOgroup->preferredGroupSize;

		/* iogroupRanks start from zero, so the final group would be n-1 */
		if( myIOgroup->iogroupRank == (myIOgroup->numIOgroups-1) )
			/* last group has remainder of processes*/
			myIOgroup->iogroupSize = myIOgroup->globalsize % myIOgroupConf->preferredGroupSize;
		else
			myIOgroup->iogroupSize = myIOgroupConf->preferredGroupSize;

	}
	PRINT_MSG(( SCORPIO_INFO, "numIOgroups: %d, iogroupSize: %d, iogroupRank: %d", myIOgroup->numIOgroups, myIOgroup->iogroupSize, myIOgroup->iogroupRank));


	if (myIOgroup->globalrank == 0 )
	{	
		fprintf(stderr, "SCORPIO_Info: Preferred group size is set to %d \n",myIOgroup->preferredGroupSize);
	}


	ierr = MPI_Comm_split(myIOgroup->globalcomm, myIOgroup->iogroupRank, myIOgroup->globalrank, &myIOgroup->localcomm); 

	if ( ierr != MPI_SUCCESS)
	{
		PRINT_MSG(( SCORPIO_ERROR, "MPI Subcommunicator creation failed."));
		exit(102);
	}

	MPI_Comm_rank(myIOgroup->localcomm, &myIOgroup->localrank);
	MPI_Comm_size(myIOgroup->localcomm, &myIOgroup->localsize);

	/* Create a communicator just for the IO nodes */
	if (myIOgroup->localrank == 0)
		ierr = MPI_Comm_split(myIOgroup->globalcomm, 1 , myIOgroup->globalrank, &myIOgroup->iocomm); 
	else
		ierr = MPI_Comm_split(myIOgroup->globalcomm, MPI_UNDEFINED, myIOgroup->globalrank, &myIOgroup->iocomm); 

	if ( ierr != MPI_SUCCESS)
	{
		PRINT_MSG(( SCORPIO_ERROR, "MPI Subcommunicator creation failed."));
		exit(102);
	}

	/*
	 * OLD METHOD TO CREATE INTER-GROUP COMMUNICATOR
	 *
	MPI_Comm_group(myIOgroup->globalcomm, &globalMPIgroup);

	ioNodes = (int *) malloc(sizeof(int)*myIOgroup->numIOgroups);
	for (i = 0; i < myIOgroup->numIOgroups; i++)
		ioNodes[i] = i*myIOgroup->iogroupSize;

	MPI_Group_incl(globalMPIgroup, myIOgroup->numIOgroups, ioNodes, &ioMPIgroup);
	MPI_Comm_create(myIOgroup->globalcomm, ioMPIgroup, &myIOgroup->iocomm);     
	free(ioNodes);
	*/

	// Initially there are no open files.
	myIOgroup->numFiles = 0; 
	// initially set the file pointer to NULL
	myIOgroup->file = NULL;

	PRINT_MSG(( SCORPIO_INFO, "Local rank: %d, Global Rank: %d, IO Rank: %d, \n Local size: %d, Global Size: %d", myIOgroup->localrank, myIOgroup->globalrank, myIOgroup->iogroupRank, myIOgroup->localsize, myIOgroup->globalsize)); 

	if ( myIOgroup->localrank == 0 )
	{ 
		MPI_Comm_rank(myIOgroup->iocomm, &ioNodeRank); 
		MPI_Comm_size(myIOgroup->iocomm, &numIOnodes);
		PRINT_MSG((SCORPIO_VERBOSE, "IO node rank: %d, IO comm size: %d ", ioNodeRank, numIOnodes));
	}

	return 0;
	PRINT_MSG(( SCORPIO_VERBOSE, "Parallel I/O Init done.")); 

}


/**
 * Cleanup of parallel framework.
 */
int scorpio_IOgroup_cleanup(iogroup_t *myIOgroup)
{

	int i;

	/* Free IO communicator only on IO nodes */
	if (myIOgroup->localrank == 0)
	{
		for( i =0; i < myIOgroup->numFiles; i++)
			free(myIOgroup->file[i]);
		free(myIOgroup->file);
		MPI_Comm_free(&myIOgroup->iocomm);
	}

	MPI_Comm_free(&myIOgroup->globalcomm);
	MPI_Comm_free(&myIOgroup->localcomm);

	PRINT_MSG(( SCORPIO_VERBOSE, "Parallel I/O group cleanup done.")); 

	return 0;
}


int scorpio_open_file(const char *filename, iogroup_t *myIOgroup, file_mode_t mode  )
{
	herr_t ret;         	/* Generic return value */

	/* only do this for IO nodes */
	if ( myIOgroup->localrank == 0)
	{	 
		int fhandle;
		iofile_t *currfile;

		// Update counter for the number of files associated with this IO group.
		myIOgroup->numFiles++;
		// Create a placeholder for the new file.
		myIOgroup->file = (iofile_t **) realloc( myIOgroup->file,  myIOgroup->numFiles * sizeof(iofile_t *) );

		// Handle indices start from 0. So reduce 1 from the current number of files to get the file handle.
		fhandle = myIOgroup->numFiles - 1 ;
		myIOgroup->file[fhandle] = (iofile_t *) calloc( 1 ,  sizeof(iofile_t) );

		currfile = myIOgroup->file[fhandle];

		/* just storing this handle in the file struct again for bookkeeping */
		currfile->handle = fhandle; 

		strncpy(currfile->filename, filename, (long) ( strlen(filename) < MAX_FILENAME) ? strlen(filename) : MAX_FILENAME ); 


		MPI_Info info = MPI_INFO_NULL;

		/* Initially set the groupid to -1, only useful when intermediate dataset groups are created explicitly */
		/* currently useful while writing block data in Pflotran, needed because of the way that codebase is structured */
		currfile->groupid = -1;

		/* setup file access template with parallel IO access. */
		/* hid_t acc_tpl;		 */
		currfile->acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
		assert(currfile->acc_tpl != FAILURE);
		PRINT_MSG(( SCORPIO_VERBOSE, "H5Pcreate access succeed"));
		/* set Parallel access with communicator */
#if !defined(SERIAL_HDF5)
		ret = H5Pset_fapl_mpio(currfile->acc_tpl, myIOgroup->iocomm, info);
		assert(ret != FAILURE);
		PRINT_MSG(( SCORPIO_VERBOSE, "H5Pset_fapl_mpio succeed"));
#endif		

		if ( mode == SCORPIO_FILE_CREATE)
		{
			/* create the file collectively */
			currfile->fid = H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,currfile->acc_tpl);
			/* assert(currfile->fid != FAILURE); */
			/* PRINT_MSG(( SCORPIO_VERBOSE, "H5Fcreate succeed")); */
		}
		else if (mode == SCORPIO_FILE_READWRITE)
		{
			currfile->fid = H5Fopen(filename, H5F_ACC_RDWR, currfile->acc_tpl);
			/* assert(currfile->fid != FAILURE); */
			/* PRINT_MSG(( SCORPIO_VERBOSE, "H5Fopen succeed")); */
		}
		else
		{
			// H5F_ACC_RDONLY
			currfile->fid = H5Fopen(filename,H5F_ACC_RDONLY, currfile->acc_tpl);
			/* assert(currfile->fid != FAILURE); */
			/* PRINT_MSG(( SCORPIO_VERBOSE, "H5Fopen succeed")); */
		}

		/* Release file-access template */
		ret=H5Pclose(currfile->acc_tpl);
		assert(ret != FAILURE);

		/* Check open error code */
		if ( currfile->fid < 0)
			ret = FAILURE;
		else
			ret = fhandle; // Return current index as the handle
	}

	/* Send return value (fhandle or error code from 'open' call above) */
	MPI_Bcast(&ret, 1, MPI_INT,  0, myIOgroup->localcomm);

	return ret;
}

int scorpio_close_file(int fhandle, iogroup_t *myIOgroup )
{

	/* only do this for IO nodes */
	if ( myIOgroup->localrank == 0)
	{	 
		iofile_t *currfile;
		currfile = myIOgroup->file[fhandle];

		/* close the file collectively */
		H5Fclose(currfile->fid);
	}

	return 0;

}

/** 
 * @brief Check if a HDF5 dataset group exists in the given file. 
 * Assumes that the file is already open 
 * @param group_name
 * @param fhandle file handle
 * @param myIOgroup IOgroup used to open this file
 * 
 * @return 0 if the group doesn't exist, 1 if the group exists
 */
int scorpio_group_exists(char *group_name, int fhandle, iogroup_t *myIOgroup )
{
	int flag;

	/* One process can check for everyone */
	if ( myIOgroup->globalrank == 0)
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
			/* group doesn't exist, so return 0 */
			flag = 0;
		}
		else
		{
			/* Group exists, so return 1 */
			flag = 1;
			H5Gclose(currfile->groupid);
		}

		H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t) H5Eprint2, stderr);
	}

	/* printf("Checking if group %s exists : Flag: %d \n", group_name, flag); */

	/* Send return value to everyone */
	MPI_Bcast(&flag, 1, MPI_INT,  0, myIOgroup->globalcomm);
	PRINT_MSG((SCORPIO_VERBOSE, "Checking if group %s exists : Flag: %d \n", group_name, flag)); 
	return flag;
}

/** 
 * @brief Check if a HDF5 dataset exists in the given file. 
 * Assumes that the file is already open 
 * @param dataset_name
 * @param fhandle file handle
 * @param myIOgroup IOgroup used to open this file
 * 
 * @return 0 if the dataset doesn't exist, 1 if the group exists
 */
int scorpio_dataset_exists(char *dataset_name, int fhandle, iogroup_t *myIOgroup )
{
	int flag;

	/* One process can check for everyone */
	if ( myIOgroup->globalrank == 0)
	{
		iofile_t *currfile;
		currfile = myIOgroup->file[fhandle];

		H5Eset_auto2(H5E_DEFAULT, NULL, stderr);

		/* Check if dataset exists by trying to open it */
		currfile->dataset = H5Dopen2(currfile->fid, dataset_name, H5P_DEFAULT); 
		if (currfile->dataset == FAILURE)
		{
			/* dataset doesn't exist, so return 0 */
			flag = 0;
		}
		else
		{
			/* dataset exists, so return 1 */
			flag = 1;
			H5Dclose(currfile->dataset);
		}

		H5Eset_auto2(H5E_DEFAULT, (H5E_auto2_t) H5Eprint2, stderr);
	}

	/* printf("Checking if dataset %s exists : Flag: %d \n", dataset_name, flag); */

	/* Send return value to everyone */
	MPI_Bcast(&flag, 1, MPI_INT,  0, myIOgroup->globalcomm);
	PRINT_MSG((SCORPIO_VERBOSE, "Checking if dataset %s exists : Flag: %d \n", dataset_name, flag)); 
	return flag;
}


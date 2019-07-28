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
#include "scorpioftn.h"

iogroup_t **iogroups = NULL;
int numIOgroups = 0;

/* Newer init interface declaration */
void F2C(scorpio_iogroup_init2) ( int *ppreferredGroupSize, int *pnumSubgroups, void *mpicomm_in, int *pgid, int *ierr);

/* Original init call interface */
void F2C(scorpio_iogroup_init) ( int *pnumSubgroups, void *mpicomm_in, int *pgid, int *ierr)
{
	/* printf("Got communicator: %d\n", *(int *)mpicomm_in); */
	int dummysize = 0;
	F2C(scorpio_iogroup_init2) ( &dummysize, pnumSubgroups, mpicomm_in, pgid, ierr);
}

/* Init call interface supporting preferredGroupSize */
void F2C(scorpio_iogroup_init2) ( int *ppreferredGroupSize, int *pnumSubgroups, void *mpicomm_in, int *pgid, int *ierr)
{
	iogroup_conf_t myIOconfig;
	int mpi_rank;
	int curr_gid;
	int myerr;

	setvbuf(stdout, (char *) NULL, _IONBF, 0);

	/* One of the following should be set. */

	if ( *pnumSubgroups > 0) 
		myIOconfig.numIOgroups = *pnumSubgroups;

	if ( *ppreferredGroupSize > 0) 
		myIOconfig.preferredGroupSize = *ppreferredGroupSize;

	if ( *pnumSubgroups == 0 && ppreferredGroupSize == 0) 
	{
		fprintf(stderr, "SCORPIO_Error: Parallel IO group initialization requires at least one parameter set (numIOgroups or preferredGroupSize)\n");
		exit(101);
	}

	/* printf("Got communicator: %d\n", *(int *)mpicomm_in); */
	/* printf("C MPI_COMM_WORLD: %d\n", MPI_COMM_WORLD); */
	/* printf("f2c(0) where 0=fortran MPI_COMM_WORLD in openmpi: %d\n", MPI_Comm_f2c(0) ); */
	/* printf("f2c incoming communicator: %d\n", MPI_Comm_f2c( *(int *)mpicomm_in) ); */

	/* You cannot use or assign incoming mpi communicator directly because OpenMPI uses a custom struct unlike MPICH. 
	 * So that caused problems with openmpi but worked with MPICH 
	 * Following two lines will work with MPICH but not with OpenMPI. So use Comm_f2c to get the right communicator */
	/* myerr = MPI_Comm_dup( * ((MPI_Comm *)mpicomm_in) , &(myIOconfig.commIncoming)); */
	/* myIOconfig.commIncoming = * ((MPI_Comm *)mpicomm_in); */
	myerr = MPI_Comm_dup( MPI_Comm_f2c( *(MPI_Fint *)mpicomm_in) , &(myIOconfig.commIncoming));
	if ( myerr != MPI_SUCCESS)
	{
		PRINT_MSG(( SCORPIO_ERROR, "MPI_Comm_dup failed."));
		exit(102);
	}
	/* printf("Created communicator: %d\n", myIOconfig.commIncoming); */

	MPI_Comm_rank(myIOconfig.commIncoming, &mpi_rank);
	scorpio_debuglog_init(mpi_rank);

	/* indexing starts from 0 */
	curr_gid = numIOgroups;
	/* send this back to Fortran calling function */
	*pgid = curr_gid;

	PRINT_MSG((SCORPIO_INFO, "Init Subgroups: %d \t MPI Comm: %d \t IOgroup ID: %d ", *pnumSubgroups, * (int *)mpicomm_in, *pgid));
	numIOgroups++;

	/* Always make sure this pointer is initialized to NULL the first time it is called */
	iogroups = (iogroup_t **) realloc( iogroups, numIOgroups * sizeof(iogroup_t *));
	iogroups[curr_gid] = (iogroup_t *) calloc(1, sizeof(iogroup_t));

	*ierr = scorpio_IOgroup_init(&myIOconfig, iogroups[curr_gid]);

}

void F2C(scorpio_cleanup) ( int *ierr)
{
	PRINT_MSG((SCORPIO_INFO, "Cleaning up residual IO groups pointer"));
	free(iogroups);
	*ierr = 0;
}


void F2C(scorpio_iogroup_cleanup) ( int *pgid, int *ierr)
{
	PRINT_MSG((SCORPIO_INFO, "Cleaning up IO group: %d ", *pgid));

	*ierr = scorpio_IOgroup_cleanup(iogroups[*pgid]);

	free(iogroups[*pgid]);
}

void F2C(print_string) ( char *filename, int *x, int *ierr, int strlen)
{
	fprintf(stdout, "Inside print_string"); fflush(stdout);
	fprintf(stdout, "\n String in C: %s Strlen: %d \n", filename, strlen);
	fflush(stdout);
	fprintf(stdout, "Extra arg: %d\n", *x);
	fflush(stdout);
	*ierr = 0;
}

void F2C(print_array2) ( void *xin, int *ierr)
{
	int *x = xin;

	fprintf(stdout, "Address: %p \n", x) ; fflush(stdout);
	fprintf(stdout, "Value: %d \n", *x) ; fflush(stdout);

	fprintf(stdout, "Vector: %d %d\n", x[0], x[1]);
	/* fprintf(stdout, "Vector(double): %f %f\n", x[0], x[1]); */
	fflush(stdout);
	*ierr = 0;
}

void F2C(print_array) ( int *x, int *ierr)
{
	fprintf(stdout, "Address: %p \n", x) ; fflush(stdout);
	fprintf(stdout, "Value: %d \n", *x) ; fflush(stdout);

	fprintf(stdout, "Vector: %d %d\n", x[0], x[1]);
	fflush(stdout);
	*ierr = 0;
}


/* Note: Strlen argument for a character string has to go at the end when interfacing with Fortran */
void F2C(scorpio_open_file) ( char *filename, int *pgid, int *filemode, int *pfid, int *ierr, int strlen)
/* void F2C(scorpio_open_file) ( int *pgid, int *filemode, int *pfid, char *filename, int strlen) */
{

	PRINT_MSG((SCORPIO_INFO, "Inside File open."));
	PRINT_MSG((SCORPIO_INFO, "Filename: %s Strlen: %d ", filename, strlen));
	PRINT_MSG((SCORPIO_INFO, "GroupID: %d Filemode: %d ", *pgid, *filemode));

	PRINT_MSG((SCORPIO_INFO, "Globalrank: %d Localrank: %d ", iogroups[*pgid]->globalrank, iogroups[*pgid]->localrank ));

	*pfid = scorpio_open_file(filename, iogroups[*pgid], (file_mode_t) *filemode);

	PRINT_MSG((SCORPIO_INFO, "File open finished."));
	if ( *pfid == FAILURE )
		*ierr = FAILURE;
	else	
		*ierr = 0;
}
void F2C(scorpio_close_dataset_group)( int *groupid, int *pfid, int *pgid, int *ierr)
{
	*ierr = scorpio_close_dataset_group(*groupid, *pfid, iogroups[*pgid]);
}

/** 
 * @brief Check if a HDF5 dataset group exists in the given file. 
 * Assumes that the file is already open
 * @param group_name
 * @param pfid file id 
 * @param pfid ptr to id of IOgroup used to open this file
 * @return 0 if the group doesn't exist, 1 if the group exists
 */
void F2C(scorpio_group_exists)( char *group_name, int *pfid, int *pgid, int *ierr, int strlen)
{
	*ierr = scorpio_group_exists(group_name, *pfid, iogroups[*pgid]);
}

/** 
 * @brief Check if a HDF5 dataset group exists in the given file. 
 * Assumes that the file is already open
 * @param group_name
 * @param pfid file id 
 * @param pfid ptr to id of IOgroup used to open this file
 * @return 0 if the group doesn't exist, 1 if the group exists
 */
void F2C(scorpio_dataset_exists)( char *dataset_name, int *pfid, int *pgid, int *ierr, int strlen)
{
	*ierr = scorpio_dataset_exists(dataset_name, *pfid, iogroups[*pgid]);
}

void F2C(scorpio_create_dataset_group)( int *mygroupid, char *group_name, int *pfid, int *pgid, int *ierr, int strlen)
/* int scorpio_create_dataset_group( char *group_name, int fhandle, iogroup_t *myIOgroup) */
{
	*mygroupid = scorpio_create_dataset_group(group_name, *pfid, iogroups[*pgid]);
	if ( *mygroupid == FAILURE)
		*ierr = FAILURE;
	else 
		*ierr = 0;
}

void F2C(scorpio_close_file) ( int *pfid, int *pgid, int *ierr)
{
	*ierr = scorpio_close_file(*pfid, iogroups[*pgid]);
}


void F2C(scorpio_read_dataset)( void *vector, int *pmytype, int *pndims, int *globaldims, int *localdims, int *pfhandle, char *dset_name, int *pgid, int *pmypattern, int *ierr, int strlen)
{
	PRINT_MSG((SCORPIO_INFO, "Type: %d NDims: %d Fhandle: %d Dataset: %s IOgroup: %d IOpattern: %d", \
				*pmytype, *pndims, *pfhandle, dset_name, *pgid, *pmypattern));
					  
	*ierr = scorpio_read_dataset(vector, (datatype_t) *pmytype, *pndims, globaldims, localdims, *pfhandle, dset_name, iogroups[*pgid], (iopattern_t) *pmypattern);
}

void F2C(scorpio_write_dataset)( void *vector, int *pmytype, int *pndims, int *globaldims, int *localdims, int *pfhandle, char *dset_name, int *pgid, int *pmypattern, int *ierr, int strlen)
{
	PRINT_MSG((SCORPIO_INFO, "Type: %d NDims: %d Fhandle: %d Dataset: %s IOgroup: %d IOpattern: %d", \
				*pmytype, *pndims, *pfhandle, dset_name, *pgid, *pmypattern));
					  
	*ierr = scorpio_write_dataset(vector, (datatype_t) *pmytype, *pndims, globaldims, localdims, *pfhandle, dset_name, iogroups[*pgid], (iopattern_t) *pmypattern);
}

void F2C(scorpio_write_dataset_block)( void *vector, int *pmytype, int *pndims, int *globaldims, int *localdims, int *localstarts, int *pfhandle, char *dset_name, int *pgid, int *ierr, int strlen)
{
	PRINT_MSG((SCORPIO_INFO, "Type: %d NDims: %d Fhandle: %d Dataset: %s IOgroup: %d Block IOpattern", \
				(datatype_t) *pmytype, *pndims, *pfhandle, dset_name, *pgid));
					  
	PRINT_MSG((SCORPIO_INFO, "Block IOpattern ndims: %d globaldims: %d %d %d ", *pndims, globaldims[0], globaldims[1], globaldims[2]));
	PRINT_MSG((SCORPIO_INFO, "Block IOpattern ndims: %d localdims: %d %d %d ", *pndims, localdims[0], localdims[1], localdims[2]));

	*ierr = scorpio_write_dataset_block(vector, (datatype_t) *pmytype, *pndims, globaldims, localdims, localstarts, *pfhandle, dset_name, iogroups[*pgid] );
}

void F2C(scorpio_get_dataset_ndims)( int *pnumdims, int *pfhandle, char *dset_name,  int *pgid, int *ierr, int strlen)
{
	*ierr = scorpio_get_dataset_ndims(pnumdims, *pfhandle, dset_name, iogroups[*pgid]);
}

void F2C(scorpio_get_dataset_dims)( int *pdims, int *pfhandle, char *dset_name,  int *pgid, int *ierr, int strlen)
{
	*ierr = scorpio_get_dataset_dims(pdims, *pfhandle, dset_name, iogroups[*pgid]);
}

void F2C(scorpio_get_dataset_size)( int *pdataset_size, int *pfhandle, char *dset_name,  int *pgid, int *ierr, int strlen)
{

	*ierr = scorpio_get_dataset_size(pdataset_size, *pfhandle, dset_name, iogroups[*pgid]);
}

void F2C(scorpio_read_same_sub_dataset)( void *vector, int *pmytype, int *pndims, int *localdims, int *localstarts, int *pfhandle, char *dset_name, int *pgid, int *ierr, int strlen)
{
	PRINT_MSG((SCORPIO_INFO, "Type: %d NDims: %d Fhandle: %d Dataset: %s IOgroup: %d IOpattern: Everyone reads same part of dataset", \
				*pmytype, *pndims, *pfhandle, dset_name, *pgid));
					  
	*ierr = scorpio_read_same_sub_dataset(vector, (datatype_t) *pmytype, *pndims, localdims, localstarts, *pfhandle, dset_name, iogroups[*pgid]);
}


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
 * @file strex.c
 * @brief Driver for testing string support in parallel I/O library
 * @author Sarat Sreepathi
 * @version 1.0
 * @date 2012-05-25
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "scorpio.h"


#define MAX_STR_LEN 100

int main(int argc, char **argv)
{
	int errcode;
	int fhandle;
	char filename[80] = "strings.h5";
	char tmpstr[80];
	int i,j;
	int nstrings;
	char **strArray;
	int nprocs;
	int rank;

	iogroup_conf_t myIOconfig;
	iogroup_t myIOgroup;

	errcode = MPI_Init(&argc, &argv);
	if (errcode != MPI_SUCCESS) { fprintf(stderr, "Error initializing MPI\n"); exit(2); }

	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	myIOconfig.numIOgroups = 1;
	myIOconfig.commIncoming = MPI_COMM_WORLD;

	scorpio_IOgroup_init(&myIOconfig, &myIOgroup);

	nstrings = 5;
	strArray = (char **) malloc ( nstrings * sizeof(char *));

	for (i=0; i < nstrings ; i++)
		strArray[i] = (char *) malloc( MAX_STR_LEN * sizeof(char) );	

	sprintf(strArray[0],"Calcium");
	sprintf(strArray[1],"Magnesium");
	sprintf(strArray[2],"Uranium");
	sprintf(strArray[3],"Unobtainium");
	sprintf(strArray[4],"My Favorite Mineral in the Whole World");


	/* --------------------------------------------------------------------------- */
	/* Parallel I/O Calls  */
	/* --------------------------------------------------------------------------- */

	fhandle = scorpio_open_file(filename, &myIOgroup, SCORPIO_FILE_CREATE);
	if (fhandle == -1)
	{
		fprintf(stderr, "Opening file for writing failed.\n");
		exit(-1);
	}

	/* int scorpio_write_str_array(char *strarray[], int nstrings, int fhandle, char *dset_name,  iogroup_t *myIOgroup) */
	scorpio_write_str_array(strArray, nstrings, fhandle, "Mineral names", &myIOgroup);

	scorpio_close_file( fhandle, &myIOgroup);
	/* free memory before reading back data */
	free(strArray);

	/* Open file for reading now */
	fhandle = scorpio_open_file(filename, &myIOgroup, SCORPIO_FILE_READONLY);
	if (fhandle == -1)
	{
		fprintf(stderr, "Opening file for reading failed.\n");
		exit(-1);
	}

	scorpio_read_str_array(&strArray, &nstrings, fhandle, "Mineral names", &myIOgroup);

	for (i=0; i < nstrings; i++)
		printf("R%d-String[%d]: %s$\n", rank, i, strArray[i]);

	free(strArray);
	scorpio_close_file( fhandle, &myIOgroup);

	/* Cleanup */
	scorpio_IOgroup_cleanup(&myIOgroup);
	MPI_Finalize();

}



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
 * @file scorpio.h 
 * @brief C Header file for the parallel I/O framework interface.
 * @author Sarat Sreepathi (sarat@computer.org)
 * @version 0.01
 * @date 2010-09-15
 */

#ifndef _PARALLELIO_H
#define _PARALLELIO_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>	
#include <hdf5.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>

#include "scorpio_util.h"

#define MAX_FILENAME 100	

#define FAILURE -1

	/** 
	 * Struct containing configuration(input) for a parallel IO group
	 * A configuration structure is used to represent incoming configuration to avoid function signature changes over time.
	 */	
	typedef struct
	{
		/** Number of IO groups to partition the global MPI communicator set */
		int numIOgroups;

		/** This is the preferred group size for each I/O group. 
		 * The user is recommended to use this instead of numIOgroups(deprecated now).
		 * If the total number of processes is not exactly divisible by this number, 
		 * the last group would have fewer than preferredGroupSize processes in it */
		int preferredGroupSize;

		/** Incoming MPI communicator that is used to split into IO groups */
		MPI_Comm commIncoming;
		/** Placeholder in case we wish to use custom group sizes instead of an uniform IO group size */
		// int *groupsizes;

		/** Additional parameters that can be used to specify initial configuration. */
	}
	iogroup_conf_t;

	/** 
	 * @brief IO File structure
	 */
	typedef struct
	{
		int handle;
		char filename[MAX_FILENAME];
		hid_t fid;			/* HDF5 file ID */
		hid_t acc_tpl;		/* File access templates */
		hid_t xfer_plist;	/* Dataset transfer properties list */
		hid_t link_plist;	/* Link properties list (used to create intermediate groups when creating a dataset */
		hid_t sid;   		/* Dataspace ID */
		hid_t file_dataspace;	/* File dataspace ID */
		hid_t mem_dataspace;	/* memory dataspace ID */
		hid_t dataset;	/* Dataset ID */
		hid_t groupid;
	}
	iofile_t; 


	/** 
	 * Struct containing state of a parallel I/O group 
	 * Nice placeholder encapsulating all the requisite data
	 * without heavily polluting the global namespace.
	 */	
	typedef struct
	{
		/** Rank in the local group */
		int localrank; 
		/** Rank in the global group */
		int globalrank;
		/** Number of processes in local IO group */
		int localsize; 
		/** Number of processes in global group */
		int globalsize;
		/** ID of the current group this process belongs to */
		int iogroupRank; 
		/** Size of the current IO group (limited use currently,can be used in future for different group sizes etc) */
		int iogroupSize;
		/** Number of groups to partition the global MPI communicator set */
		int numIOgroups;

		/** This is the preferred group size for each I/O group. 
		 * The user is recommended to use this instead of numIOgroups(deprecated now).
		 * If the total number of processes is not exactly divisible by this number, 
		 * the last group would have fewer than preferredGroupSize processes in it */
		int preferredGroupSize;

		/** MPI communicators for the local group and global group */
		MPI_Comm localcomm, globalcomm, iocomm;
		/** Files opened using this IO group */
		iofile_t **file;
		/** Number of files used with this iogroup, this is a monotonically increasing counter that just keeps track of all files opened with this group */
		int numFiles;

		// Temporary variables to store current dataset datatype and its corresponding MPI and HDF5 datatypes
		int datatype_size;
		MPI_Datatype mpi_type;
		hid_t hdf_type;

	}
	iogroup_t;

	typedef enum
	{
		SCORPIO_UNIFORM_CONTIGUOUS_READ, /* Each PE reads the same contiguous amount of data */
		SCORPIO_NONUNIFORM_CONTIGUOUS_READ,
		SCORPIO_NONCONTIGUOUS_READ,      /* Each PE reads noncontiguous data ... */
		SCORPIO_EVERYONE_ENTIRE_DATASET_READ,
		SCORPIO_UNIFORM_CONTIGUOUS_WRITE,
		SCORPIO_NONUNIFORM_CONTIGUOUS_WRITE,
		SCORPIO_NONCONTIGUOUS_WRITE
	} iopattern_t; 

	typedef enum 
	{
		SCORPIO_FILE_CREATE,
		SCORPIO_FILE_READONLY,
		SCORPIO_FILE_READWRITE
	} file_mode_t ;

	/* The prefix PIO refers to Parallel I/O, prefix used to avoid conflicts with any existing constants */
	typedef enum 
	{
		SCORPIO_INTEGER,
		SCORPIO_DOUBLE,
		SCORPIO_FLOAT,
		SCORPIO_LONG,
		SCORPIO_CHAR,
		SCORPIO_BYTE,
		SCORPIO_STRING
	} datatype_t ;

	/** 
	 * @brief Function to initialize type info based on user input
	 *  This is used to customize the read/write calls based on current datatype info.
	 * @param mytype Current datatype (INTEGER or DOUBLE)
	 * @param myIOgroup current IOgroup
	 */
	void initialize_datatype(datatype_t mytype, iogroup_t *myIOgroup);

	/** 
	 * Externally callable functions
	 */

	/** 
	 * @brief Parallel IO framework initialization
	 * Splits the given MPI communicator into requisite number of groups.
	 * A configuration structure is used to represent incoming configuration to avoid function signature changes over time.
	 * 
	 * @param myIOgroupConf pointer to the IO group configuration
	 * @param myIOgroup
	 * 
	 * @return error code
	 */
	int scorpio_IOgroup_init(iogroup_conf_t *myIOgroupConf /* in */, iogroup_t *myIOgroup);

	/** 
	 * @brief Cleans up I/O group and any internal datastructures.
	 * @param myIOgroup IOgroup to clean up
	 * @return error code
	 */
	int scorpio_IOgroup_cleanup(iogroup_t *myIOgroup);

	/** 
	 * @brief Open a file in parallel
	 * 
	 * @param filename name of file
	 * @param myIOgroup IOgroup used 
	 * @param mode SCORPIO_FILE_CREATE, SCORPIO_FILE_READONLY or SCORPIO_FILE_READWRITE
	 * 
	 * @return file handle
	 */
	int scorpio_open_file(const char *filename,iogroup_t *myIOgroup, file_mode_t mode );

	/** 
	 * @brief Close file
	 * 
	 * @param fhandle file handle
	 * @param myIOgroup IOgroup used to open the file
	 * 
	 * @return error code
	 */
	int scorpio_close_file(int fhandle, iogroup_t *myIOgroup );

	/** 
	 * @brief Helper function to obtain number of dimensions of a given dataset
	 * 
	 * @param pndims Number of dimensions (pointer) returned
	 * @param fhandle file handle
	 * @param dset_name dataset name
	 * @param myIOgroup IOgroup used to open the file containing this dataset
	 * 
	 * @return error code
	 */
	int scorpio_get_dataset_ndims( int *pndims, int fhandle, char *dset_name,  iogroup_t *myIOgroup);

	/** 
	 * @brief Helper function to obtain details of dimensions of a given dataset
	 * 
	 * @param dataset_dims Detailed dimensions of the dataset (array)
	 * @param fhandle file handle
	 * @param dset_name dataset name
	 * @param myIOgroup IOgroup used to open the file containing this dataset
	 * 
	 * @return error code
	 */
	int scorpio_get_dataset_dims( int *dataset_dims, int fhandle, char *dset_name,  iogroup_t *myIOgroup);

	/** 
	 * @brief Helper function to obtain total size of a given dataset
	 * 
	 * @param mydataset_size Size of entire dataset(pointer) returned
	 * @param fhandle file handle
	 * @param dset_name dataset name
	 * @param myIOgroup IOgroup used to open the file containing this dataset
	 * 
	 * @return error code
	 */
	int scorpio_get_dataset_size( int *mydataset_size, int fhandle, char *dset_name,  iogroup_t *myIOgroup);

	/** 
	 * @brief Read a dataset from a previously opened file
	 * 
	 * @param vector pointer to which data is read
	 * @param mytype datatype of data being read
	 * @param ndims number of dimensions in target dataset
	 * @param globaldims dimensions of dataset in file
	 * @param localdims dimensions of portion of dataset that is being read
	 * @param fhandle file handle
	 * @param dset_name dataset name
	 * @param myIOgroup IOgroup used to open this file
	 * @param mypattern IO pattern used
	 * 
	 * @return error code
	 */
	int scorpio_read_dataset( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int fhandle, char *dset_name,  iogroup_t *myIOgroup, iopattern_t mypattern);

	/** 
	 * @brief Write a dataset to a previously opened file
	 * 
	 * @param vector pointer that contains data to be written
	 * @param mytype datatype of data being written
	 * @param ndims number of dimensions in target dataset
	 * @param globaldims dimensions of dataset in file
	 * @param localdims dimensions of portion of dataset that is being written
	 * @param fhandle file handle
	 * @param dset_name dataset name
	 * @param myIOgroup IOgroup used to open this file
	 * @param mypattern IO pattern used
	 * 
	 * @return error code 
	 */
	int scorpio_write_dataset( void *vector, datatype_t mytype, int ndims, int *globaldims, int *localdims, int fhandle, char *dset_name,  iogroup_t *myIOgroup, iopattern_t mypattern);

	/** 
	 * @brief Write an array of strings to a previously opened file
	 * 
	 * @param strarray array of strings to write
	 * @param nstrings number of strings to write
	 * @param fhandle file handle
	 * @param dset_name dataset name
	 * @param myIOgroup IOgroup used to open this file
	 * 
	 * @return error code 
	 */
	int scorpio_write_str_array(char **strarray, int nstrings, int fhandle, char *dset_name,  iogroup_t *myIOgroup);

	/** 
	 * @brief Read an array of strings to a previously opened file
	 * 
	 * @param pstrarray pointer to array of strings to read (storage allocated within function as needed)
	 * @param pnstrings pointer to number of strings to read
	 * @param fhandle file handle
	 * @param dset_name dataset name
	 * @param myIOgroup IOgroup used to open this file
	 * 
	 * @return error code 
	 */
	int scorpio_read_str_array(char ***pstrarray, int *pnstrings, int fhandle, char *dset_name,  iogroup_t *myIOgroup);

	/**
	 * @brief Reads a simple attribute by its name, including atomic data types like int, double or even a string.
	 *
	 * @param attr_name Name of attribute
	 * @param pattr_data OUTGOING: Pointer to attribute data (storage allocated as needed). User needs to free the space allocated
	 * @param mytype Datatype of attribute
	 * @param fhandle file handle
	 * @param primary_obj_name this can be a dataset or group or committed datatype name to which this attribute is related
	 * @param myIOgroup IOgroup used to open this file
	 *
	 * @return error code
	 */
	int scorpio_read_simple_attr(char *attr_name, void **pattr_data, datatype_t mytype, int fhandle, char *primary_obj_name,  iogroup_t *myIOgroup);

	/**
	 * @brief Reads an attribute by its name, typically used for complex attributes like arrays of atomic data types like int, double or even a string.
	 *
	 * @param attr_name Name of attribute
	 * @param pattr_data OUTGOING: Pointer to attribute data (storage allocated as needed). User needs to free the space allocated
	 * @param mytype Datatype of attribute
	 * @param pndims OUTGOING: Pointer to Number of dims (storage allocated as needed). User needs to free the space allocated
	 * @param padims OUTGOING: Pointer to array with detailed attribute dimensions (storage allocated as needed). User needs to free the space allocated
	 * @param fhandle file handle
	 * @param primary_obj_name this can be a dataset or group or committed datatype name to which this attribute is related
	 * @param myIOgroup IOgroup used to open this file
	 *
	 * @return error code
	 */
	int scorpio_read_attr(char *attr_name, void **pattr_data, datatype_t mytype, int *pndims, int **padims, int fhandle, char *primary_obj_name,  iogroup_t *myIOgroup);

	/**
	 * @brief Writes a simple attribute by its name, including atomic data types like int, double or even a string.
	 *
	 * @param attr_name Name of attribute
	 * @param attr_data Attribute data to write
	 * @param mytype Datatype of attribute
	 * @param fhandle file handle
	 * @param primary_obj_name this can be a dataset or group or committed datatype name to which this attribute is related
	 * @param myIOgroup IOgroup used to open this file
	 *
	 * @return error code
	 */
	int scorpio_write_simple_attr(char *attr_name, void *attr_data, datatype_t mytype, int fhandle, char *primary_obj_name,  iogroup_t *myIOgroup);

	/**
	 * @brief Writes an attribute by its name, typically used for complex attributes like arrays of atomic data types like int, double or even a string.
	 *
	 * @param attr_name Name of attribute
	 * @param pattr_data Pointer to attribute data (storage allocated as needed). User needs to free the space allocated
	 * @param mytype Datatype of attribute
	 * @param ndims Number of dims in attribute
	 * @param adims Array with detailed attribute dimensions
	 * @param fhandle file handle
	 * @param primary_obj_name this can be a dataset or group or committed datatype name to which this attribute is related
	 * @param myIOgroup IOgroup used to open this file
	 *
	 * @return error code
	 */	
	int scorpio_write_attr(char *attr_name, void *attr_data, datatype_t mytype, int ndims, int *adims, int fhandle, char *primary_obj_name,  iogroup_t *myIOgroup);

	/** 
	 * @brief Create a HDF5 dataset group in the given file. 
	 * 
	 * @param group_name
	 * @param fhandle file handle
	 * @param myIOgroup IOgroup used to open this file
	 * 
	 * @return groupid for the new group
	 */
	int64_t scorpio_create_dataset_group( char *group_name, int fhandle, iogroup_t *myIOgroup);

	/** 
	 * @brief Check if a HDF5 dataset group exists in the given file. 
	 * Assumes that the file is already open 
	 * @param group_name
	 * @param fhandle file handle
	 * @param myIOgroup IOgroup used to open this file
	 * 
	 * @return 0 if the group doesn't exist, 1 if the group exists
	 */
	int scorpio_group_exists(char *group_name, int fhandle, iogroup_t *myIOgroup);

	/** 
	 * @brief Check if a HDF5 dataset exists in the given file. 
	 * Assumes that the file is already open 
	 * @param dataset_name
	 * @param fhandle file handle
	 * @param myIOgroup IOgroup used to open this file
	 * 
	 * @return 0 if the dataset doesn't exist, 1 if the group exists
	 */
	int scorpio_dataset_exists(char *dataset_name, int fhandle, iogroup_t *myIOgroup );
	
	/** 
	 * @brief Closes a HDF5 dataset group previously opened
	 * 
	 * @param groupid groupid to be closed
	 * @param fhandle file handle
	 * @param myIOgroup IOgroup used to open this file
	 * 
	 * @return error code
	 */
	int scorpio_close_dataset_group( int64_t groupid, int fhandle, iogroup_t *myIOgroup);

        /**
         * @brief Creates a soft link to an object in the given file
         *
         * @param target path to the target object
         * @param link_name name of the new soft link
         * @param fhandle file handle
         * @param myIOgroup IOgroup used to open the file
         *
         * @return error code
         */
        int scorpio_create_link(char *target, int64_t link_loc_id, char *link_name, int fhandle, iogroup_t *myIOgroup);

#ifdef __cplusplus
}
#endif

#endif 

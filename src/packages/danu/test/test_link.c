/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
#include <stdlib.h>
#include <stdio.h>

#include <hdf5.h>

#include <danu_error.h>
#include <danu_h5_error.h>
#include <danu_types.h>
#include <danu_utils.h>
#include <danu_file.h>
#include <danu_group.h>
#include <danu_memory.h>
#include <danu_dataset.h>

#include <danu_link.h>

#define TEST_FILE "deleteme.h5"
#define LINK_FILE "link.h5"

#define GRPA_NAME  "GroupA"
#define GRPB_NAME  "GroupB"

#define FOO_NAME  "foo"
#define BAR_NAME  "bar"

#define DATA_NAME "Data"
#define DATA_DIM  2
#define DATA_SIZE 10

int main(int argc, char **argv)
{

   hid_t fid, lid;
   hid_t grp_a, grp_b;
   hid_t data,link;
   hid_t foo_id, bar_id;
   int *stuff;
   hsize_t *size,i,num;
   char obj_name[128],link_name[128];

   /* Genreate data */
   size = DANU_MALLOC(hsize_t,DATA_DIM);
   for(i=0;i<DATA_DIM;i++)
       size[i] = DATA_SIZE;
   num = danu_compute_sizeh(DATA_DIM,size);    
   printf("Generate data for %lu elements\n",(long unsigned) num);
   stuff = DANU_MALLOC(int,num);
   danu_rand_data_int(0,DATA_SIZE,num,stuff);

   fid = danu_file_create(TEST_FILE);
   lid = danu_file_create(LINK_FILE);
   grp_a = danu_group_create(fid,GRPA_NAME);

   /* Create a link foo->bar */

#if 0
   /* SOFT LINK EXAMPLE */
   sprintf(obj_name,"/%s/%s",GRPB_NAME,BAR_NAME);
   danu_link_create_soft(grp_a,FOO_NAME,obj_name);

   /* Create the GroupB/bar */
   grp_b = danu_group_create(fid,GRPB_NAME,FALSE);
   bar_id = danu_group_create(grp_b,BAR_NAME,FALSE);
#endif

#if 0
   /* Create the GroupB/bar */
   grp_b = danu_group_create(fid,GRPB_NAME,FALSE);
   bar_id = danu_group_create(grp_b,BAR_NAME,FALSE);

   /* HARD LINK EXAMPLE */
   sprintf(obj_name,"/%s/%s", GRPB_NAME,BAR_NAME);
   danu_link_create_hard(grp_a,FOO_NAME,grp_b,BAR_NAME);
#endif

   /* EXTERNAL LINK EXAMPLE */
   
   /* Create the GroupB/bar */
   grp_b = danu_group_create(lid,GRPB_NAME);
   bar_id = danu_group_create(grp_b,BAR_NAME);
   sprintf(obj_name,"/%s/%s",GRPB_NAME,BAR_NAME);
   H5Lcreate_external(LINK_FILE,obj_name,grp_a,FOO_NAME,H5P_DEFAULT,H5P_DEFAULT);


   /* Write using the link name */
   sprintf(link_name,"/%s/%s",GRPA_NAME,FOO_NAME);
   printf("Opening link %s\n", link_name);
   link = danu_group_open(fid,link_name);
   if( H5_ISA_INVALID_ID(link)) {
       DANU_ERROR_MESS("Failed to open link");
       goto FAIL_EXIT;
   }
   else {
       data = danu_dataset_create(link,DATA_NAME,H5T_STD_I32LE,DATA_DIM,size,FALSE);
       danu_dataset_write(data,NULL,H5T_NATIVE_INT,DATA_DIM,size,stuff);
       danu_dataset_close(data);
   }


   //danu_group_close(grp_a);
   danu_file_close(fid);
   danu_file_close(lid);
   

   return DANU_SUCCESS;

FAIL_EXIT:

   return DANU_FAIL;

}


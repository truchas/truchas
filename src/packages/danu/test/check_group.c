/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
/* ****************************************************************************
*                                                                             *
*  Danu Unit Test                                                             *
*  Unit: Group                                                                  *
*                                                                             *
* Requires Check Unit Test software package                                   *
* http://check.sourceforge.net/                                               *
*                                                                             *
* *************************************************************************** */ 
#include <stdlib.h>
#include <check.h>

#include <hdf5.h>

#include <danu_types.h>
#include <danu_h5_error.h>
#include <danu_memory.h>
#include <danu_error.h>
#include <danu_utils.h>
#include <danu_file.h>

#include <danu_group.h>

#define FILENAME "deleteme-group.h5"


/* BEGIN TESTS */

START_TEST (test_group_list)
{
    const char test_file[] = FILENAME;
    char group_name[16];
    hid_t  fid,gid,gidA,gidB;
    herr_t status;
    int num_found, i;
    char **names;
    int num = 3;

    fid = danu_file_create(test_file);

    names = DANU_MALLOC(char *,num);
    
    gidA = danu_group_create(fid, "GroupA");

    status = danu_group_get_subgroups(gidA,num,names,&num_found);
    fail_unless(H5_RETURN_OK(status),
		"Failed to search for subgroups");
    fail_unless(num_found == 0,
		"Failed to return the correct number of subgroups (0)");

    for(i=0;i<num;i++) {
      sprintf(group_name,"Group%04d",i);
      gid = danu_group_create(gidA,group_name);
      danu_group_close(gid);
    }

    status = danu_group_get_subgroups(gidA,num,names,&num_found);
    fail_unless(H5_RETURN_OK(status),
		"Failed to search for subgroups");
    fail_unless(num_found == num,
		"Failed to return the correct number of subgroups");

    danu_group_close(gidA);
    danu_file_close(fid);

}
END_TEST

Suite *
group_suite (void)
{
    Suite *s = suite_create("Danu Group");

    TCase *group_utils = tcase_create("Group Utils");
    tcase_add_test(group_utils,test_group_list);
    suite_add_tcase(s,group_utils);

    return s;
}

int main(void)
{
  int number_failed;
  Suite *s = group_suite ();
  SRunner *sr = srunner_create (s);

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}




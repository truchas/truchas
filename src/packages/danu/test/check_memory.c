/* ****************************************************************************
*                                                                             *
*  Danu Unit Test                                                             *
*  Unit: File                                                                  *
*                                                                             *
* Requires Check Unit Test software package                                   *
* http://check.sourceforge.net/                                               *
*                                                                             *
* *************************************************************************** */ 
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include <check.h>

#include <danu_memory.h>

/* BEGIN TESTS */

START_TEST (test_bad_bytes)
{

  char *p;
  size_t bytes;

  bytes = -1;

  p = (char * ) danu_malloc(bytes);

  fail_unless(p == NULL,
	      "Failed to return bad pointer with bad input");

  free(p);
    
}
END_TEST

START_TEST(test_bad_free)
{
  char *p = NULL;

  /* This won't fail but should see a warning message */
  danu_free(p);

}
END_TEST

START_TEST(test_calloc)
{
  int *p;
  int i;
  int num = 10 ;
  size_t bytes = sizeof(int)*num;

  p = (int *) danu_calloc(bytes);
  fail_unless(p != NULL,
	      "Failed to calloc an array"); 
  for(i=0; i< num; i++) {
    fail_unless(p[i] == 0,
		"Failed to populate array with zeros");
  }
  danu_free(p);

  /* Now test the macro */
  p = DANU_CALLOC(int,10);
  fail_unless(p != NULL,
	      "Failed to calloc an array"); 
  for(i=0; i< num; i++) {
    fail_unless(p[i] == 0,
		"Failed to populate array with zeros");
  }

  DANU_FREE(p);

   
}
END_TEST

START_TEST(test_realloc)
{
  int *p, *n;
  int i;
  int num = 100 ;
  int new_num = 110;
  size_t bytes = sizeof(int)*num;

  p = (int *) danu_malloc(bytes);
  fail_unless(p != NULL,
	      "Failed to malloc an array"); 

  n = danu_realloc(p,new_num);
  fail_unless(n != NULL,
	      "Failed to realloc");

  /* This will SEGV if something bad happens */
  for(i=0; i<new_num; i++) {
    printf("i=%d n[%d]=%d\n", i, i, n[i]);
  }

  danu_free(n);
}
END_TEST

START_TEST(test_malloc)
{
  float *p;
  int i;
  int num = 100 ;

  p = DANU_MALLOC(float,num);
  fail_unless(p != NULL,
	      "Failed to malloc an array"); 

  for(i=0; i<num; i++) {
    printf("i=%d n[%d]=%f\n", i, i, p[i]);
  }

  DANU_FREE(p); 
}
END_TEST



Suite *
memory_suite (void)
{
    Suite *s = suite_create("Danu Memory Tools");
	
    TCase *bad_tests = tcase_create("Tests Designed to Fail");
    tcase_add_test(bad_tests,test_bad_bytes);
    tcase_add_test(bad_tests,test_bad_free);
    suite_add_tcase(s,bad_tests);
    
    TCase *mem_tests = tcase_create("Memory Tests");
    tcase_add_test(mem_tests,test_calloc);
    tcase_add_test(mem_tests,test_realloc);
    tcase_add_test(mem_tests,test_malloc);
    suite_add_tcase(s,mem_tests);
    return s;
}

int main(void)
{
  int number_failed;
  Suite *s = memory_suite ();
  SRunner *sr = srunner_create (s);

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}




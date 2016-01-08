/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
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

#include "random_generators.h"

#define MAX_LOOP 10
#define TEST_SEED 10647859

/* Test suite utils for this series of tests */
int define_test_seed(void);


int define_test_seed()
{
	return -1 * TEST_SEED;
}

/* BEGIN TESTS */

START_TEST (test_random_int)
{
	int seed;
	int diff_seed;
	time_t now;
	int *num, *compare;
	int i;
	
	num = (int *) malloc(sizeof(int)*MAX_LOOP);
	compare = (int *) malloc(sizeof(int)*MAX_LOOP);
	
	seed = define_test_seed();
	for(i = 0; i < MAX_LOOP; i++ ) {
		num[i] = generate_random_int(&seed);
	}
	
	/* Reset the seed to get the same string of numbers */
	seed = define_test_seed();
	for(i=0; i < MAX_LOOP; i++ ) {
		compare[i] = generate_random_int(&seed);
	}
	
	/* Should be identical */
	for(i=0;i<MAX_LOOP;i++) {
		//printf("i=%d num=%d compare=%d\n", i, num[i], compare[i]);
		fail_unless( num[i] == compare[i],
					"Did not generate identical random numbers in int case");
	}
	
	
	/* Now generate a different list of numbers */
	memset(compare,0,sizeof(int)*MAX_LOOP);
	
	time(&now);
	diff_seed = (int) now * -1;
	for(i=0; i < MAX_LOOP; i++) {
		compare[i] = generate_random_int(&diff_seed);
	}
	
	/* Should be different! */
	for(i = 0; i < MAX_LOOP; i ++) {
		//printf("i=%d num=%d compare=%d\n", i, num[i], compare[i]);
		fail_unless( num[i] != compare[i],
					"Failed to generate different random numbers in int case");
	}
	
	free(num);
	free(compare);
    
}
END_TEST

START_TEST(test_random_int_data)
{
	int iseed = define_test_seed();
	int diff_seed = (int) time(NULL);
	int * data, *compare;
	int size = generate_random_bound_int(64,2048,&diff_seed);
	size_t nbytes = size*sizeof(int);
	
	data = (int * ) malloc(nbytes);
	compare = (int *) malloc(nbytes);
	
	iseed = define_test_seed();
	generate_random_int_data(data,size,&iseed);
	iseed = define_test_seed();
	generate_random_int_data(compare,size,&iseed);
	
	fail_unless( memcmp(data,compare,nbytes) == 0,
				"Failed to generate identical int data");
	
	/* Now generate other data */
	memset(compare,0,nbytes);
	generate_random_int_data(compare,size,&diff_seed);
	
	fail_unless( memcmp(data,compare,nbytes) != 0,
				"Failed to generate different int data");
	
	free(data);
	free(compare);
	 
	
}
END_TEST

START_TEST ( test_random_bound_int)
{
	int tmp, min, max;
	int * num, *compare;
	int i;
	int iseed;
	time_t now;
	
	/* Determine a min max bound randonly */
	time(&now);
	iseed = (int) now;
	min = generate_random_int(&iseed);
	max = generate_random_int(&iseed);
	
	fail_unless( min != max,
				"Successive calls to generator produced identical results");
	
	num = (int *) malloc(sizeof(int)*MAX_LOOP);
	iseed = define_test_seed();
	for(i=0; i < MAX_LOOP; i++ ) {
		num[i] = generate_random_bound_int(min,max,&iseed);
	}
	
	/* Now check each number */
	if ( min > max ) {
		tmp = min;
		min = max;
		max = tmp;
	}
	for(i=0; i < MAX_LOOP; i++ ) {
		//printf("i=%d min=%d num=%d max=%d\n", i, min, num[i], max);
		fail_unless( ( num[i] >= min && num[i] <= max ),
					"Returned random number out of range");
	}
	
	free(num);
	
}
END_TEST

START_TEST(test_random_bound_int_data)
{
	int iseed = define_test_seed();
	int diff_seed = (int) time(NULL);
	int * data, *compare;
	int size = generate_random_bound_int(64,2048,&diff_seed);
	size_t nbytes = size*sizeof(int);
	int min, max;
	
	min = generate_random_bound_int(64,4096,&diff_seed);
	max = generate_random_bound_int(4097, 10000, &diff_seed);
	
	data = (int * ) malloc(nbytes);
	compare = (int *) malloc(nbytes);
	
	iseed = define_test_seed();
	generate_random_bound_int_data(min,max,data,size,&iseed);
	iseed = define_test_seed();
	generate_random_bound_int_data(min,max,compare,size,&iseed);
	
	fail_unless( memcmp(data,compare,nbytes) == 0,
				"Failed to generate identical int data");
	
	/* Now generate other data */
	memset(compare,0,nbytes);
	generate_random_bound_int_data(min,max,compare,size,&diff_seed);
	
	fail_unless( memcmp(data,compare,nbytes) != 0,
				"Failed to generate different int data");
	
	
	free(data);
	free(compare);
	
	
}
END_TEST

START_TEST (test_random_float)
{
	int seed;
	int diff_seed;
	time_t now;
	float *num, *compare;
	int i;
	
	num = (float *) malloc(sizeof(float)*MAX_LOOP);
	compare = (float *) malloc(sizeof(float)*MAX_LOOP);
	
	seed = define_test_seed();
	for(i = 0; i < MAX_LOOP; i++ ) {
		num[i] = generate_random_float(&seed);
	}
	
	/* Reset the seed to get the same string of numbers */
	seed = define_test_seed();
	for(i=0; i < MAX_LOOP; i++ ) {
		compare[i] = generate_random_float(&seed);
	}
	
	/* Should be identical */
	for(i=0;i<MAX_LOOP;i++) {
		//printf("i=%d num=%d compare=%d\n", i, num[i], compare[i]);
		fail_unless( num[i] == compare[i],
					"Did not generate identical random numbers in float case");
	}
	
	
	/* Now generate a different list of numbers */
	memset(compare,0,sizeof(float)*MAX_LOOP);
	
	time(&now);
	diff_seed = (int) now * -1;
	for(i=0; i < MAX_LOOP; i++) {
		compare[i] = generate_random_float(&diff_seed);
	}
	
	/* Should be different! */
	for(i = 0; i < MAX_LOOP; i ++) {
		//printf("i=%d num=%d compare=%d\n", i, num[i], compare[i]);
		fail_unless( num[i] != compare[i],
					"Failed to generate different random numbers in float case");
	}
	
	free(num);
	free(compare);
    
}
END_TEST

START_TEST(test_random_float_data)
{
	int iseed = define_test_seed();
	int diff_seed = (int) time(NULL);
	float * data, *compare;
	int size = generate_random_bound_int(64,2048,&diff_seed);
	size_t nbytes = size*sizeof(float);
	
	data = (float * ) malloc(nbytes);
	compare = (float *) malloc(nbytes);
	
	iseed = define_test_seed();
	generate_random_float_data(data,size,&iseed);
	iseed = define_test_seed();
	generate_random_float_data(compare,size,&iseed);
	
	fail_unless( memcmp(data,compare,nbytes) == 0,
				"Failed to generate identical float data");
	
	/* Now generate other data */
	memset(compare,0,nbytes);
	generate_random_float_data(compare,size,&diff_seed);
	
	fail_unless( memcmp(data,compare,nbytes) != 0,
				"Failed to generate different float data");
	
	free(data);
	free(compare);
	
	
}
END_TEST

START_TEST (test_ran3)
{
	long seed;
	long diff_seed;
	time_t now;
	float * num, *compare;
	int i;
	
	
	num = (float *) malloc(sizeof(float)*MAX_LOOP);
	compare = (float *) malloc(sizeof(float)*MAX_LOOP);
	
	seed = -1 * TEST_SEED;
	for(i=0; i < MAX_LOOP; i++ ) {
		num[i] = ran3(&seed);
	}
	
	/* ran3 changes the value of seed */
	seed = -1 * TEST_SEED;
	for (i=0; i < MAX_LOOP; i++ ) {
		compare[i] = ran3(&seed);
	}
	
	for(i=0; i < MAX_LOOP; i++ ) {
		fail_unless( num[i] == compare[i],
					"Did not generate identical sequence");
	}
		
    /* Generate with a different seed */
	time(&now);
    diff_seed = (long) now;	
	for(i=0; i< MAX_LOOP; i++) {
		compare[i] = ran3(&diff_seed);
	}
	
	for(i=0; i < MAX_LOOP; i++ ) {
		fail_unless( num[i] != compare[i],
					"Did not generate a different sequence");
	}
	
	free(num);
	free(compare);
	
}
END_TEST

START_TEST (test_random_double)
{
	int seed;
	int diff_seed;
	time_t now;
	double *num, *compare;
	int i;
	
	num = (double *) malloc(sizeof(double)*MAX_LOOP);
	compare = (double *) malloc(sizeof(double)*MAX_LOOP);
	
	seed = define_test_seed();
	for(i = 0; i < MAX_LOOP; i++ ) {
		num[i] = generate_random_double(&seed);
	}
	
	/* Reset the seed to get the same string of numbers */
	seed = define_test_seed();
	for(i=0; i < MAX_LOOP; i++ ) {
		compare[i] = generate_random_double(&seed);
	}
	
	/* Should be identical */
	for(i=0;i<MAX_LOOP;i++) {
		//printf("i=%d num=%1.9e compare=%1.9e\n", i, num[i], compare[i]);
		fail_unless( num[i] == compare[i],
					"Did not generate identical random numbers in double case");
	}
	
	
	/* Now generate a different list of numbers */
	memset(compare,0,sizeof(float)*MAX_LOOP);
	
	time(&now);
	diff_seed = (int) now * -1;
	for(i=0; i < MAX_LOOP; i++) {
		compare[i] = generate_random_double(&diff_seed);
	}
	
	/* Should be different! */
	for(i = 0; i < MAX_LOOP; i ++) {
		//printf("i=%d num=%d compare=%d\n", i, num[i], compare[i]);
		fail_unless( num[i] != compare[i],
					"Failed to generate different random numbers in double case");
	}
	
	free(num);
	free(compare);
	
	
    
}
END_TEST

START_TEST(test_random_double_data)
{
	int iseed = define_test_seed();
	int diff_seed = (int) time(NULL);
	double * data, *compare;
	int size = generate_random_bound_int(64,2048,&diff_seed);
	size_t nbytes = size*sizeof(double);
	
	data = (double * ) malloc(nbytes);
	compare = (double *) malloc(nbytes);
	
	iseed = define_test_seed();
	generate_random_double_data(data,size,&iseed);
	iseed = define_test_seed();
	generate_random_double_data(compare,size,&iseed);
	
	fail_unless( memcmp(data,compare,nbytes) == 0,
				"Failed to generate identical float data");
	
	/* Now generate other data */
	memset(compare,0,nbytes);
	generate_random_double_data(compare,size,&diff_seed);
	
	fail_unless( memcmp(data,compare,nbytes) != 0,
				"Failed to generate different float data");
	
	free(data);
	free(compare);
	
	
}
END_TEST

START_TEST(test_random_string)
{
	int iseed = define_test_seed();
	int diff_seed;
	char * string, *compare;
	int size;
	int i;
	
	size = generate_random_bound_int(64,2048,&iseed);
	string = (char *) malloc(size);
	compare = (char *) malloc(size);
	
	iseed = define_test_seed();
	generate_random_string(string,size,&iseed);
	iseed = define_test_seed();
	generate_random_string(compare,size,&iseed);
	
	fail_unless(strcmp(string,compare) == 0,
			    "Failed to generate identical string");
	
	/* Now generate a differnt string */
	memset(compare,0,size);
	
	diff_seed = (int) time(NULL) * -1;
	generate_random_string(compare,size,&diff_seed);
	
	fail_unless( strcmp(string,compare) != 0,
				"Failed to generate a different sequence string");
	
	free(string);
	free(compare);

	
}
END_TEST
START_TEST(test_random_fortran_string)
{
	int iseed = define_test_seed();
	int diff_seed;
	char * string, *compare;
	int size, trim_len;
	int i;
	
	size = generate_random_bound_int(64,2048,&iseed);
	string = (char *) malloc(size);
	compare = (char *) malloc(size);
	
	/* Pad strings with white space */
	trim_len = generate_random_bound_int(size/3,size/2,&iseed);
	
	memset(string,' ',size);
	memset(compare,' ',size);
	iseed = define_test_seed();
	generate_random_string(string,trim_len,&iseed);
	iseed = define_test_seed();
	generate_random_string(compare,trim_len,&iseed);
	
	fail_unless(strncmp(string,compare,size) == 0,
			    "Failed to generate identical fortran strings");
	
	/* Now generate a differnt string */
	memset(compare,' ',size);
	
	diff_seed = (int) time(NULL) * -1;
	generate_random_string(compare,trim_len,&diff_seed);
	
	fail_unless( strncmp(string,compare,size) != 0,
				"Failed to generate a different sequence fortran string");
	
	free(string);
	free(compare);
	
	
}
END_TEST

START_TEST(test_pi_compute)
{
	
	static double PI_D = 3.1415926535897932384626433832795028841971693993751058209749445923078;
	static float PI_F = 3.1415926535897932384626433832795028841971693993751058209749445923078;
	int num_pts = 1000000;
	int i;
	double xd,yd, rd, ed;
	float xf,yf,rf, ef;
	long lseed;
	int iseed;
	int dcnt, fcnt;
	float pi_f_est;
	double pi_d_est;
	
	time_t now;
	
	time(&now);
	iseed = (int) now;
	lseed = (long) now;
	
	dcnt =0;
	fcnt = 0;
	for(i=1; i<=num_pts; i++) {
		
      xd = generate_random_double(&iseed);
	  yd = generate_random_double(&iseed);
	  rd = (xd*xd) + (yd*yd);
	  if ( rd <= 1.0 ) {
		  dcnt++;
	  }
			
		
	  xf = ran3(&lseed);
	  yf = ran3(&lseed);
	  rf = (xf*xf) + (yf*yf);
	  if ( rf <= 1.0 ) {
		fcnt++;
	  }
		
	}
	
	/* Now compare to PI values */
	pi_d_est = (double)dcnt/num_pts * 4.0;
	pi_f_est = (float)fcnt/num_pts * 4.0;
	
	ed = fabs(pi_d_est - PI_D)/PI_D * 100.0;
	ef = fabsf(pi_f_est - PI_F)/PI_F * 100.0;
	
	printf("PI Compute Results\n");
	printf("My Double Generator pi_est=%1.9e error=%1.4e%%\n", pi_d_est, ed);
	printf("ran3 Generator pi_est=%1.9e error=%1.4e%%\n", pi_f_est, ef);
}
END_TEST	

Suite *
random_suite (void)
{
    Suite *s = suite_create("Random Generator");
	
    TCase *int_utils = tcase_create("Random Int Utils");
    tcase_add_test(int_utils,test_random_int);
	tcase_add_test(int_utils,test_random_bound_int);
	tcase_add_test(int_utils,test_random_int_data);
	tcase_add_test(int_utils,test_random_bound_int_data);
    suite_add_tcase(s,int_utils);
    
    TCase *float_utils = tcase_create("Random Float Utils");
    tcase_add_test(float_utils,test_random_float);
	tcase_add_test(float_utils,test_ran3);
	tcase_add_test(float_utils,test_random_float_data);
    suite_add_tcase(s,float_utils);

	TCase *double_utils = tcase_create("Random Double Utils");
    tcase_add_test(double_utils,test_random_double);
	tcase_add_test(double_utils,test_random_double_data);
    suite_add_tcase(s,double_utils);
	
	TCase *char_utils = tcase_create("Random Char Utils");
    tcase_add_test(char_utils,test_random_string);
	tcase_add_test(char_utils,test_random_fortran_string);
    suite_add_tcase(s,char_utils);
	
	TCase *pi_compute = tcase_create("PI Compute Test");
	tcase_add_test(pi_compute,test_pi_compute);
	suite_add_tcase(s,pi_compute);
	
    return s;
}

int main(void)
{
  int number_failed;
  Suite *s = random_suite ();
  SRunner *sr = srunner_create (s);

  srunner_run_all (sr, CK_NORMAL);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}




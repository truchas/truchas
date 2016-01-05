/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

#ifndef RAND_GEN_H
#define RAND_GEN_H


/* From Numerical Recipes */
float ran3(long *idum);

/* Genrators */
int    generate_random_int(int *seed);
void   generate_random_int_data(int * data, int size, int *iseed);
void   generate_random_bound_int_data(int min, int max, int * data, int size, int *iseed);
int    generate_random_bound_int(int min,int, int *seed);
float  generate_random_float(int *seed);
void   generate_random_float_data(float * data, int size, int *iseed);
double generate_random_double(int *seed);
void   generate_random_double_data(double * data, int size, int *iseed);
void   generate_random_string(char *, int, int *seed);
void   generate_random_fortran_string(char *, int, int *iseed);

#endif

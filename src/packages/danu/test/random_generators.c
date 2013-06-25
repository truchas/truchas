/* *************************************************************************** *
*                                                                              *
*                                                                              *
*                             Copyright  (C) 20xx,                             *
*                      Los Alamos National Security, LLC                       *
*                                                                              *
*                             LA-CC-xxxxxx                                     *
*                                                                              *
* **************************************************************************** */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "random_generators.h"


#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)


/* From numerical recipes */
float ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=labs(MSEED-labs(*idum));
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

/* Generators */
int generate_random_int(int *iseed)
{
  int tmp;
  static int flag = 0;

  if ( *iseed < 0 || flag == 0 ) {
	  tmp = abs(*iseed);
	  srand(tmp);
	  *iseed = tmp;
	  flag = 1;
  }

  return rand();
}

void generate_random_int_data(int *data, int size, int *iseed)
{
	int i;
	
	for(i=0;i<size;i++) {
		data[i] = generate_random_int(iseed);
	}
	
}
		
void generate_random_bound_int_data(int min, int max, int *data, int size, int *iseed)
{
	int i;
	
	for(i=0;i<size;i++) {
		data[i] = generate_random_bound_int(min,max,iseed);
	}
	
}

int generate_random_bound_int(int min, int max, int *iseed)
{
  int r = generate_random_int(iseed);
  int shift = max > min ? min : max;
  float scale = max > min ? (max-min+1) : (min-max+1);
  int result = shift + (int) ( scale*r/(RAND_MAX+1.0) );
	
  return result;

}


float generate_random_float(int *iseed)
{
   int r = generate_random_int(iseed);
   float d = (float)RAND_MAX + 1.0;
   float result;
   result = (float)r/d;
   
   return result;
   
}

void generate_random_float_data(float *data, int size, int *iseed)
{
	int i;
	
	for(i=0;i<size;i++) {
		data[i] = generate_random_float(iseed);
	}
	
}


double generate_random_double(int *iseed)
{

  int r = generate_random_int(iseed);
  double d = (double)RAND_MAX + 1.0;
  double result = (double)r/d;

  //printf("result=%1.8e\n",result);

  return result;

}

void generate_random_double_data(double *data, int size, int *iseed)
{
	int i;
	
	for(i=0;i<size;i++) {
		data[i] = generate_random_double(iseed);
	}
	
}



void generate_random_string(char * string, int size, int *iseed)
{


  static const char *data =
     "0123456789"
     "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
     "abcdefghijklmnopqrstuvwxyz"
     "!@#$%^&*()_-+={}[]|,.<>?/:;";

  static int data_size = 10 + 2*26 + 27;
  int i,idx;

  for(i=0; i<size; i++) {
    idx = generate_random_bound_int(0,data_size-1,iseed);
    string[i] = data[idx];
  }

  string[size-1] = '\0';
}



void generate_random_fortran_string(char *f_string, int size, int *iseed)
{
  char tmp[2];

  /* These routines NULL terminate the string */
  generate_random_string(f_string,size,iseed);
  generate_random_string(tmp,2,iseed);

  f_string[size-1] = tmp[0];

}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

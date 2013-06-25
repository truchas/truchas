
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "vf.h"

double VF_QMC (Poly *poly_i, Poly *poly_j, Sampling *sampling);

main()
{
  Poly   poly_i, poly_j;
  int    Comm=0;
  int	 i, n, nsegi, nsegj, ncnt=10000;
  int	 type[] = { VF_RANDOM_SAMPLE, VF_UNIFORM_SAMPLE,  
		    VF_JITTER_SAMPLE, VF_HALTON_SAMPLE }; 
  double start, stop, vfa, vf, d=-1.0, r_over_d;
  Sampling sampling;

  VF_Initialize(1, 2, 1, 1, Comm);
  VF_DefineEnclosure("testFF", 1, 0, 0.0, 0, NULL, 0);
  
  poly_i.np	= 3;
  
  poly_i.p[0].x = 0.0;
  poly_i.p[0].y = 0.0;
  poly_i.p[0].z = 0.0;
  
  poly_i.p[1].x = 0.0;
  poly_i.p[1].y = 0.0;
  poly_i.p[1].z = 3.0;
  
  poly_i.p[2].x = 3.0;
  poly_i.p[2].y = 0.0;
  poly_i.p[2].z = 0.0;
  
  poly_j.np = 3;
  
  poly_j.p[0].x = 0.0;
  poly_j.p[0].y = 4.0;
  poly_j.p[0].z = 0.0;
  
  poly_j.p[1].x = 0.0;
  poly_j.p[1].y = 4.0;
  poly_j.p[1].z =-2.0;
  
  poly_j.p[2].x = 2.0;
  poly_j.p[2].y = 4.0;
  poly_j.p[2].z = 0.0;
  
  VF_PolyNormal(&poly_i,&(poly_i.normal));
  VF_PolyNormal(&poly_j,&(poly_j.normal));
  
  r_over_d = d/(2.0*sqrt(1.0/M_PI));
  
  fprintf(stderr,"1 x 1 parallel plates, seperation = %f, r/d = %f\n",
	  d,r_over_d);
  
  start = VF_Clock();
  for (vfa=0.0, n=0; n<1; n++) {
      vfa += VF_CalcVF_Analytic(&poly_i, &poly_j);
  }
  stop = VF_Clock();
  //vfa  /= (double)ncnt;
  fprintf(stderr,"VF_Analytic() = %8.4f  (vf = %.10f)\n\n",stop-start,vfa);

  for (nsegj=0; nsegj<4; nsegj++) {
    sampling.method = type[nsegj];
    for (nsegi=30; nsegi<=30; nsegi++) {
      sampling.n = nsegi*nsegi;
      sampling.index = 11;
      start = VF_Clock();
      for (vf=0.0, n=0; n<ncnt; n++) {
        vf += VF_QMC(&poly_i, &poly_j, &sampling);
      }
      stop = VF_Clock();
      vf  /= (double)ncnt;
      fprintf(stderr,"VF_QMC()      = %8.4f  (vf = %.10f, err = %.4f)  (nsamples = %d)\n",
    	      stop-start,vf,100.0*(vfa-vf)/vfa,nsegi*nsegi);
    }
  }

  for (nsegi=1; nsegi<=5; nsegi++) {
    start = VF_Clock();
    for (vf=0.0, n=0; n<ncnt; n++) {
      vf += VF_CalcVF_Contour(&poly_i, &poly_j, nsegi, nsegi);
    }
    stop = VF_Clock();
    vf  /= (double)ncnt;
    fprintf(stderr,"\n");
    fprintf(stderr,"VF_Contour()  = %8.4f  (vf = %.10f, err = %.4f)  (nseg = %d)\n",
            stop-start,vf,100.0*fabs(vfa-vf)/vfa,nsegi);
    start = VF_Clock();
    for (vf=0.0, n=0; n<ncnt; n++) {
      vf += VF_CalcVF_Gauss(&poly_i, &poly_j, nsegi, nsegi);
    }
    stop = VF_Clock();
    vf  /= (double)ncnt;
    fprintf(stderr,"VF_Gauss()    = %8.4f  (vf = %.10f, err = %.4f)  (nseg = %d)\n",
            stop-start,vf,100.0*fabs(vfa-vf)/vfa,nsegi);
  }
  /*
  for (nsegi=1; nsegi<=3; nsegi++) {
    for (nsegj=1; nsegj<=3; nsegj++) {
      start = VF_Clock();
      for (vf=0.0, n=0; n<ncnt; n++) {
    	vf += VF_CalcVF_Contour(&poly_i, &poly_j, nsegi, nsegj);
      }
      stop = VF_Clock();
      vf  /= (double)ncnt;
      fprintf(stderr,"\n");
      fprintf(stderr,"VF_Contour()  = %8.4f  (vf = %.10f)  (nseg = %d, %d)\n",
    	      stop-start,vf,nsegi,nsegj);
      start = VF_Clock();
      for (vf=0.0, n=0; n<ncnt; n++) {
    	vf += VF_CalcVF_Gauss(&poly_i, &poly_j, nsegi, nsegj);
      }
      stop = VF_Clock();
      vf  /= (double)ncnt;
      fprintf(stderr,"VF_Gauss()    = %8.4f  (vf = %.10f)  (nseg = %d, %d)\n",
    	      stop-start,vf,nsegi,nsegj);
    }
  }
  */
}  

double VF_QMC (Poly *poly_i, Poly *poly_j, Sampling *sampling)
{
  int	  n;
  double  rsq, cos_theta_i, cos_theta_j;
  double  area, da, delta_ff, vf=0.0, tmp;
  Ray	  ray;
  //Point3  samples_i[100000], samples_j[100000];
  Point  samples_i[100000], samples_j[100000];
  //Vector3 normal_i, normal_j, r;
  Vector normal_i, normal_j, r;
  static  int first=1, NS=0;
  //static  Point3 uv_samples[200000];
  static  Point uv_samples[200000];

  if (first || sampling->n!=NS) {
    VF_SetupSampling(sampling, uv_samples);
    first = 0;
    NS    = sampling->n;
  }
  VF_SamplesUVtoXYZ(poly_i, sampling, uv_samples, samples_i);
  VF_SamplesUVtoXYZ(poly_j, sampling, uv_samples, samples_j);
  VF_PolyNormal(poly_i,&normal_i);
  VF_PolyNormal(poly_j,&normal_j);
  area = VF_PolyArea(poly_j);
  da   = area/(double)(sampling->n);
  for (n=0; n<sampling->n; n++) {
    V3_Sub(&samples_j[n],&samples_i[n],&r);
    ray.O = samples_i[n];
    ray.D = r;
    V3_Normalize(&ray.D,tmp);
    rsq 	= V3_Dot(&r,&r);
    cos_theta_i = V3_Dot(&ray.D,&normal_i);
    cos_theta_j = V3_Dot(&ray.D,&normal_j);
    delta_ff	= -(cos_theta_i*cos_theta_j)/(M_PI*rsq+da);
    if (delta_ff>0.0) vf += delta_ff;
  }
  return (da*vf);
}

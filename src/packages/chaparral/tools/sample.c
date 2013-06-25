
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <values.h>

#include "vf_types.h"

main(int argc, char *argv[])
{
    FILE   *fp;
    char   filename[80];
    char   *name[] = { "random", 
                       "uniform", 
                       "jitter", 
                       "halton" };
    int    type[] = { RANDOM_SAMPLE, 
                      UNIFORM_SAMPLE, 
                      JITTER_SAMPLE, 
                      HALTON_SAMPLE };
    int    i, j, n, ns0=100, ns1, option, errflag=0;
    Sampling sampling;
    Point3 samples_tri[100000];
    Point3 samples_quad[100000];
    Point3 p0={0.0,0.0,0.0};
    Point3 p1={1.0,0.0,0.0};
    Point3 p2={1.0,1.0,0.0};
    Point3 p3={0.0,1.0,0.0};
    Point3 p4={0.5,0.866,0.0};
    Poly   poly_tri, poly_quad;
    
    while ((option=getopt(argc, argv, "n:")) != EOF) {
        switch (option) {
        case 'n':
            ns0 = atoi(optarg);
            break;
        case '?':
            errflag++;
            break;
        }
    }
    if (errflag) {
        fprintf(stderr,
                "usage: sample [-n #_of_samples]\n");
        exit (2);
    }
    
    poly_tri.np    = 3;
    poly_tri.p[0]  = p0;
    poly_tri.p[1]  = p1;
    poly_tri.p[2]  = p4;
    poly_quad.np   = 4;
    poly_quad.p[0] = p0;
    poly_quad.p[1] = p1;
    poly_quad.p[2] = p2;
    poly_quad.p[3] = p3;
    
    InitRandomSeed();
    for (i=0; i<5; i++) {
        if (type[i]==UNIFORM_SAMPLE || type[i]==JITTER_SAMPLE) {
            n   = (int)(sqrt((double)ns0)+0.5);
            ns1 = n*n;
        } else {
            ns1 = ns0;
        }
        sampling.method = type[i];
        sampling.n      = ns1;
        sampling.index  = 3;
        SamplePoly(&poly_tri, &sampling, samples_tri);
        SamplePoly(&poly_quad, &sampling, samples_quad);
        sprintf(filename,"sample_tri_%s.dat",name[i]);
        fp = fopen(filename,"w");
        for (j=0; j<ns1; j++) {
            fprintf(fp,"%g %g\n",samples_tri[j].x,samples_tri[j].y);
        }
        fclose(fp);
        sprintf(filename,"sample_quad_%s.dat",name[i]);
        fp = fopen(filename,"w");
        for (j=0; j<ns1; j++) {
            fprintf(fp,"%g %g\n",samples_quad[j].x,samples_quad[j].y);
        }
        fclose(fp);
    }
}  

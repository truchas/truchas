#ifndef VF_NO_MPI
#include <mpi.h>
#endif

#ifndef WIN32
#include <unistd.h>
#endif
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "vf.h"

#define NO_SUPPORT_HYBRID
#define NO_RADIOSITY_OUTPUT

#undef TRUE
#include "exodusII.h"

#if (__STDC__ == 1)     
int getopt(int argc, char * const *argv,  const  char  *optstring);
extern char *optarg;
extern int optind, opterr, optopt;
#endif

#define NextLine { while (1) { \
                     if (fgets(input_line,130,fp) == NULL) { \
                       fprintf(stderr,"FATAL ERROR: Premature EOF detected\n"); \
                       exit(0); \
                     } \
                     if (input_line[0]!='&') break; \
                   } \
                 }

void AllCaps(char*);

void 
Read_Mesh_Ascii(char *filename, int *ndim, int *nnodes, int *nfacets,
                int *npatches, double **x, double **y, double **z, int **c, 
                int **gid, int **map, int mysize, int myrank, int partial); 

void 
Read_Mesh_Exodus(char *filename, int *ndim, int *nnodes, int *nfacets,
                 int *npatches, double **x, double **y, double **z, int **c, 
                 int **gid, int **map, int mysize, int myrank, int partial); 

typedef struct ParamsStruct {
    char    enclID[100];
    int     geom;
    int     nonblocking;
    int     partial;
    double  asink;
    double  spatial_tol;
    int     xmirror;
    int     ymirror;
    int     zmirror;
    int     nrotations;
    int     bsp_depth;
    int     bsp_length;
    int     vf_method;
    int     vf_output;
    int     hc_res;
    int     hc_subs;
    double  hc_mind;
    int     vis_nsamples;
    int     vis_sampling;
    int     mc_nsamples;
    int     mc_sampling;
    double  mc_tol1;
    double  mc_tol2;
    int     smooth;
    double  sm_wt;
    double  sm_tol;
    int     sm_iter;
    int     sm_symm;
    int     sm_output;
    int     rad_solve;
    int     rad_solver;
    int     rad_dlb;
    int     rad_output;
    int     vf_nsteps;
    char**  mesh_file;
    int*    mesh_format;
    int*    vf_file_output;
    int*    vf_file_format;
    int*    gids;
    double* tsurf;
    double* radq;
    double* radj;
    double* eps;
} Params;

int main(int argc,char **argv)
{
  Params *params=NULL;
  FILE   *fp=NULL;
  char   title[250],tmp[250],*path=NULL;
  char   dir[132], param_file[250];
  char   gen_file[250], exo_file[250];
  char   inp_file[250], out_file[250];
  char   input_line[250], token0[40], token1[40], token2[40];
  int    mysize=1, myrank=0;
  int    jitter=1, randomize=1;
  int    i, j, k, n, cnt, proc, size=0, offset, opt, errflag=0;
  int    staged_io=0, ncontrollers=1, nsteps=0, step=0;
  int    nenclosures, *enclosure=NULL;
  int    ndim, nnodes, nfacets, npatches, npatches_g, npatches_r, max_patches;
  int    vertex_offset=1, exoid1;
  int    *c=NULL, *gid=NULL, *map=NULL;
  int    max_iter=500, num_iter;
  double sigma=5.7e-5, tol=1.0e-6;
  double *x=NULL, *y=NULL, *z=NULL;
  double radq, radj, tsurf, eps;
  extern char *optarg;
  extern int  optind;

#ifndef VF_NO_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&mysize);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  VF_Setup(staged_io, ncontrollers, MPI_COMM_WORLD);
#else
  VF_Setup(staged_io, ncontrollers, 0);
#endif

  if (argc==1) {
    errflag++;
  } else {
    while ((opt=getopt(argc, argv, "j:r:")) != EOF) {
      switch (opt) {
      case 'j':
        jitter = atoi(optarg);
        break;
      case 'r':
        randomize = atoi(optarg);
        break;
      case '?':
        errflag++;
        break;
      }
    }
  }
  if (errflag) {
    printf("\n");
    printf("usage:\tcmd [-j hemicube_jitter_flag] [-r random_surface_flag] ");
    printf("directory\n");
    printf("\n");
    printf("\t    hemicube_jitter_flag:  0=no jitter, 1=jitter\n");
    printf("\t    random_surface_flag:   0=no, 1=yes\n");
    printf("\n");
    exit(2);
  }
  sprintf(dir,"%s/",argv[optind]);

  if (mysize==1) {
    path = getenv("PWD");
  } else {
    memset(tmp,0,250);
    if (myrank==0) {
      path = getenv("PWD");
      strcpy(tmp,path);
      size = strlen(tmp);
      printf("<%s>\n",path);
    }
    VF_BroadcastInt(&size, 1, 0);
    VF_BroadcastChar(tmp, size, 0);
    path = tmp;
  }
  sprintf(param_file,"%s/%sparams.txt",path,dir);
  if ((fp=fopen(param_file,"r")) == NULL) {
    fprintf(stderr,"FATAL ERROR: Cannot open parameter file %s\n",param_file);
    exit(0);
  }
  
  NextLine;
  strcpy(title,input_line);
 
  NextLine;
  sscanf(input_line,"%d%d%d",&nsteps,&nenclosures,&max_patches);
  
  VF_SetNumEnclosures(nenclosures);
  VF_SetMaxSurfaces(max_patches);
  VF_OutputInitBanner();
  if (!jitter) {
    if (myrank==0) printf("   Turning hemicube jitter off\n\n");
    VF_JitterOff();
  }
  if (!randomize) {
    if (myrank==0) printf("   Turning surface randomization off\n\n");
    VF_RandomizeSurfacesOff();
  }
 
  enclosure = VF_Newi(nenclosures);
  params    = (Params*)VF_Newc(nenclosures*sizeof(Params));
  for (n=0; n<nenclosures; n++) {
    params[n].vf_nsteps      = 0;
    params[n].mesh_file      = NULL;
    params[n].mesh_format    = NULL;
    params[n].vf_file_output = NULL;
    params[n].vf_file_format = NULL;
    params[n].gids           = NULL;
    params[n].tsurf          = NULL;
    params[n].radq           = NULL;
    params[n].radj           = NULL;
    params[n].eps            = NULL;
  }

  for (n=0; n<nenclosures; n++) {
  
    NextLine;
    sscanf(input_line,"%s",params[n].enclID);
    
    NextLine;
    sscanf(input_line,"%s%s%s%lf%lf",
           token0,token1,token2,
           &params[n].asink,&params[n].spatial_tol);
    AllCaps(token0);
    AllCaps(token1);
    AllCaps(token2);
    if      (strcmp(token0,"2DAXISYM")==0) params[n].geom = VF_2Daxisym;
    else if (strcmp(token0,"2DPLANAR")==0) params[n].geom = VF_2Dplanar;
    else if (strcmp(token0,"3D"      )==0) params[n].geom = VF_3D;
    else {
      printf("FATAL ERROR: unrecognized geometry, %s\n",token0);
      exit(0);
    }
    if      (strcmp(token1,"NO" )==0) params[n].nonblocking = 0;
    else if (strcmp(token1,"YES")==0) params[n].nonblocking = 1;
    else {
      printf("FATAL ERROR: unrecognized nonblocking choice, %s\n",token1);
      exit(0);
    }
    if      (strcmp(token2,"NO" )==0) params[n].partial = 0;
    else if (strcmp(token2,"YES")==0) params[n].partial = 1;
    else {
      printf("FATAL ERROR: unrecognized partial choice, %s\n",token2);
      exit(0);
    }
    
    NextLine;
    sscanf(input_line,"%s%s%s%d",
           token0,token1,token2,&params[n].nrotations);
    AllCaps(token0);
    AllCaps(token1);
    AllCaps(token2);
    if      (strcmp(token0,"NO" )==0) params[n].xmirror = 0;
    else if (strcmp(token0,"YES")==0) params[n].xmirror = 1;
    else {
      printf("FATAL ERROR: unrecognized xmirror choice, %s\n",token0);
      exit(0);
    }
    if      (strcmp(token1,"NO" )==0) params[n].ymirror = 0;
    else if (strcmp(token1,"YES")==0) params[n].ymirror = 1;
    else {
      printf("FATAL ERROR: unrecognized ymirror choice, %s\n",token1);
      exit(0);
    }
    if      (strcmp(token2,"NO" )==0) params[n].zmirror = 0;
    else if (strcmp(token2,"YES")==0) params[n].zmirror = 1;
    else {
      printf("FATAL ERROR: unrecognized zmirror choice, %s\n",token2);
      exit(0);
    }
           
    
    NextLine;
    sscanf(input_line,"%d%d",&params[n].bsp_depth,&params[n].bsp_length);
 
    NextLine;
    sscanf(input_line,"%s%d",token0,&params[n].vf_output);
    AllCaps(token0);
    if      (strcmp(token0,"READ"    )==0) params[n].vf_method = VF_FILE_READ;
    else if (strcmp(token0,"PAIRWISE")==0) params[n].vf_method = VF_PAIRWISE;
    else if (strcmp(token0,"HEMICUBE")==0) params[n].vf_method = VF_HEMICUBE;
#ifdef SUPPORT_HYBRID
    else if (strcmp(token0,"HYBRID"  )==0) params[n].vf_method = VF_HYBRID;
#endif
    else {
      printf("FATAL ERROR: unrecognized VF calc method %s\n",token0);
      exit(0);
    }
 
    switch (params[n].vf_method) {
    case VF_HEMICUBE:
      NextLine;
      sscanf(input_line,"%d%d%lf",&params[n].hc_res,
             &params[n].hc_subs,&params[n].hc_mind);
      break;
    case VF_PAIRWISE:
      NextLine;
      sscanf(input_line,"%d%s%d%s%lf%lf",
             &params[n].vis_nsamples,token0,
             &params[n].mc_nsamples,token1,
             &params[n].mc_tol1,&params[n].mc_tol2);
      AllCaps(token0);
      AllCaps(token1);
      if      (strcmp(token0,"RANDOM" )==0) params[n].vis_sampling = VF_RANDOM_SAMPLE;
      else if (strcmp(token0,"UNIFORM")==0) params[n].vis_sampling = VF_UNIFORM_SAMPLE;
      else if (strcmp(token0,"JITTER" )==0) params[n].vis_sampling = VF_JITTER_SAMPLE;
      else if (strcmp(token0,"HALTON" )==0) params[n].vis_sampling = VF_HALTON_SAMPLE;
      else {
        printf("FATAL ERROR: unrecognized vis sampling method, %s\n",token0);
        exit(0);
      }
      if      (strcmp(token1,"RANDOM" )==0) params[n].mc_sampling  = VF_RANDOM_SAMPLE;
      else if (strcmp(token1,"UNIFORM")==0) params[n].mc_sampling  = VF_UNIFORM_SAMPLE;
      else if (strcmp(token1,"JITTER" )==0) params[n].mc_sampling  = VF_JITTER_SAMPLE;
      else if (strcmp(token1,"HALTON" )==0) params[n].mc_sampling  = VF_HALTON_SAMPLE;
      else {
        printf("FATAL ERROR: unrecognized mc sampling method, %s\n",token1);
        exit(0);
      }
      break;
#ifdef SUPPORT_HYBRID
    case VF_HYBRID:
      NextLine;
      sscanf(input_line,"%d%d%lf",&params[n].hc_res,
             &params[n].hc_subs,&params[n].hc_mind);
      NextLine;
      sscanf(input_line,"%d%s%d%s%lf%lf",
             &params[n].vis_nsamples,token0,
             &params[n].mc_nsamples,token1,
             &params[n].mc_tol1,&params[n].mc_tol2);
      AllCaps(token0);
      AllCaps(token1);
      if      (strcmp(token0,"RANDOM" )==0) params[n].vis_sampling = VF_RANDOM_SAMPLE;
      else if (strcmp(token0,"UNIFORM")==0) params[n].vis_sampling = VF_UNIFORM_SAMPLE;
      else if (strcmp(token0,"JITTER" )==0) params[n].vis_sampling = VF_JITTER_SAMPLE;
      else if (strcmp(token0,"HALTON" )==0) params[n].vis_sampling = VF_HALTON_SAMPLE;
      else {
        printf("FATAL ERROR: unrecognized vis sampling method, %s\n",token0);
        exit(0);
      }
      if      (strcmp(token1,"RANDOM" )==0) params[n].mc_sampling  = VF_RANDOM_SAMPLE;
      else if (strcmp(token1,"UNIFORM")==0) params[n].mc_sampling  = VF_UNIFORM_SAMPLE;
      else if (strcmp(token1,"JITTER" )==0) params[n].mc_sampling  = VF_JITTER_SAMPLE;
      else if (strcmp(token1,"HALTON" )==0) params[n].mc_sampling  = VF_HALTON_SAMPLE;
      else {
        printf("FATAL ERROR: unrecognized mc sampling method, %s\n",token1);
        exit(0);
      }
      break;
#endif
    }
           
    NextLine;
    sscanf(input_line,"%s%lf%lf%d%s%d",
           token0,&params[n].sm_wt,&params[n].sm_tol,
           &params[n].sm_iter,token1,&params[n].sm_output);
    AllCaps(token0);
    AllCaps(token1);
    if      (strcmp(token0,"NO" )==0) params[n].smooth = 0;
    else if (strcmp(token0,"YES")==0) params[n].smooth = 1;
    else {
      printf("FATAL ERROR: unrecognized smoothing choice, %s\n",token0);
      exit(0);
    }
    if      (strcmp(token1,"NONE")==0) params[n].sm_symm = VF_SYMMETRIC_NONE;
    else if (strcmp(token1,"SUB" )==0) params[n].sm_symm = VF_SYMMETRIC_SUB;
    else if (strcmp(token1,"ADD" )==0) params[n].sm_symm = VF_SYMMETRIC_ADD;
    else if (strcmp(token1,"AVG" )==0) params[n].sm_symm = VF_SYMMETRIC_AVG;
    else {
      printf("FATAL ERROR: unrecognized symmetry method, %s\n",token1);
      exit(0);
    }
    
    NextLine;
    sscanf(input_line,"%d",&params[n].vf_nsteps);
    
    params[n].mesh_format    = VF_Newi(params[n].vf_nsteps);
    params[n].mesh_file      = VF_NewcPtr(params[n].vf_nsteps);
    params[n].vf_file_output = VF_Newi(params[n].vf_nsteps);
    params[n].vf_file_format = VF_Newi(params[n].vf_nsteps);
    for (i=0; i<params[n].vf_nsteps; i++) {
      params[n].mesh_file[i] = VF_Newc(200);
      NextLine;
      sscanf(input_line,"%s%s%s%s",title,token0,token1,token2);
      AllCaps(token0);
      AllCaps(token1);
      AllCaps(token2);
      sprintf(params[n].mesh_file[i],"%s/%s%s",path,dir,title);
      if      (strcmp(token0,"ASCII" )==0) params[n].mesh_format[i] = VF_ASCII;
      else if (strcmp(token0,"EXODUS")==0) params[n].mesh_format[i] = VF_BINARY;
      else {
        printf("FATAL ERROR: unrecognized mesh file format, %s\n", token1);
        exit(0);
      }
      if      (strcmp(token1,"NO" )==0) params[n].vf_file_output[i] = 0;
      else if (strcmp(token1,"YES")==0) params[n].vf_file_output[i] = 1;
      else {
        printf("FATAL ERROR: unrecognized vf output file choice, %s\n", token1);
        exit(0);
      }
      if      (strcmp(token2,"ASCII" )==0) params[n].vf_file_format[i] = VF_ASCII;
      else if (strcmp(token2,"BINARY")==0) params[n].vf_file_format[i] = VF_BINARY;
      else if (strcmp(token2,"XDR"   )==0) params[n].vf_file_format[i] = VF_XDRFMT;
      else {
        printf("FATAL ERROR: unrecognized mesh file format, %s\n", token1);
        exit(0);
      }
    }
    
    NextLine;
    sscanf(input_line,"%s%s%s%d",
           token0,token1,token2,&params[n].rad_output);
    AllCaps(token0);
    AllCaps(token1);
    AllCaps(token2);
    if      (strcmp(token0,"NO" )==0) params[n].rad_solve = 0;
    else if (strcmp(token0,"YES")==0) params[n].rad_solve = 1;
    else {
      printf("FATAL ERROR: unrecognized radiosity solve choice\n");
      exit(0);
    }
    if      (strcmp(token1,"CG"     )==0) params[n].rad_solver = VF_RADSOLVE_CG;
    else if (strcmp(token1,"GMRES"  )==0) params[n].rad_solver = VF_RADSOLVE_GMRES;
    else if (strcmp(token1,"AZ_CG"  )==0) params[n].rad_solver = VF_RADSOLVE_AZTEC_CG;
    else if (strcmp(token1,"AZGMRES")==0) params[n].rad_solver = VF_RADSOLVE_AZTEC_GMRES;
    else {
      printf("FATAL ERROR: unrecognized radiosity solver\n");
      exit(0);
    }
    if      (strcmp(token2,"NO" )==0) params[n].rad_dlb = 0;
    else if (strcmp(token2,"YES")==0) params[n].rad_dlb = 1;
    else {
      printf("FATAL ERROR: unrecognized dlb choice\n");
      exit(0);
    }
    if (mysize==1) params[n].rad_dlb = 0;
    
  }  /* end of enclosure param loop */
   
  srand(325467);
  for (step=0; step<nsteps; step++) { 
   
    if (myrank==0) {
      printf("\n=======================================");
      printf("=======================================\n");
      printf(" S T E P :   %d\n",step+1);
      printf("=======================================");
      printf("=======================================\n");
    }
    for (n=0; n<nenclosures; n++) { 
      
      if (step<params[n].vf_nsteps) {
      
        if (params[n].mesh_format[step]==0) {
          Read_Mesh_Ascii(params[n].mesh_file[step], &ndim, &nnodes, &nfacets, 
                          &npatches, &x, &y, &z, &c, &gid, &map, mysize, myrank, params[n].partial);
        } else {
          Read_Mesh_Exodus(params[n].mesh_file[step], &ndim, &nnodes, &nfacets, 
                           &npatches, &x, &y, &z, &c, &gid, &map, mysize, myrank, params[n].partial);
        }
        npatches_g = npatches;
        VF_GlobalSumInt(&npatches_g);
        if (step==0) {
          params[n].gids  = VF_Newi(npatches_g);
          params[n].tsurf = VF_Newd(npatches_g);
          params[n].radq  = VF_Newd(npatches_g);
          params[n].radj  = VF_Newd(npatches_g);
          params[n].eps   = VF_Newd(npatches_g);
        } else {
          params[n].gids  = VF_ReNewi(params[n].gids, npatches_g);
          params[n].tsurf = VF_ReNewd(params[n].tsurf,npatches_g);
          params[n].radq  = VF_ReNewd(params[n].radq, npatches_g);
          params[n].radj  = VF_ReNewd(params[n].radj, npatches_g);
          params[n].eps   = VF_ReNewd(params[n].eps,  npatches_g);
        }
        if (mysize>1) {
          for (size=0, proc=0; proc<mysize; proc++) {
            if (proc==myrank) {
              cnt = npatches;
              for (i=0; i<npatches; i++) {
                params[n].gids[size+i] = gid[i];
              }
            }
            VF_BroadcastInt(&cnt, 1, proc);
            VF_BroadcastInt(&(params[n].gids[size]), cnt, proc);
            size += cnt;
          }
          /*VF_SortIntegerArray(params[n].gids, npatches_g);*/
        } else {
          for (i=0; i<npatches; i++) {
            params[n].gids[i] = gid[i];
          }
        }

        sprintf(gen_file,"%s/%s%s_step%d_topology.gen",
                path,dir,params[n].enclID,step);
        sprintf(exo_file,"%s/%s%s_step%d_enclosure.exo",
                path,dir,params[n].enclID,step);
        if (mysize==1) {
          sprintf(inp_file,"%s/%s%s_step%d_input_matrix.vf",
                  path,dir,params[n].enclID,step);
          sprintf(out_file,"%s/%s%s_step%d_output_matrix.vf",
                  path,dir,params[n].enclID,step);
        } else {
          sprintf(inp_file,"%s/%s%s_step%d_input_matrix.vf.%d.%d",
                  path,dir,params[n].enclID,step,mysize,myrank);
          sprintf(out_file,"%s/%s%s_step%d_output_matrix.vf.%d.%d",
                  path,dir,params[n].enclID,step,mysize,myrank);
        }
                         
        if (params[n].vf_method==VF_FILE_READ) {
          enclosure[n] = VF_NewEnclosure(params[n].enclID,
                                         params[n].vf_output);
          VF_CalcRead(enclosure[n], inp_file);
        } else {
          enclosure[n] = VF_DefineEnclosure(params[n].enclID,
                                            params[n].nonblocking, 
                                            params[n].partial, 
                                            params[n].asink, 
                                            npatches, gid, 
                                            params[n].vf_output);
          VF_DefineTopology(enclosure[n], params[n].geom, nfacets, nnodes, 
                            x, y, z, c, vertex_offset, map, 
                            params[n].nrotations, 
                            params[n].xmirror, 
                            params[n].ymirror, 
                            params[n].zmirror, 
                            params[n].bsp_depth, 
                            params[n].bsp_length, 
                            params[n].spatial_tol, 
                            params[n].vf_output);
          switch (params[n].vf_method) {
          case VF_PAIRWISE:
            VF_CalcPairwise(enclosure[n], 
                            params[n].vis_nsamples, params[n].vis_sampling,
                            params[n].mc_nsamples, params[n].mc_sampling, 
                            params[n].mc_tol1, params[n].mc_tol2);
            break;
          case VF_HEMICUBE:
            VF_CalcHemicube(enclosure[n], 
                            params[n].hc_subs, 
                            params[n].hc_res, 
                            params[n].hc_mind);
            break;
#ifdef SUPPORT_HYBRID
	  case VF_HYBRID:
            VF_CalcHybrid(enclosure[n], 
                          params[n].hc_subs, 
                          params[n].hc_res, 
                          params[n].hc_mind,
                          params[n].vis_nsamples, params[n].vis_sampling,
                          params[n].mc_nsamples, params[n].mc_sampling, 
                          params[n].mc_tol1, params[n].mc_tol2);
            break;
#endif
	  }
        }
        if (params[n].smooth) {
          VF_SmoothMatrix(enclosure[n], 
                          params[n].sm_wt, 
                          params[n].sm_tol,
                          params[n].sm_iter, 
                          params[n].sm_symm, 
                          params[n].sm_output);
        }
        VF_OutputMatrixSummaryBanner();
        if (params[n].vf_method!=VF_FILE_READ) {
          if (myrank==0) printf("   VF_WriteGenesis(%s)\n",gen_file);
          VF_WriteGenesis(gen_file);
          if (myrank==0) printf("   VF_WriteExodus (%s)\n",exo_file);
#ifdef RADIOSITY_OUTPUT
          exoid1 = VF_WriteExodusAux1(exo_file, 1, 0.0);
#else
          VF_WriteExodus(exo_file);
#endif
        }
        if (params[n].vf_method!=VF_FILE_READ) {
          VF_ResetTopology(enclosure[n]);
        }
        
        VF_Free(x);
        VF_Free(y);
        VF_Free(z);
        VF_Free(c);
        VF_Free(gid);
        VF_Free(map);
        if (params[n].vf_file_output[step]) {
          if (myrank==0) printf("   VF_MatrixWrite (%s, %d)\n",
                                out_file, params[n].vf_file_format[step]);
          VF_MatrixWrite(enclosure[n], out_file, params[n].vf_file_format[step]);
        }
      }
      if (params[n].rad_solve) {
        if (step>0 && params[n].rad_dlb) {
          npatches_g = VF_NumPatches_g(enclosure[n]);
          for (i=npatches_g-1; i>1; i--) {
            j                  = (i-1)*((double)rand()/(double)RAND_MAX);
            k                  = params[n].gids[j];
            tsurf              = params[n].tsurf[j];
            radq               = params[n].radq[j];
            radj               = params[n].radj[j];
            eps                = params[n].eps[j];
            params[n].gids[j]  = params[n].gids[i];
            params[n].tsurf[j] = params[n].tsurf[i];
            params[n].radq[j]  = params[n].radq[i];
            params[n].radj[j]  = params[n].radj[i];
            params[n].eps [j]  = params[n].eps[i];
            params[n].gids[i]  = k;
            params[n].tsurf[i] = tsurf;
            params[n].radq[i]  = radq;
            params[n].radj[i]  = radj;
            params[n].eps[i]   = eps;
          }
          if (myrank%2) {
            npatches += 5;
          }else {
            if (myrank!=mysize-1) npatches -= 5;
          }
          for (offset=0, proc=0; proc<mysize-1; proc++) {
            cnt = npatches;
            VF_BroadcastInt(&cnt, 1, proc);
            if (myrank>proc) offset+=cnt;
          }
          VF_RemapEnclosure (enclosure[n], npatches, &(params[n].gids[offset]));
        } else {
          npatches   = VF_NumPatches_l(enclosure[n]);
          npatches_g = VF_NumPatches_g(enclosure[n]);
          npatches_r = npatches_g/4;
          for (j=0; j<npatches_r; j++) {
            params[n].tsurf[j] = 400.0;
          }
          for (j=npatches_r; j<npatches_g; j++) {
            params[n].tsurf[j] = 300.0;
          }
          for (offset=0, proc=0; proc<mysize-1; proc++) {
            cnt = npatches;
            VF_BroadcastInt(&cnt, 1, proc);
            if (myrank>proc) offset+=cnt;
          }
        }
        for (j=0; j<npatches_g; j++) {
          params[n].radq[j] = 0.0;
        }
        for (j=0; j<npatches_g; j++) {
          params[n].eps[j] = 0.8;
        }
#ifdef RADIOSITY_OUTPUT
        VF_RadSolveAux(enclosure[n], 
                       &(params[n].radq[offset]), 
                       &(params[n].radj[offset]),
                       &(params[n].tsurf[offset]), 
                       &(params[n].eps[offset]), 
                       sigma, tol, 
                       max_iter, &num_iter, 
                       params[n].rad_solver, 0, 
                       params[n].rad_output);
        VF_WriteExodusAux2(enclosure[n],exoid1,1,0.0,
                           &(params[n].radq[offset]),
                           &(params[n].radj[offset]));
#else
        VF_RadSolve(enclosure[n], 
                    &(params[n].radq[offset]), 
                    &(params[n].tsurf[offset]), 
                    &(params[n].eps[offset]), 
                    sigma, tol, 
                    max_iter, &num_iter, 
                    params[n].rad_solver, 0, 
                    params[n].rad_output);
#endif
      }
    }
  }
  VF_CleanUp();
  
  for (n=0; n<nenclosures; n++) {
    VF_Free(params[n].mesh_format);
    VF_Free(params[n].vf_file_output);
    VF_Free(params[n].vf_file_format);
    VF_Free(params[n].gids);
    VF_Free(params[n].tsurf);
    VF_Free(params[n].radq);
    VF_Free(params[n].radj);
    VF_Free(params[n].eps);
    for (i=0; i<params[n].vf_nsteps; i++) {
      VF_Free(params[n].mesh_file[i]);
    }
    VF_Free(params[n].mesh_file);
  }
  VF_Free(params);
  VF_Free(enclosure);

#ifndef VF_NO_MPI
  fflush(stdout);
  sleep(3);
  MPI_Finalize();
#endif
  exit(0);
}

void AllCaps(char *s)
{
  int i, n;

  n = strlen(s);
  for (i=0; i<n; i++) {
    s[i] = toupper(s[i]);
  }
}

void 
Read_Mesh_Ascii(char *filename, int *ndim, int *nnodes, int *nfacets,
                int *npatches, double **x, double **y, double **z, int **c, 
                int **gid, int **map, int mysize, int myrank, int partial)
    
{
  FILE   *fp;
  char   input_line[250], *token;
  int    i, j, m, n1, n2;
  int    patch_id, num_nodes, num_facets;
  int    *cc=NULL, *gg=NULL, *mm=NULL;
  double *xx=NULL, *yy=NULL, *zz=NULL;
  
  if ((fp=fopen(filename,"r")) == NULL) {
    fprintf(stderr,"FATAL ERROR: Cannot open mesh file %s\n",filename);
    exit(0);
  }
  NextLine;
  sscanf(input_line,"%d%d%d",ndim,nnodes,nfacets);
  num_nodes  = *nnodes;
  num_facets = *nfacets;

  xx = VF_Newd(num_nodes);
  yy = VF_Newd(num_nodes);
  zz = VF_Newd(num_nodes);
  for (i=0; i<num_nodes; i++) {
    xx[i] = 0.0;
    yy[i] = 0.0;
    zz[i] = 0.0;
  }

  /* READ THE NODE COORDINATES */
  for (i=0; i<num_nodes; i++) {
    NextLine;
    token = strtok(input_line, " \t"); /* node_id */
    if (token=strtok(NULL, " \t")) xx[i]=atof(token);else continue;
    if (token=strtok(NULL, " \t")) yy[i]=atof(token);else continue;
    if (token=strtok(NULL, " \t")) zz[i]=atof(token);else continue;
  }
  *x = xx;
  *y = yy;
  *z = zz;

  /* READ THE FACET CONNECTIVITIES */
  cc = VF_Newi(4*num_facets);
  gg = VF_Newi(1+num_facets);
  mm = VF_Newi(1+num_facets);
  for (i=0; i<4*num_facets; i++) cc[i] = -1;
  for (i=0; i<num_facets; i++) {
    m = 0;
    j = 4*i;
    NextLine;
    token    = strtok(input_line, " \t"); /* facet_id */
    token    = strtok(NULL, " \t");       /* patch_id */
    patch_id = atoi(token);
    token    = strtok(NULL, " \t");
    while (token!=NULL) {
      cc[j+m] = atoi(token); m++;
      token  = strtok(NULL, " \t");
    }
    gg[i] = patch_id;
    mm[i] = patch_id;
  }
  fclose(fp);
  
  *npatches = 1;
  VF_SortIntegerArray(gg,num_facets);
  for (j=1, i=1; i<num_facets; i++) {
    if (gg[i]!=gg[i-1]) {
      if (i!=j) gg[j] = gg[i];
      *npatches += 1;
      j++;
    }
  }
  if (partial) {
    gg[*npatches] = 9999999;
    *npatches += 1;
  }
  if (mysize>1) {
    n1 = *npatches/mysize;
    n2 = *npatches%mysize;
    if (myrank==mysize-1) {
      *npatches = n1+n2;
    } else {
      *npatches = n1;
    }
    *gid = &gg[myrank*n1];
    memcpy(gg,&gg[myrank*n1],*npatches*sizeof(int));
    
    n1 = num_facets/mysize;
    n2 = num_facets%mysize;
    if (myrank==mysize-1) {
      *nfacets = n1+n2;
    } else {
      *nfacets = n1;
    }
    *map = &mm[myrank*n1];
    *c   = &cc[myrank*n1*4];
    memcpy(mm,&mm[myrank*n1],  *nfacets*sizeof(int));
    memcpy(cc,&cc[myrank*n1*4],*nfacets*4*sizeof(int));
  }
  *c   = cc;
  *gid = gg;
  *map = mm;
} 

void 
Read_Mesh_Exodus(char *fname, int *ndim, int *nnodes, int *nfacets,
                 int *npatches, double **x, double **y, double **z, int **c, 
                 int **gid, int **map, int mysize, int myrank, int partial)
{
  char   title[250], filename[256],*meshfile,*dirname;
  int    ws1, ws2, exoid, exoerr;
  int    i, j, k, m, nblocks, nns, nss, nelems, npe, nattrs, *c0, *blk_id;
  float  version;
  int    *cc, *gg, *mm, num_nodes, num_facets;
  double *xx, *yy, *zz;
  
  if (mysize>1) {
    dirname     = fname;
    meshfile    = strrchr(fname,'/');
    *meshfile++ = (char)NULL;
    meshfile[strlen(meshfile)-2] = (char)NULL;
    sprintf(filename,"%s/01/%s.par.%d.%d",dirname,meshfile,mysize,myrank);
  } else {
    strcpy(filename, fname);
  }
  
  ws1    = sizeof(double);
  ws2    = 0;
  exoid  = ex_open(filename, EX_READ, &ws1, &ws2, &version);
  exoerr = ex_get_init(exoid,title,ndim,nnodes,
                       nfacets,&nblocks,&nns, &nss);
  num_nodes  = *nnodes;
  num_facets = *nfacets;

  /* READ THE NODE COORDINATES */
  xx = VF_Newd(num_nodes);
  yy = VF_Newd(num_nodes);
  zz = VF_Newd(num_nodes);
  for (i=0; i<num_nodes; i++) {
    xx[i] = 0.0;
    yy[i] = 0.0;
    zz[i] = 0.0;
  }
  exoerr = ex_get_coord(exoid,xx,yy,zz);
  
  gg     = VF_Newi(1+num_facets);
  mm     = VF_Newi(1+num_facets);
  exoerr = ex_get_elem_num_map(exoid, gg);
  exoerr = ex_get_elem_num_map(exoid, mm);

  /* READ THE FACET CONNECTIVITIES */
  cc     = VF_Newi(4*num_facets);
  c0     = VF_Newi(4*num_facets);
  blk_id = VF_Newi(nblocks);
  exoerr = ex_get_elem_blk_ids (exoid, blk_id);
  for (i=0; i<4*num_facets; i++) cc[i] = -1;
  for (m=0, i=0; i<nblocks; i++) {
    exoerr = ex_get_elem_block (exoid, blk_id[i], title,
                                &nelems, &npe, &nattrs);
    exoerr = ex_get_elem_conn(exoid,blk_id[i],c0);
    for (j=0; j<nelems; j++) {
      for (k=0; k<npe; k++) {
        cc[4*m+k] = c0[npe*j+k];
      }
      m++;
    }
  }
  VF_Free(c0);
  VF_Free(blk_id);
  ex_close(exoid);
  *x   = xx;
  *y   = yy;
  *z   = zz;
  *c   = cc;
  *gid = gg;
  *map = mm;
  if (exoerr) exoerr=0;
  
  *npatches = *nfacets;
  if (partial && myrank==mysize-1) {
    gg[*npatches] = 9999999;
    *npatches += 1;
  }
  
}

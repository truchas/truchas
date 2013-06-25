
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "netcdf.h"
#include "exodusII.h"

void main(int argc, char *argv[])
{

    FILE  *fp;
    char  *infile=NULL, *outfile=NULL;
    char  title[81], elem_type[81];
    int   cpuws=0, iows=0, size;
    int   exoid, num_dim, num_nodes, num_elems;
    int   num_elem_blocks, num_node_sets, num_side_sets;
    int   num_elem_in_block, nnodes_elem, num_attr;
    int   i, j, k, elem, block, *blk_id, *connect, conn_size=0;
    float *x, *y, *z, version;

    if (argc==3) {
        infile  = argv[1];
        outfile = argv[2];
    } else {
        printf("Usage:  exo2inp exo_file inp_file\n");
        exit(1);
    }
    if ((exoid=ex_open(infile,EX_READ,&cpuws,&iows,&version)) < 0) {
        printf("Cannot open .exo file %s\n",infile);
        exit(1);
    }
    if ((fp=fopen(outfile,"w")) == NULL) {
        fprintf(stderr,"cannot open .inp file %s",outfile);
        exit(1);
    }
    ex_get_init(exoid, title, &num_dim, &num_nodes, &num_elems,
                &num_elem_blocks, &num_node_sets, &num_side_sets);
    
    printf("In File: %s\n",infile);
    printf("cpuws:   %d\n",cpuws);
    printf("iows:    %d\n",iows);
    printf("version: %f\n",version);
    
    printf("Title: %s\n",title);
    printf("num_dim:          %d\n",num_dim);
    printf("num_nodes:        %d\n",num_nodes);
    printf("num_elems:        %d\n",num_elems);
    printf("num_elem_blocks:  %d\n",num_elem_blocks);
    printf("num_node_sets:    %d\n",num_node_sets);
    printf("num_side_sets:    %d\n",num_side_sets);
                
    fprintf(fp,"%s\n",title);
    fprintf(fp,"&\n");
    fprintf(fp,"& method: PAIRWISE=0, HEMICUBE=1\n");
    fprintf(fp,"& sample: RANDOM=0, UNIFORM=1, JITTER=2, HALTON=3\n");
    fprintf(fp,"& smooth: NONE=0, SIMPLE=1, LSQ=2, NORMALIZE=3\n");
    fprintf(fp,"& symmetry:  AUTO=-1, NONE=0, SUB=1, ADD=2, AVG=3\n");
    fprintf(fp,"&\n");
    fprintf(fp,"& ndim nnodes nfacets\n");
    fprintf(fp,"   %d  %d  %d\n",num_dim,num_nodes,num_elems);
    fprintf(fp,"& nodes:  node_id  coordinates\n");
 
    x = (float *)malloc(sizeof(float) * num_nodes);
    y = (float *)malloc(sizeof(float) * num_nodes);
    z = (float *)malloc(sizeof(float) * num_nodes);
    ex_get_coord(exoid, x, y, z);
    
    
    fprintf(fp,"&\n");
    for (i=0; i<num_nodes; i++) {
        switch (num_dim) {
        case 1:
            fprintf(fp," %6d\t%.10e\n",i+1,x[i]);
            break;
        case 2:
            fprintf(fp," %6d\t%.10e\t%.10e\n",i+1,x[i],y[i]);
            break;
        case 3:
            fprintf(fp," %6d\t%.10e\t%.10e\t%.10e\n",i+1,x[i],y[i],z[i]);
            break;
        }    
    }
    fprintf(fp,"&\n");
        
        
    blk_id = (int *)malloc(sizeof(int) * num_elem_blocks);
    ex_get_elem_blk_ids (exoid, blk_id);
    for (elem=1, block=0; block<num_elem_blocks; block++) {
        ex_get_elem_block (exoid, blk_id[block], elem_type,
                           &num_elem_in_block, &nnodes_elem, &num_attr);
        size = nnodes_elem*num_elem_in_block*sizeof(int);
        if (conn_size==0) {
            conn_size = size;
            connect   = (int *)malloc(conn_size);
        } else if (conn_size<size) {
            conn_size = size;
            connect   = (int *)realloc(connect,conn_size);
        }
        ex_get_elem_conn(exoid,blk_id[block],connect);
        
        
        for (i=0; i<num_elem_in_block; i++) {
            fprintf(fp," %6d",elem);
            for (j=0; j<nnodes_elem; j++) {
                fprintf(fp,"\t%6d",connect[i*nnodes_elem+j]);
            }
            fprintf(fp,"\n");
            elem++;
        }
    }
    
    ex_close(exoid);
    fclose(fp);
}

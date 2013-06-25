/* Header file for all files in scan (PREFIX and SUFFIX) routines */

/* $Id: scan.h,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */

/* Prototypes for scan routines */

#ifndef SCAN_H__
#define SCAN_H__

void off_node_sum_prefix_int_c   (int *Dest_Data, int *Src_Data, int *Count);
void off_node_sum_prefix_float_c (float *Dest_Data, float *Src_Data, int *Count);
void off_node_sum_prefix_double_c(double *Dest_Data, double *Src_Data, int *Count);

void off_node_sum_suffix_int_c   (int *Dest_Data, int *Src_Data, int *Count);
void off_node_sum_suffix_float_c (float *Dest_Data, float *Src_Data, int *Count);
void off_node_sum_suffix_double_c(double *Dest_Data, double *Src_Data, int *Count);

#endif

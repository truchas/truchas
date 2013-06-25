/* Header file for all files in shift */

/* $Id: shift.h,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */

/* Prototypes for shift routines */

#ifndef SHIFT_H__
#define SHIFT_H__

void physical_c_shift_up_int_c   (int *Dest_Data, int *Src_Data, int *Count);
void physical_c_shift_up_float_c (float *Dest_Data, float *Src_Data, int *Count);
void physical_c_shift_up_double_c(double *Dest_Data, double *Src_Data, int *Count);
void physical_c_shift_up_log_c   (C_LOG_TYPE *Dest_Data, C_LOG_TYPE *Src_Data, int *Count);

void physical_c_shift_down_int_c   (int *Dest_Data, int *Src_Data, int *Count);
void physical_c_shift_down_float_c (float *Dest_Data, float *Src_Data, int *Count);
void physical_c_shift_down_double_c(double *Dest_Data, double *Src_Data, int *Count);
void physical_c_shift_down_log_c   (C_LOG_TYPE *Dest_Data, C_LOG_TYPE *Src_Data, int *Count);

#endif

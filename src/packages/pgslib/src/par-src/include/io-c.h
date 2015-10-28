/* Header for io-c.c */

/* $Id: io-c.h,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $*/

#ifndef IO_H__
#define IO_H__

/* Prototypes for routines in io-c.c */
void pgslib_bcast_int8_scalar_c(int8_t *);
void pgslib_bcast_int_scalar_c(int *);
void pgslib_bcast_real_scalar_c(float *);
void pgslib_bcast_double_scalar_c(double *);
void pgslib_bcast_log_scalar_c(C_LOG_TYPE *);
void pgslib_bcast_char_scalar_c(char *);

void pgslib_bcast_int8_vector_c(int8_t *, int *);
void pgslib_bcast_int_vector_c(int *, int *);
void pgslib_bcast_real_vector_c(float *, int *);
void pgslib_bcast_double_vector_c(double *, int *);
void pgslib_bcast_log_vector_c(C_LOG_TYPE *, int *);
void pgslib_bcast_char_vector_c(char *, int *);

void pgslib_dist_int8_scalar_c(int8_t *, int8_t *);
void pgslib_dist_int_scalar_c(int *, int *);
void pgslib_dist_real_scalar_c(float *, float *);
void pgslib_dist_double_scalar_c(double *, double *);
void pgslib_dist_log_scalar_c(C_LOG_TYPE *, C_LOG_TYPE *);

void pgslib_dist_int8_vector_c(int8_t *, int *, int8_t *, int *);
void pgslib_dist_int_vector_c(int *, int *, int *, int *);
void pgslib_dist_real_vector_c(float *, int *, float *, int *);
void pgslib_dist_double_vector_c(double *, int *, double *, int *);
void pgslib_dist_log_vector_c(C_LOG_TYPE *, int *, C_LOG_TYPE *, int *);

void pgslib_collate_int8_scalar_c(int8_t *, int8_t *);
void pgslib_collate_int_scalar_c(int *, int *);
void pgslib_collate_real_scalar_c(float *, float *);
void pgslib_collate_double_scalar_c(double *, double *);
void pgslib_collate_log_scalar_c(C_LOG_TYPE *, C_LOG_TYPE *);

void pgslib_collate_int8_vector_c(int8_t *, int *, int8_t *, int *);
void pgslib_collate_int_vector_c(int *, int *, int *, int *);
void pgslib_collate_float_vector_c(float *, int *, float *, int *);
void pgslib_collate_double_vector_c(double *, int *, double *, int *);
void pgslib_collate_log_vector_c(C_LOG_TYPE *, int *, C_LOG_TYPE *, int *);

#endif

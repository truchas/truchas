/* *************************************************************************** *
*                                                                              *
*                                                                              *
*                             Copyright  (C) 20xx,                             *
*                      Los Alamos National Security, LLC                       *
*                                                                              *
*                             LA-CC-xxxxxx                                     *
*                                                                              *
* **************************************************************************** */

/*
* danu_attribute.h
*
*  DANU scalar attributes
*
*/
#include<string.h>
#include <ctype.h>

#include <hdf5.h>

#include <danu_error.h>
#include <danu_fort_strings.h>

/* Fortran String Data */
int compute_fortran_trim_len(const char * fort_data, const int f_len)
{
    int trim_len = f_len;
    const char *ptr;
    char c;

    ptr=fort_data;
    ptr+=trim_len-1;
    c = *ptr;
    while ( ptr >= fort_data && isblank(c) ) {
	ptr--;
	trim_len--;
	c = *ptr;
    }

    return trim_len;
}

danu_err_t convert_string_f2c(const char *fort_str, int f_len, char * c_str, int c_len)
{
    danu_err_t status = DANU_SUCCESS;
    char *p;
    int cp_len;

    cp_len = MIN(f_len,c_len);
    memcpy(c_str,fort_str,cp_len);
    p = c_str;
    p+=cp_len;
    *p = '\0';
    p--;
    while ( p >= c_str && isblank(*p) ) {
	*p = '\0';
	p--;
    }

    if ( c_len < (f_len+1) ) {
	DANU_WARN_MESS("C string buffer too small to hold fortran string");
	status = DANU_FAILURE;
    }

    return status;
}

danu_err_t convert_string_c2f(const char *c_str, char *fort_str, int f_len)
{
    danu_err_t status = DANU_SUCCESS;
    int cp_len;
    char *p;
    int c_len = strlen(c_str);

    /* Wipe the fortran string with spaces frst */
    memset(fort_str,' ',f_len);

    /* Now copy */
    cp_len = MIN(c_len,f_len);
    memcpy(fort_str,c_str,cp_len);

    if ( c_len > f_len ) {
      DANU_WARN_MESS("Fortran character string buffer too small");
      status = DANU_FAILURE;
    }

    return status;
}

danu_err_t copy_c_array_to_fortran_array(const char * const *c_array, int num,
	                                 char *fort_data, const int *fnum, const int *flen)
{
    danu_err_t status = DANU_SUCCESS;
    int cp_len;
    int str_len;
    int cp_num;
    char *p;
    int i;
    size_t bytes;

    /* Initialize the fort_data to whitespace */
    bytes = (*flen)*(*fnum)*sizeof(char);
    memset(fort_data,' ',bytes);

    /* Now copy data to each point in the fort_data array */
    p = fort_data;
    cp_num = MIN(*fnum,num);
    for(i=0; i< cp_num; i++) {
        str_len = strlen((const char *)c_array[i]);
	cp_len = MIN(str_len,*flen);
	memcpy(p,c_array[i],cp_len);
	p+=(*flen);
    }

    return status;

}

char ** copy_fortran_array_to_c_array(const char *fort_data, const int *fnum, const int *flen)
{
    danu_err_t status = DANU_SUCCESS;
    char **ptr;
    const char *p;
    int i;
    int trim_len, cp_len;
    
    ptr = DANU_MALLOC(char *, *fnum);
    p = fort_data;
    for(i=0; i<*fnum; i++) {

	trim_len = compute_fortran_trim_len(p,*flen);
	ptr[i] = DANU_MALLOC(char, *flen+1);
	cp_len = MIN(trim_len,*flen);
	memset(ptr[i],'\0',*flen+1);
	memcpy(ptr[i],p,cp_len);

	p+=(*flen);
    }

    return ptr;

}




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
* danu_errno.h
*
*  DANU Error Codes
*
*/

#ifndef DANU_ERRNO_H
#define DANU_ERRNO_H

enum {
   /* DANU error codes */
   DANU_SUCCESS    = 0  ,  /* Success */

   DANU_FAILURE         ,  /* Generic failure */
   DANU_EFAULT          ,  /* invalid pointer */
   DANU_EINVAL          ,  /* invalid argument supplied by user */
   DANU_ENOMEM          ,  /* malloc failed */
   DANU_ETABLE		    ,  /* Table error, i.e. exceed a table limit */
   DANU_EBADO           ,  /* Invalid object (group, dataset, attribute) */
   DANU_EACCES          ,  /* Required permissions are denied */
   DANU_EEXIST          ,  /* Existence error, i.e. attempt access object or file
							  that DNE */
   DANU_EBADF           ,  /* Invalid file descriptor */
   DANU_EOF             ,  /* end of file */

};

#endif /* DANU_ERRNO_H */


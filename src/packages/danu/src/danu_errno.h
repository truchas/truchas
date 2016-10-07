/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/
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


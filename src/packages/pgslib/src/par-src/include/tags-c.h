/* These routines and constants support the tags used in the MPI calls
   in PGSLib.*/

/* $Id: tags-c.h,v 1.1.1.1 2000/10/11 22:44:24 ferrell Exp $ */

#ifndef TAGS_H__
#define TAGS_H__

int initializeTagBits();
int constrainedSendRcv_TagBits();
int scan_TagBits();

#endif

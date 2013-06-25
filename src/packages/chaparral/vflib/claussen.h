/*
 * $Id: claussen.h,v 1.1 2005/04/08 20:17:03 mwglass Exp $
 * $Source: /cvsroot/chaparral/chaparral/vflib/claussen.h,v $
 *
 * Revision 1000.2  93/12/05  15:03:07  ps
 * replaced plane test with a new one and made all mach[] tests much tighter.
 * Also fixed a but in evalint{m|p}()
 * 
 * Revision 1000.1  93/06/18  06:56:08  ps
 * Bump to release, no file changes
 * 
 * Revision 1.2  93/04/14  22:04:04  ps
 * header file for claussen's function
 * 
 * Revision 1.1  93/04/06  21:25:48  ps
 * Initial revision
 * 
 */
 
#ifndef _CLAUSSEN_H_

#define _CLAUSSEN_H_

static  char    rcs_id_claussen_h[] = "$Header: /cvsroot/chaparral/chaparral/vflib/claussen.h,v 1.1 2005/04/08 20:17:03 mwglass Exp $";

extern  double  mach[5];
extern  double  claussen( double x );

#endif  /* CLAUSSEN_H */

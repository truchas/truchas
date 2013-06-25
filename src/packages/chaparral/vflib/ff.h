/*
 * $Id: ff.h,v 1.1 2005/04/08 20:17:03 mwglass Exp $
 * $Source: /cvsroot/chaparral/chaparral/vflib/ff.h,v $
 *
 * Revision 1000.4  93/12/28  18:23:17  ps
 * removed not used warning by wrapping up rcs string
 * 
 * Revision 1000.3  93/12/05  15:03:25  ps
 * replaced plane test with a new one and made all mach[] tests much tighter.
 * Also fixed a but in evalint{m|p}()
 * 
 * Revision 1000.2  93/11/28  13:14:37  ps
 * moved ALLCOEFF to ff.h
 * 
 * Revision 1000.1  93/06/18  06:56:09  ps
 * Bump to release, no file changes
 * 
 * Revision 1.1  93/04/14  22:03:29  ps
 * Initial revision
 * 
 */
#ifndef _FF_H_

#define _FF_H_

static char
__rcs_ref_ff_h( void )
{
  static char rcs_id_h[] = "$Header: /cvsroot/chaparral/chaparral/vflib/ff.h,v 1.1 2005/04/08 20:17:03 mwglass Exp $";
  return rcs_id_h[0];
}

extern  double  IDilog( double z[2] );
extern  double  G( double a, double b, double c, double t );
extern  double  H( double a, double b, double c, double t );
extern  void    Pair( double ( *c )[2],
                      double p1[3], double p2[3],
                      double q1[3], double q2[3] );
extern  int     Bilinearf( float ( *c )[2] );
extern  int     Bilinear( double ( *c )[2] );
extern  double  IntegralPlanar( double ( *c )[2] );
extern  int     Lcisf( float xax[2], float yax[2],
                      float rad[2], float cnt[2] );
extern  int     Lcis( double xax[2], double yax[2],
                      double rad[2], double cnt[2] );
extern  int     LogSelect( double rad[2], double cnt[2], double psi[2] );
extern  double  RM( double z[2] );
extern  double  IDilogPath( int k, double rad[2], double cnt[2],
                            double ( *c )[2] );
extern  double  ILogPart( int p1, double z[2], double ( *c )[2] );
extern  double  ILogIntegral( double ( *c )[2], double s );
extern  double  Integral( double ( *c )[2] );
extern  double  FormFactor( double ( *p )[3], int np,
                            double ( *q )[3], int nq );

extern  int     fferror;
#define NO_FF_ERROR     0
#define INCONSISTENT_K  1
#define OUT_OF_RANGE_K  2

/*
 * this is really defined in ffp.c but moved
 * out here to support custom calling sequences
 */
#define ALLCOEFF        21      /* length of c array */

#endif  /* FF_H */

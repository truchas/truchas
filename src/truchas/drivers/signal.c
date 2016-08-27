/* signal.C */

/*==============================================================================

  Copyright (c) Los Alamos National Security, LLC.  This file is part of the
  Truchas code (LA-CC-15-097) and is subject to the revised BSD license terms
  in the LICENSE file found in the top-level directory of this distribution.

==============================================================================*/

/*------------------------------------------------------------------------------
 * catch BSD style signals, and provide inquire routine
 * (to be called once/cycle, or when necessary) to see if a signal has
 * been caught.
 *
 * see src/drivers/driver.F90 for how this gets used
 *
 * signals caught:
 *
 *    SIGHUP
 *    SIGUSR2
 *    SIGURG
 *
 * Do _NOT_ catch SIGUSR1.  The mpich ch_p4 device uses SIGUSR1 to open
 * communications sockets.  The resulting hang is hard to diagnose.  And
 * it's likely that we will use mpich and ch_p4 sometime.
 *----------------------------------------------------------------------------*/

#include <signal.h>

static int HUP  = 0;
void sighup(int s)
{
  signal (SIGHUP,  &sighup);
  HUP++;
}

static int USR2 = 0;
void sigusr2(int s)
{
  signal (SIGUSR2, &sigusr2);
  USR2++;
}

static int URG = 0;
void sigurg(int s)
{
  signal (SIGURG, &sigurg);
  URG++;
}

void signalset()
{
  signal (SIGHUP,  &sighup );
  signal (SIGUSR2, &sigusr2);
  signal (SIGURG,  &sigurg);
}

void signalinquire(int *HUP_f90, int *USR2_f90, int *URG_f90)
{
  /* return any caught signals and clear the state */
  *HUP_f90  = HUP;  HUP  = 0;
  *USR2_f90 = USR2; USR2 = 0;
  *URG_f90  = URG;  URG  = 0;
}

/* signal.C end */

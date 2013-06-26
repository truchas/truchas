#!/usr/bin/env python
if __name__=='__main__':

    import os, sys

    #developers : please specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '/../..')

    errstring  = "\n\n    truchasdir is set to %s. No tools directory exists in this directory \n" %(truchasdir)
    errstring += "           Please alter your choice of truchasdir in your CapabilityTest script. \n"

    assert 'tools' in os.listdir(truchasdir), errstring

    extension   = '/tools/PythonPackages/TestingInfrastructure/'
    if os.path.isabs(truchasdir):
        testdir = os.path.abspath(truchasdir + extension)
    else:
        testdir = os.path.expanduser(truchasdir + extension)
    sys.path.append(testdir)
    sys.path.append(testdir+'/Runners')

    from RunThisCapability import RunThisCapability

    #developers : can specify clean=1 to remove all outputs and logdirs
    RunThisCapability(testdir,clean=0)

from TestCases import TruchasCapabilityTest

class CapabilityTest(TruchasCapabilityTest):
  "EM1 CapabilityTest"

  def setDataStores(self):
    "defines testdata and goldendata output directories needed by this Capability TestCase"

    self.testdirs   = ['em1_output']
    self.goldendirs = ['em1_golden']

  def setDefinitions(self):
    "Define the free-space and conducting regions"
    self.free_space = self.testdata[0].getRegion(['MB',10])
    self.conducting = self.testdata[0].getRegion(['MB',11,12])

  def testJouleHeat(self):
    "Verifying the Joule heat fields"

    tol = 1.0e-8  # Tolerance for the max relative error norm.
    tol2 = 1.0e-15  # Tight tolerance for scaled Joule heat check

    self.logger.write("Verifying the Joule heat fields...\n")

    fail = 0

    # Initial Joule heat calculation
    n = 0
    q = self.testdata[0].getField(field='Joule_P', cycle=n, region=self.conducting)
    qref = self.goldendata[0].getField(field='Joule_P', cycle=n, region=self.conducting)
    error = max(abs(q.data-qref.data)/qref.data)
    self.logger.write("  Cycle %2d: maximum relative error = %10.4e" %(n,error))
    if error > tol:
      fail = fail + 1
      self.logger.write(": FAIL (tolerance = %3.2e)\n" %(tol))
    else:
      self.logger.write(": PASS (tolerance = %8.2e)\n" %(tol))
    for j in range(1,6):
      qnext = self.testdata[0].getField(field='Joule_P', cycle=j, region=self.conducting)
      error = max(abs(qnext.data-q.data))
      if error != 0:
        fail = fail + 1
        self.logger.write("  Unexpected change in the Joule heat for cycle %2d: FAIL\n" %(j))

    # Scaled Joule heat
    n = 6
    qref.data = 4*q.data
    q = self.testdata[0].getField(field='Joule_P', cycle=n, region=self.conducting)
    error = max(abs(q.data-qref.data)/qref.data)
    self.logger.write("  Cycle %2d: maximum relative error = %10.4e" %(n,error))
    if error > tol2:
      fail = fail + 1
      self.logger.write(": FAIL (tolerance = %3.2e)\n" %(tol2))
    else:
      self.logger.write(": PASS (tolerance = %8.2e)\n" %(tol2))
    for j in range(7,11):
      qnext = self.testdata[0].getField(field='Joule_P', cycle=j, region=self.conducting)
      error = max(abs(qnext.data-q.data))
      if error != 0:
        fail = fail + 1
        self.logger.write("  Unexpected change in the Joule heat for cycle %2d: FAIL\n" %(j))

    # Zero Joule heat
    for j in range(11,16):
      q = self.testdata[0].getField(field='Joule_P', cycle=j)
      error = max(abs(q.data))
      if error != 0:
        fail = fail + 1
        self.logger.write("  Unexpected non-zero Joule heat for cycle %2d: FAIL\n" %(j))

    # Final Joule heat calculation
    n = 16
    q = self.testdata[0].getField(field='Joule_P', cycle=n, region=self.conducting)
    qref = self.goldendata[0].getField(field='Joule_P', cycle=n, region=self.conducting)
    error = max(abs(q.data-qref.data)/qref.data)
    self.logger.write("  Cycle %2d: maximum relative error = %10.4e" %(n,error))
    tol = 1.0e-8
    if error > tol:
      fail = fail + 1
      self.logger.write(": FAIL (tolerance = %3.2e)\n" %(tol))
    else:
      self.logger.write(": PASS (tolerance = %8.2e)\n" %(tol))
    for j in range(17,21):
      qnext = self.testdata[0].getField(field='Joule_P', cycle=j, region=self.conducting)
      error = max(abs(qnext.data-q.data))
      if error != 0:
        fail = fail + 1
        self.logger.write("  Unexpected change in the Joule heat for cycle %2d: FAIL\n" %(j))

    # Verify no Joule heat in the free-space region.
    for j in range(0,21):
      q = self.testdata[0].getField(field='Joule_P', cycle=j, region=self.free_space)
      error = max(abs(q.data))
      if error != 0:
        fail = fail + 1
        self.logger.write("  Unexpected non-zero free-space Joule heat for cycle %2d: FAIL\n" %(j))

    if fail: self.fail("Incorrect Joule heat fields: %2d failed tests" %(fail))

  def testFinalTemperature(self):
    "Verifying the final temperature field"

    tol = 1.0e-8  # Tolerance for the max relative error norm.

    self.logger.write("Verifying the final temperature field...\n")
    t = self.testdata[0].getField(field='T', cycle=20, region=self.conducting)
    tref = self.goldendata[0].getField(field='T', cycle=20, region=self.conducting)
    error = max(abs(t.data-tref.data)/tref.data)
    self.logger.write("  Maximum relative error = %10.4e" %(error))
    if error > tol:
      self.logger.write(": FAIL (tolerance = %3.2e)\n" %(tol))
      self.fail("  Error exceeds the tolerance %10.4e" %(tol))
    else:
      self.logger.write(": PASS (tolerance = %8.2e)\n" %(tol))


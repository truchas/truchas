#!/usr/bin/env python
if __name__=='__main__':

    import os, sys

    #developers : please specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '/../../')

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
  "DS9CapabilityTest"

  def setDataStores(self):
    "defines testdata and goldendata output directories needed by this Capability TestCase"

    self.testdirs   = ['ds9_output']
    #self.goldendirs = ['ds9_golden', 'ds9_pgolden']

#  def setDefinitions(self):
#    self.reg = self.testdata[0].getRegion(['Default','all'])

  def testTemp(self):
    "Verifying the temperature field"

    self.logger.write("Verifying the final temperature field ...\n")
    
    coord = self.testdata[0].getField(field='CENTROIDS', cycle=0)
    xcoord = coord.data[:,0]
    ycoord = coord.data[:,1]
    Tref = 9 + 6*xcoord*ycoord - xcoord*xcoord - ycoord*ycoord
    T = self.testdata[0].getField(field='T', cycle=67)
    error = max(abs(T.data-Tref)/Tref)

    tol = 5.0e-4
    self.logger.write("Maximum relative error = %10.4e" %(error))
    if error > tol:
        self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
        self.fail("  Error exceeds the tolerance")
    else:
        self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))

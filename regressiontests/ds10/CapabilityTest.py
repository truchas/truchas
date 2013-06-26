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

    self.testdirs   = ['ds10_output']
    #self.goldendirs = ['ds9_golden', 'ds9_pgolden']

#  def setDefinitions(self):
#    self.reg = self.testdata[0].getRegion(['Default','all'])

  def testTemp(self):
    "Verifying the temperature field"

    self.logger.write("Verifying the final temperature field ...\n")
    
    coord = self.testdata[0].getField(field='CENTROIDS', cycle=0)
    xcoord = coord.data[:,0]
    ycoord = coord.data[:,1]
    
    Tref = xcoord + ycoord
    for j in range(Tref.shape[0]):
      x = xcoord[j]
      y = ycoord[j]
      if x < 0:
        if y < 0:
          Tref[j] += 2.5
        else:
          Tref[j] += 2.7
      else:
        if y < 0:
          Tref[j] += 3.3
        else:
          Tref[j] += 3.5

    T = self.testdata[0].getField(field='T', cycle=233)
    
    # Take care of gap cells so they don't enter into error calc
    for j in range(Tref.shape[0]):
      x = xcoord[j]
      y = ycoord[j]
      if abs(x) < 1.0e-3 or abs(y) < 1.0e-3:
        Tref[j] = T.data[j]
        
    error = max(abs(T.data-Tref)/Tref)

    tol = 1.0e-7
    self.logger.write("Maximum relative error = %10.4e" %(error))
    if error > tol:
        self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
        self.fail("  Error exceeds the tolerance")
    else:
        self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))

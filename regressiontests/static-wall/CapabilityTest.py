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
  "static-wall CapabilityTest"

  def setDataStores(self):
    "defines testdata and goldendata output directories needed by this Capability TestCase"

    self.testdirs   = ['static-wall_output']
    
  def testVelocity(self):
    "Verifying maximum velocity"
    
    from math import sqrt
    
    tol = 1.0e-11
    
    V = self.testdata[0].getField(field='Velocity',cycle=1)
    
    a = V.data[:,0] # making an array of the correct type and size
    for j in range(a.shape[0]):
      a[j] = sqrt(V.data[j,0]**2 + V.data[j,1]**2 + V.data[j,2]**2)
    error = max(a)
    if error > tol:
        self.fail (" Maximum velocity exceeds the tolerance")


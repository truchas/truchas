#!/usr/bin/env python
try:
    import numpy.oldnumeric as Numeric
except ImportError:
    import Numeric
except:
    raise


def AnalyticSolution(Ncellx,Ncelly):

    U = Numeric.zeros((Ncellx,Ncelly), 'd')

    for j in range(0,Ncelly):
        y = (j + 0.5)/Ncelly  # the y-coordinate of the cell center
        if y < 0.5:
            U[0][j]=100.0/500.5*y
        else:
            U[0][j]=100000.0/500.5*y - 49950.0/500.5

    for i in range(1,Ncellx):
        for j in range(0,Ncelly):
            U[i][j] = U[0][j]
        
    return U


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
  "HT CapabilityTest"

  def setDataStores(self):
    "defines testdata  directory needed by this Capability TestCase"

    self.testdirs   = ['conduction_90_output']

  def setTolerances(self):
    "define tolerances"

    self.tol['temperature l infinity error'] = 1e-8


  def setDefinitions(self):
    "Define the regions"
    self.reg = self.testdata[0].getRegion(['Default'])

  def testHTcond90(self):
    "Verifying Heat Transfer with ortho and a vertical conductivity jump"

    self.logger.write("Verifying Heat Transfer with ortho and a vertical conductivity Jump... ")

    Ncellx = 10
    Ncelly = 10
    
    temp_analytic = Numeric.ravel(AnalyticSolution(Ncellx,Ncelly)) # flatten the anaytical solution array

    temp = self.testdata[0].getField(field='T',
                                      cycle=11,
                                      region=self.reg)

    TLinferr = self.measure.linfError(temp.data,temp_analytic)
    
    if TLinferr > self.tol['temperature l infinity error']:
        self.logger.write(": FAIL (tolerance = %12.7e)\n" %(TLinferr))
        self.fail("The temperature does not match the analytic solution.")
    else:
        self.logger.write(": PASS (tolerance = %12.7e)\n" %(TLinferr))
       

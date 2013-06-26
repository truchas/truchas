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
  "Grain Growth RNT CapabilityTest"

  def setTolerances(self):
    "Defines error tolerances"

    self.tol['R_error']  = 1.0e-04
    self.tol['UC_error'] = 1.0e-01
    
  def setDataStores(self):
    "Defines testdata output directories needed by this Capability TestCase"

    self.testdirs   = ['grain_growth_rnt_output']

  def setDefinitions(self):
    "Define the radii and undercooling at cell 40 (for t=15.0)" 

    self.radiusc40    = 1.137e-09
    self.undercoolc40 = 0.261

    self.c40          = self.testdata[0].getRegion(selector = ['ID','CELL','40'],
                                                   desc     = 'Cell 40')
    
  def testRadiusAndUC(self):
    "verifying grain growth radius and undercooling at cell 40"

    #get list of timesteps from Truchas testdata output that satisfy 1.5<=t<=3.5 
    these_timesteps = self.testdata[0].getTimeSteps(timerange = [15,16],
                                                    dtstep    = 1.0)

    for timestep in these_timesteps:

        thistime   = timestep.time
        index      = these_timesteps.index(timestep)

        if index == 0:
            rc40        = self.testdata[0].getField(field    = 'Grain_R',
                                                    time     = thistime,
                                                    region   = self.c40)

            print 'XXX..'
            print rc40.data
            err         = self.measure.percentError(rc40.data[0],self.radiusc40)
            #fail if err > R_error
            failmessage = 'In %s: R percent error = %10.3e is > %10.3e' %(self.methodName, err, self.tol['R_error'])
	    self.failIf(err>self.tol['R_error'],msg = failmessage)

            ucc40        = self.testdata[0].getField(field    = 'Max_Underc',
                                                     time     = thistime,
                                                     region   = self.c40)
            
            err        = self.measure.percentError(ucc40.data[0],self.undercoolc40)
            #fail if err > UC_error
            failmessage = 'In %s: UC percent error = %10.3e is > %10.3e' %(self.methodName, err, self.tol['UC_error'])
	    self.failIf(err>self.tol['UC_error'],msg = failmessage)

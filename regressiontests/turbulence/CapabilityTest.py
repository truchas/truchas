if __name__=='__main__':

    import os, sys
    # specify the location of the truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '../../../')

    errstring  = "\n\n    truchasdir is set to %s.  \
                  No tools directory exists in this directory \n" %(truchasdir)
    errstring += \
       "           Please alter truchasdir in your CapabilityTest script. \n"
    
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
    "Turbulence test"

    def setTolerances(self):
        "define tolerances here"
        
        # Setup pressure tolerance a little tighter since we're really
        # just checking a linear profile
        self.tol['PL1_error'] = 1.0e-09
        self.tol['UL1_error'] = 1.2e-06

    def setDataStores(self):
        "defines all datastores to be used throughout this TestCase"

        #define testdata output directories
        self.testdirs = ['turbulence_output']

    def setDefinitions(self):
        "define analytic velocity at probes"

        # Empirically determined velocities :( MAC
        # The pressure values are from the linear pressure profile
        # The min/max velocity values are essentially "golden" values
	self.umin =  632.876
	self.umax = 2062.680
        self.pmin = 25000.0
	self.pmax = 75000.0

    def testMinMaxProbes(self):
        "Tests velocity & pressure at probes 'umin/umax'"

        umin = self.testdata[0].getProbe(name='umin',
                                         field='Velocity',
                                         timerange=[0.250,0.251]) 

        umax = self.testdata[0].getProbe(name='umax',
                                         field='Velocity',
                                         timerange=[0.250,0.251]) 

        pmax = self.testdata[0].getProbe(name='umax',
                                         field='P',
                                         timerange=[0.250,0.251]) 

	pmin = self.testdata[0].getProbe(name='pdown',
                                         field='P',
                                         timerange=[0.250,0.251]) 


        #velocity probe data is always stored [cycle,time,(u,v,w)] 
        umin_chk = umin.data[0,2]
        umax_chk = umax.data[0,2]
        #p probe data is always stored [cycle,time,p] 
        pmin_chk = pmin.data[0,2]
        pmax_chk = pmax.data[0,2]

        #get % error = |umax - self.umax|/|self.umax|
 	umin_err = self.measure.percentError(umin_chk,self.umin)
 	umax_err = self.measure.percentError(umax_chk,self.umax)
 	pmin_err = self.measure.percentError(pmin_chk,self.pmin)
 	pmax_err = self.measure.percentError(pmax_chk,self.pmax)

        logmessage = '\n Min. Velocity = %g\n' %(umin_chk)
        self.logger.write(logmessage)
        logmessage = '\n Max. Velocity = %g\n' %(umax_chk)
        self.logger.write(logmessage)
        logmessage = '\n Min. Pressure = %g\n' %(pmin_chk)
        self.logger.write(logmessage)
        logmessage = '\n Max. Pressure = %g\n' %(pmax_chk)
        self.logger.write(logmessage)

        # Failure criteria 
        # umin
        failmessage  = 'In %s: umin percent_error = %10.3e is > %10.3e' \
                        %(self.methodName, umin_err, self.tol['UL1_error'])
	self.failIf(umin_err > self.tol['UL1_error'], msg = failmessage)

 	# umax
        failmessage  = 'In %s: umax percent_error = %10.3e is > %10.3e' \
	                %(self.methodName, umax_err, self.tol['UL1_error'])
	self.failIf(umax_err > self.tol['UL1_error'], msg = failmessage)

	# pmin
        failmessage  = 'In %s: pmin percent_error = %10.3e is > %10.3e' \
	                %(self.methodName, pmin_err, self.tol['PL1_error'])
	self.failIf(pmin_err > self.tol['PL1_error'], msg = failmessage)

	# pmax
        failmessage  = 'In %s: pmax percent_error = %10.3e is > %10.3e' \
	                %(self.methodName, pmax_err, self.tol['PL1_error'])
	self.failIf(pmax_err > self.tol['PL1_error'], msg = failmessage)

if __name__=='__main__':

    import os, sys
    #developers : please specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '../../../')

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

    def setTolerances(self):
        "define tolerances here"
        
        self.tol['V_error']= 0.015

    def setDataStores(self):
        "defines all datastores to be used throughout this TestCase"

        #define testdata output directories..these result from running a BasicRun TestCase
        self.testdirs  = ['devahl_davis_restart_output']

    def setDefinitions(self):
        "define analytic velocities at probes"

        self.uanalytic = 7.585e-5
        self.wanalytic = 7.685e-5

    def testUVelocityProbe(self):
        "tests u velocity at probe 'VhMax'"

        VhMax        = self.testdata[0].getProbe(name='VhMax',field='Velocity',timerange=[60499,61010]) 
        #VhMax.str()

        #velocity probe data is always stored [cycle,time,(u,v,w)] 
        U_VhMax      = VhMax.data[0,2]

        #get % error = |U_VhMax - self.uanalytic|/|self.uanalytic|
	err          = self.measure.percentError(U_VhMax,self.uanalytic)

        logmessage   = '\n Percentage u error = %10.3e \n' %(err)
        self.logger.write(logmessage)

        #fail if err < self.tol['V_error']
        failmessage  = 'In %s: U percent_error = %10.3e is > %10.3e' %(self.methodName, err, self.tol['V_error'])

	self.failIf(err > self.tol['V_error'],
                    msg = failmessage)

    def testWVelocityProbe(self):
        "tests w velocity at probe 'VvMax'"

        VvMax        = self.testdata[0].getProbe(name='VvMax',field='Velocity',timerange=[60499,61010])
        #VvMax.str()

        #velocity probe data is always stored [cycle,time,(u,v,w)] 
        W_VvMax      = VvMax.data[0,4]

        #get % error = |W_VvMax - self.wanalytic|/|self.wanalytic|
	err          = self.measure.percentError(W_VvMax,self.wanalytic)

        logmessage   = '\n Percentage w error = %10.3e \n' %(err)
        self.logger.write(logmessage)

        #fail if err < self.tol['V_error']
        failmessage  = 'In %s: W percent_error = %10.3e is > %10.3e' %(self.methodName, err, self.tol['V_error'])

	self.failIf(err > self.tol['V_error'], msg = failmessage)

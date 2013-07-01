if __name__=='__main__':

    import os, sys
    #developers : please specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '../../../../../../')

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
        
        self.tol['P_error'] = 10e-1

    def setDataStores(self):
        "defines all datastores to be used throughout this TestCase"

        #define testdata output directories..these result from running a BasicRun TestCase
        self.testdirs     = ['1D_diverging_duct_output']

    def setDefinitions(self):
        "defines analytic expressions and regions needed by this Capability TestCase"

        X              = self.testdata[0].getField(field   = 'CENTROIDS',
                                                   cycle   = 0)
        
        #define analytic pressure..
        self.Panalytic = 0.5 - 0.5*(1.0 + 0.1*X.data[:,0])**(-2.0)

        #define 1D analytic velocity
        self.Vanalytic = (1.0 + 0.1*X.data[:,0])**(-1.0)


    def testPressure(self):
	"tests pressure "

        #get pressure field 
	P  = self.testdata[0].getField(field  = 'P',
                                       cycle  = 146)
        #P.str()
        
	#log P field
	self.logger.writeField(P, desc = 'Pressure in the entire domain')

	#get linf_error(P - Panalytic) 
	err          = self.measure.linfError(P.data,self.Panalytic)

        logmessage   = '\n linf P error = %10.3e \n' %(err)
        self.logger.write(logmessage)
	
	#fail if err < P_error
        failmessage = 'In %s: Pressure drop percent_error = %10.3e is > %10.3e' %(
                       self.methodName, err, self.tol['P_error'])

	self.failIf(err > self.tol['P_error'],
                    msg = failmessage)
 





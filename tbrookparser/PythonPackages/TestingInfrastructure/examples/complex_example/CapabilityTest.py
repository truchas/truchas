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
        
        self.tol['PD_error'] = 10e-01

    def setDataStores(self):
        "defines analytic results and testdata output directories needed by this Capability TestCase"

        #define analytic pressure drop
	self.PD_analytic  = 100.0 

        #define testdata output directories..these result from running a BasicRun TestCase
        self.testdirs     = ['complex_example_output']

        #define restartdata output directories...these result from running a RestartRun TestCase
        self.restartdirs  = ['complex_example_restart_output']

    def setDefinitions(self):
        "developer defined regions and analytic results to be used throughout this TestCase"

        #get region representing the inside of the static drop 
        self.rinside     = self.testdata[0].getRegion(selector=['VAR','VOF0001',[0.99,1.0]],
                                                      desc    = 'Region inside the drop')

        #get region representing the outside of the static drop 
        self.routside    = self.testdata[0].getRegion(selector=['VAR','VOF0001',[0.0,0.01]],
                                                     desc     = 'Region outside the drop')

        #get region representing the entire mesh from the restartdata store
        self.meshregion  = self.restartdata[0].getRegion(selector=['Default'],
                                                         desc    = 'Region representing the entire Truchas mesh')


    def testPressureDrop(self):
	"tests interface pressure drop"

        #get pressure field inside the drop at cycle 1
	Pin  = self.testdata[0].getField(field  = 'P',
                                         cycle  = 1,
                                         region = self.rinside)

	#get pressure field outside the drop at cycle 1
	Pout = self.testdata[0].getField(field  = 'P',
                                         cycle  = 1,
                                         region = self.routside)

	#log Pin and Pout fields
	self.logger.writeField(Pin,
                               region = self.rinside,
                               desc   = 'Total Measure PD : Pressure inside the drop')
        
	self.logger.writeField(Pout,
                               region = self.routside,
                               desc   = 'Total Measure PD : Pressure outside the drop')

	#get mean pressure inside the drop..
        #note all self.measure methods take a Numeric array as input (Pin.data)
	Pin_mean     = self.measure.meanValue(Pin.data)
        
	#get mean pressure outside the drop..
        #note all self.measure methods take a Numeric array as input (Pout.data)
	Pout_mean    = self.measure.meanValue(Pout.data)

	#get mean Pressure drop
	PD_mean      = Pin_mean - Pout_mean

	#get % error = |PD_mean - PD_analytic|/|PD_analytic|
	err          = self.measure.percentError(PD_mean,self.PD_analytic)

        logmessage   = '\n Percentatge PD error = %10.3e \n' %(err)
        self.logger.write(logmessage)
	
	#fail if err < PD_error
        failmessage  = 'In %s: Pressure drop percent_error = %10.3e is > %10.3e' %(self.methodName, err, self.tol['PD_error'])

	self.failIf(err > self.tol['PD_error'],
                    msg = failmessage)
 




if __name__=='__main__':

    import os, sys
    #specify the location of the truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '../../../')

    errstring  = "\n\n    truchasdir is set to %s. \
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
import os, sys

class CapabilityTest(TruchasCapabilityTest):

    def setTolerances(self):
        "define tolerances here"
        
        self.tol['eps_minmax'] = 2.0e-2
        self.tol['eps_avg']    = 2.0e-2
        self.tol['VL1_error'] = 1.0e-03

    def setDataStores(self):
        "defines all datastores to be used throughout this TestCase"

        # Golden-output is not required for this test case
        self.testdirs     = ['static_drop_output']

    def setDefinitions(self):
        "defines analytic expressions and regions needed by this Capability TestCase"

        # define analytic pressure drop based on surface tension and radius
	self.PD_analytic  = 73.0

        # Get region representing the inside of the static drop 
        # Tune inside/outside based on volume-fraction values
        Fmin = 0.99
        Fmax = 1.00
        self.rinside = \
            self.testdata[0].getRegion(selector=['VAR','VOF0001',[Fmin,Fmax]],
                                       desc    = 'Region inside the drop')

        # Get region representing the outside of the static drop 
        self.routside = \
            self.testdata[0].getRegion(selector=['VAR','VOF0002',[Fmin,Fmax]],
                                       desc     = 'Region outside the drop')

        #get region representing the entire mesh
        self.meshregion  = self.testdata[0].getRegion(selector=['Default'],
                                                      desc    = 'The entire mesh')


    def testPressureDrop(self):
	"tests interface pressure drop"

        # Case-1: Pressure drop based on min/max pressures in the domain
	P = self.testdata[0].getField(field  = 'P',
                                      cycle  = 1,
                                      region = self.meshregion)

	Pmin = self.measure.minValue(P.data)
	Pmax = self.measure.maxValue(P.data)

	err_minmax = self.measure.percentError((Pmax - Pmin),self.PD_analytic)
        logmsg = '\n Min/Max Percentage error in pressure drop = %10.3e \n' \
                 % (err_minmax)
        self.logger.write(logmsg)


        # Case-2: Use the inside/outside regions to get average pressures
        # get pressure field inside the drop at cycle 1
	Pin  = self.testdata[0].getField(field  = 'P',
                                         cycle  = 1,
                                         region = self.rinside)


	# get pressure field outside the drop at cycle 1
	Pout = self.testdata[0].getField(field  = 'P',
                                         cycle  = 1,
                                         region = self.routside)

	# log Pin and Pout fields
	# self.logger.writeField(Pin,
        #                      region = self.rinside,
        #                      desc   = 'Total Measure PD : Pressure inside the drop')
        
	# self.logger.writeField(Pout,
        #                      region = self.routside,
        #                      desc   = 'Total Measure PD : Pressure outside the drop')

	# Get mean pressure inside and outside the droplet
        # Note all self.measure methods take a Numeric array as input (Pin.data)

	Pin_avg  = self.measure.meanValue(Pin.data)

	Pout_avg = self.measure.meanValue(Pout.data)

	# Get pressure drop based on the mean pressures
	PD_avg = Pin_avg - Pout_avg

	#get % error = |PD_avg - PD_analytic|/|PD_analytic|
	err_avg = self.measure.percentError(PD_avg,self.PD_analytic)

        logmsg = '\n Percentage PD error = %10.3e \n' % (err_avg)
        self.logger.write(logmsg)
	
	# fail if err_minmax < eps_minmax
        failmessage = 'In %s: Min/Max Pressure drop error = %10.3e is > %10.3e' \
                      %(self.methodName, err_minmax, self.tol['eps_minmax'])
	self.failIf(err_minmax > self.tol['eps_minmax'], msg = failmessage)

	# fail if avg_err < PD_error
        failmessage = 'In %s: Pressure drop percent_error = %10.3e is > %10.3e' \
                      %(self.methodName, err_avg, self.tol['eps_avg'])
	self.failIf(err_avg > self.tol['eps_avg'], msg = failmessage)

 
    def testSpuriousVelocities(self):
	"Test for spurious velocity generation in the entire domain"

 	#get velocities on the entire domain
	V     = self.testdata[0].getField(field  = 'Velocity',
                                          region = self.meshregion,
                                          cycle  = 1)

	#Ucomp = V.data[:,0]
	#Vcomp = V.data[:,1]
        #V.str()

        #log V field
	self.logger.writeField(V,
                               region = self.meshregion,
                               desc   = 'Velocity at cycle 1')

        # Umin = self.measure.minValue(V.data[:,0])
	# Umax = self.measure.maxValue(V.data[:,0])
        # print "\nUmin,Umax = %10.4e, %10.4e\n" % (Umin, Umax)

	#get L1 error
	L1err    = self.measure.l1Error(V.data)

        #fail if L1err > VL1_error
        failmessage = 'In %s velocity L1_error = %10.3e is > %10.3e' \
                      %(self.methodName, L1err, self.tol['VL1_error'])

	self.failIf(L1err > self.tol['VL1_error'], msg = failmessage)

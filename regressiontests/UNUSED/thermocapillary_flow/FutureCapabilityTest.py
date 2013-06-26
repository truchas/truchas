if __name__=='__main__':

    import os, sys
    #specify the location of your truchas checkout directory 'truchasdir'
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

class CapabilityTest(TruchasCapabilityTest):

    def setTolerances(self):
        "define tolerances here"
        
        self.tol['UL1_error']   = 1.0e-5
        self.tol['VL1_error']   = 1.0e-5
        self.tol['WL1_error']   = 1.0e-12  # epsilon above 0.0
        self.tol['PL1_error']   = 1.0e-5
        self.tol['VOFL1_error'] = 1.0e-5

    def setDataStores(self):
        "defines all datastores to be used throughout this TestCase"

        #define testdata output directories
        self.testdirs  = ['broken_dam_hex_restart_output']

    def setDefinitions(self):
        "Define known values at probe"

        # velocities, pressure and volume-fraction
        self.velx = 1.15013
        self.vely = 3.33312e-2
        self.velz = 0.00000     # This is really unused
        self.p    = 5.83694
        self.vof  = 2.20233e-1
        
    def testProbeData(self):
        "Test velocity, pressure, volume-fraction at probe 'dbp'"

        dbp_vel = self.testdata[0].getProbe(name='dbp',
                                            field='Velocity',
                                            timerange=[0.080,0.081]) 

        dbp_p   = self.testdata[0].getProbe(name='dbp',
                                            field='P',
                                            timerange=[0.080,0.081]) 

        dbp_vof = self.testdata[0].getProbe(name='dbp',
                                            field='VOF0001',
                                            timerange=[0.080,0.081]) 

        velx_chk = dbp_vel.data[0,2]
        vely_chk = dbp_vel.data[0,3]
        velz_chk = dbp_vel.data[0,4]
        vof_chk  = dbp_vof.data[0,2]
        p_chk    = dbp_p.data[0,2]

        #get % error = |umax - self.umax|/|self.umax|
 	velx_err = self.measure.percentError(velx_chk,self.velx)
 	vely_err = self.measure.percentError(vely_chk,self.vely)
 	p_err    = self.measure.percentError(p_chk,   self.p   )
 	vof_err  = self.measure.percentError(vof_chk, self.vof )

        # Should be essentially zero
 	velz_err = abs(velz_chk)

        logmessage = '\n X-Velocity      = %g\n' %(velx_chk)
        self.logger.write(logmessage)
        logmessage = '\n Y-Velocity      = %g\n' %(vely_chk)
        self.logger.write(logmessage)
        logmessage = '\n Z-Velocity      = %g\n' %(velz_chk)
        self.logger.write(logmessage)
        logmessage = '\n Pressure        = %g\n' %(p_chk)
        self.logger.write(logmessage)
        logmessage = '\n Volume Fraction = %g\n' %(vof_chk)
        self.logger.write(logmessage)
    
        # Failure criteria 
        # velx
        failmessage  = 'In %s: velx percent_error = %10.3e is > %10.3e' \
                        %(self.methodName, velx_err, self.tol['UL1_error'])
	self.failIf(velx_err > self.tol['UL1_error'], msg = failmessage)

        # vely
        failmessage  = 'In %s: vely percent_error = %10.3e is > %10.3e' \
                        %(self.methodName, vely_err, self.tol['VL1_error'])
	self.failIf(vely_err > self.tol['VL1_error'], msg = failmessage)

        # velz
        failmessage  = 'In %s: velz = %10.3e is > %10.3e' \
                       %(self.methodName, velz_err, self.tol['WL1_error'])
        self.failIf(velz_err > self.tol['WL1_error'], msg = failmessage)

        # pressure
        failmessage  = 'In %s: pressure percent_error = %10.3e is > %10.3e' \
                        %(self.methodName, p_err, self.tol['PL1_error'])
	self.failIf(p_err > self.tol['PL1_error'], msg = failmessage)

        # vof
        failmessage  = 'In %s: pressure percent_error = %10.3e is > %10.3e' \
                        %(self.methodName, vof_err, self.tol['VOFL1_error'])
	self.failIf(vof_err > self.tol['VOFL1_error'], msg = failmessage)

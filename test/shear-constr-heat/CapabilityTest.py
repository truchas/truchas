try:
    import numpy.oldnumeric as Numeric
except ImportError:
    import Numeric
except:
    raise

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
    "Shear constrain heat capability test"
    
    def setTolerances(self):
        "defines all tolerances to be used in this testcase"

        self.tol['Stress_rel']    = 1.0e-8
        self.tol['Stress_abs']    = 1.0e-8
        self.tol['Strain_rel']    = 1.0e-8
        self.tol['Strain_abs']    = 1.0e-8
        self.tol['Displacement']  = 1.0e-10

    def setDataStores(self):
        "defines testdata and goldendata output directories needed by this Capability TestCase"

        self.testdirs   = ['shear_constr_heat_unstruc_output']
        self.goldendirs = ['shear_constr_heat_unstruc_golden']

    def setDefinitions(self):
        "Define the data fields"

        self.whole_mesh  = self.testdata[0].getRegion(['Default'])
        self.sigma_ss    = self.testdata[0].getField(field='sigma', cycle=0, region= self.whole_mesh)
        self.epsilon_ss  = self.testdata[0].getField(field='epsilon', cycle=0, region= self.whole_mesh)
        self.coords      = self.testdata[0].getField(field='NODES', cycle=0, region= self.whole_mesh)
        self.disp        = self.testdata[0].getField(field='disp', cycle=0, region= self.whole_mesh)


    def testSteadyStateStress(self):
        "Compare stresses at initial temperature"

        n = 0
        sigzzref = -152.533333333
        sigxzref = 500.00

        XXcomp = self.sigma_ss.data[:,0]
        YYcomp = self.sigma_ss.data[:,1]
        ZZcomp = self.sigma_ss.data[:,2]
        XYcomp = self.sigma_ss.data[:,3]
        XZcomp = self.sigma_ss.data[:,4]
        YZcomp = self.sigma_ss.data[:,5]

        errorxx = max(abs(XXcomp))
        erroryy = max(abs(YYcomp))
        errorzz = max(abs((ZZcomp - sigzzref) / sigzzref))
        errorxy = max(abs(XYcomp))
        errorxz = max(abs((XZcomp - sigxzref) / sigxzref))
        erroryz = max(abs(YZcomp))

        self.logger.write("  Cycle %2d: maximum absolute sigxx error = %10.4e\n" %(n,errorxx))    
        self.logger.write("  Cycle %2d: maximum absolute sigyy error = %10.4e\n" %(n,erroryy))    
        self.logger.write("  Cycle %2d: maximum relative sigzz error = %10.4e\n" %(n,errorzz))
        self.logger.write("  Cycle %2d: maximum absolute sigxy error = %10.4e\n" %(n,errorxy))    
        self.logger.write("  Cycle %2d: maximum relative sigxz error = %10.4e\n" %(n,errorxz))
        self.logger.write("  Cycle %2d: maximum absolute sigyz error = %10.4e\n" %(n,erroryz))

        fail = 0
        if errorxx > self.tol['Stress_abs']:
            fail = fail + 1
            self.logger.write("  Steady state sigxx not within %82e: FAIL\n" %(self.tol['Stress_abs']))
        if erroryy > self.tol['Stress_abs']:
            fail = fail + 1
            self.logger.write("  Steady state sigyy not within %82e: FAIL\n" %(self.tol['Stress_abs']))
        if errorzz > self.tol['Stress_rel']:
            fail = fail + 1
            self.logger.write("  Steady state sigxx not within %82e: FAIL\n" %(self.tol['Stress_rel']))
        if errorxy > self.tol['Stress_abs']:
            fail = fail + 1
            self.logger.write("  Steady state sigxy not within %82e: FAIL\n" %(self.tol['Stress_abs']))
        if errorxz > self.tol['Stress_rel']:
            fail = fail + 1
            self.logger.write("  Steady state sigxz not within %82e: FAIL\n" %(self.tol['Stress_rel']))
        if erroryz > self.tol['Stress_abs']:
            fail = fail + 1
            self.logger.write("  Steady state sigyz not within %82e: FAIL\n" %(self.tol['Stress_abs']))

        if fail:
            self.fail("Incorrect stress fields: %2d failed tests" %(fail))
        else:
            self.logger.write(": PASS \n")


    def testSteadyStateStrain(self):
        "Compare strains at initial temperature"

        n = 0
        epsxxref = 2.93333333333e-3 
        epsyyref = 2.93333333333e-3 
        epsxzref = 500.0 / 5.2e4

        XXcomp = self.epsilon_ss.data[:,0]
        YYcomp = self.epsilon_ss.data[:,1]
        ZZcomp = self.epsilon_ss.data[:,2]
        XYcomp = self.epsilon_ss.data[:,3]
        XZcomp = self.epsilon_ss.data[:,4]
        YZcomp = self.epsilon_ss.data[:,5]

        errorxx = max(abs((XXcomp - epsxxref) / epsxxref))
        erroryy = max(abs((YYcomp - epsyyref) / epsyyref))
        errorzz = max(abs(ZZcomp))
        errorxy = max(abs(XYcomp))
        errorxz = max(abs((XZcomp - epsxzref) / epsxzref))
        erroryz = max(abs(YZcomp))

        self.logger.write("  Cycle %2d: maximum relative epsxx error = %10.4e\n" %(n,errorxx))    
        self.logger.write("  Cycle %2d: maximum relative epsyy error = %10.4e\n" %(n,erroryy))    
        self.logger.write("  Cycle %2d: maximum absolute epszz error = %10.4e\n" %(n,errorzz))
        self.logger.write("  Cycle %2d: maximum absolute epsxy error = %10.4e\n" %(n,errorxy))    
        self.logger.write("  Cycle %2d: maximum relative epsxz error = %10.4e\n" %(n,errorxz))
        self.logger.write("  Cycle %2d: maximum absolute epsyz error = %10.4e\n" %(n,erroryz))

        fail = 0
        if errorxx > self.tol['Strain_abs']:
            fail = fail + 1
            self.logger.write("  Steady state epsxx not within %82e: FAIL\n" %(self.tol['Strain_abs']))
        if erroryy > self.tol['Strain_abs']:
            fail = fail + 1
            self.logger.write("  Steady state epsyy not within %82e: FAIL\n" %(self.tol['Strain_abs']))
        if errorzz > self.tol['Strain_rel']:
            fail = fail + 1
            self.logger.write("  Steady state epsxx not within %82e: FAIL\n" %(self.tol['Strain_rel']))
        if errorxy > self.tol['Strain_abs']:
            fail = fail + 1
            self.logger.write("  Steady state epsxy not within %82e: FAIL\n" %(self.tol['Strain_abs']))
        if errorxz > self.tol['Strain_rel']:
            fail = fail + 1
            self.logger.write("  Steady state epsxz not within %82e: FAIL\n" %(self.tol['Strain_rel']))
        if erroryz > self.tol['Strain_abs']:
            fail = fail + 1
            self.logger.write("  Steady state epsyz not within %82e: FAIL\n" %(self.tol['Strain_abs']))

        if fail:
            self.fail("Incorrect strain fields: %2d failed tests" %(fail))
        else:
            self.logger.write(": PASS \n")

    def testSteadyStateDisplacement(self):
        "Compare displacements at initial temperature"

        n = 0
        nnodes = 78

        epsilonxx = 2.93333333333e-3
        epsilonyy = 2.93333333333e-3
        epsilonxz = 500.0 / 5.2e4

        error = 0.0
        disp_ana = [0.0, 0.0, 0.0]
        disp_ana = Numeric.array(disp_ana)

        for j in range(nnodes):
            disp_ana[0] = epsilonxx * self.coords.data[j,0] + epsilonxz * 2.0 * self.coords.data[j,2]
            disp_ana[1] = epsilonyy * self.coords.data[j,1]
            for k in range(3):
                diff = abs(self.disp.data[j,k] - disp_ana[k])
                if diff > error: error = diff

        self.logger.write("  Cycle %2d, Displacement component maximum absolute error = %10.4e\n" %(n,error))

        if error > self.tol['Displacement']:
            self.fail("Incorrect displacement field")
        else:
            self.logger.write(": PASS \n")


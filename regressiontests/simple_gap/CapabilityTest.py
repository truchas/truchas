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

        self.tol['Stress_rel']    = 1.0e-6
        self.tol['Strain_rel']    = 1.0e-7

    def setDataStores(self):
        "defines testdata and goldendata output directories needed by this Capability TestCase"

        self.testdirs   = ['simple_gap_output']
        self.goldendirs = ['simple_gap_golden']

    def setDefinitions(self):
        "Define the data fields"

        self.whole_mesh  = self.testdata[0].getRegion(['Default'])
        self.sigma       = self.testdata[0].getField(field='sigma', cycle=0, region= self.whole_mesh)
        self.epsilon     = self.testdata[0].getField(field='epsilon', cycle=0, region= self.whole_mesh)
        self.ntrac_7     = self.testdata[0].getField(field='NTRAC_07', cycle=0, region= self.whole_mesh)
        self.ntrac_8     = self.testdata[0].getField(field='NTRAC_08', cycle=0, region= self.whole_mesh)
        self.sigma_ref   = self.goldendata[0].getField(field='sigma', cycle=0, region= self.whole_mesh)
        self.epsilon_ref = self.goldendata[0].getField(field='epsilon', cycle=0, region= self.whole_mesh)
        self.ntrac_7_ref = self.goldendata[0].getField(field='NTRAC_07', cycle=0, region= self.whole_mesh)
        self.ntrac_8_ref = self.goldendata[0].getField(field='NTRAC_08', cycle=0, region= self.whole_mesh)
        
    def testStress(self):
         "Compare stresses"

         n = 0
         fail = 0

         for j in range(6):
             i = j + 1
             error = max(abs((self.sigma.data[:,j]-self.sigma_ref.data[:,j])/self.sigma_ref.data[:,j]))

             self.logger.write("  Cycle %2d, Sigma %2d: maximum relative error = %10.4e\n" %(n,j,error))
             if error > self.tol['Stress_rel']:
                 self.logger.write("  Unexpected change in the stress %2d for cycle %2d: FAIL\n" %(i,n))
                 fail = fail + 1

         if fail:
             self.fail("Incorrect stress fields: %2d failed tests\n" %(fail))
         else:
             self.logger.write(": PASS \n")

    def testStrain(self):
         "Compare strains"

         n = 0
         fail = 0

         for j in range(6):
             i = j + 1
             error = max(abs((self.epsilon.data[:,j]-self.epsilon_ref.data[:,j])/self.epsilon_ref.data[:,j]))
    
             self.logger.write("  Cycle %2d, Epsilon %2d: maximum relative error = %10.4e\n" %(n,j,error))
             if error > self.tol['Strain_rel']:
                 self.logger.write("  Incorrect change in strain %2d for cycle %2d: FAIL\n" %(i,n))
                 fail = fail + 1

         if fail:
             self.fail("Incorrect strain fields: %2d failed tests\n" %(fail))
         else:
             self.logger.write(": PASS \n")

    def testNorm_Trac_07(self):
         "Compare normal traction on sideset 7"

         n = 0
         nnodes = 150
         error = 0.0
         
         for j in range(nnodes):
             if self.ntrac_7_ref.data[j] != 0.0:
                 diff = abs((self.ntrac_7.data[j]-self.ntrac_7_ref.data[j])/self.ntrac_7_ref.data[j])
                 if diff > error: error = diff
                 
         self.logger.write("  Cycle %2d, Normal traction 7: maximum relative error = %10.4e\n" %(n,error))
         if error > self.tol['Stress_rel']:
             self.fail("Incorrect normal traction on sideset 7\n" %(fail))
         else:
             self.logger.write(": PASS \n")

    def testNorm_Trac_08(self):
         "Compare normal traction on sideset 8"

         n = 0
         nnodes = 150
         error = 0.0
         
         for j in range(nnodes):
             if self.ntrac_8_ref.data[j] != 0.0:
                 diff = abs((self.ntrac_8.data[j]-self.ntrac_8_ref.data[j])/self.ntrac_8_ref.data[j])
                 if diff > error: error = diff
                 
         self.logger.write("  Cycle %2d, Normal traction 8: maximum relative error = %10.4e\n" %(n,error))
         if error > self.tol['Stress_rel']:
             self.fail("Incorrect normal traction on sideset 8\n" %(fail))
         else:
             self.logger.write(": PASS \n")

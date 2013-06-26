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
    "Viscoplastic ring capability test"
    
    def setTolerances(self):
        "defines all tolerances to be used in this testcase"

        self.tol['Stress_abs_loose'] = 5.0e+3
        self.tol['Stress_abs']       = 1.0e-1
        self.tol['Strain_abs_loose'] = 5.0e-9
        self.tol['Strain_abs']       = 1.0e-10
        self.tol['Displacement']     = 1.0e-10

    def setDataStores(self):
        "defines testdata and goldendata output directories needed by this Capability TestCase"

        self.testdirs   = ['viscoplastic_ring_output']
        self.goldendirs = ['viscoplastic_ring_golden']

    def setDefinitions(self):
        "Define the data fields"

        self.whole_mesh      = self.testdata[0].getRegion(['Default'])
        self.sigma_ini       = self.testdata[0].getField(field='sigma', cycle=0, region= self.whole_mesh)
        self.epsilon_ini     = self.testdata[0].getField(field='epsilon', cycle=0, region= self.whole_mesh)
        self.epsdot_ini      = self.testdata[0].getField(field='epsdot', cycle=0, region= self.whole_mesh)
        self.ntrac_ini       = self.testdata[0].getField(field='NTRAC_04', cycle=0, region= self.whole_mesh)
        self.sigma_ini_ref   = self.goldendata[0].getField(field='sigma', cycle=0, region= self.whole_mesh)
        self.epsilon_ini_ref = self.goldendata[0].getField(field='epsilon', cycle=0, region= self.whole_mesh)
        self.epsdot_ini_ref  = self.goldendata[0].getField(field='epsdot', cycle=0, region= self.whole_mesh)
        self.ntrac_ini_ref   = self.goldendata[0].getField(field='NTRAC_04', cycle=0, region= self.whole_mesh)
        self.sigma_tr        = self.testdata[0].getField(field='sigma', cycle=5, region= self.whole_mesh)
        self.epsilon_tr      = self.testdata[0].getField(field='epsilon', cycle=5, region= self.whole_mesh)
        self.epsplas_tr      = self.testdata[0].getField(field='e_plastic', cycle=5, region= self.whole_mesh)
        self.sigma_tr_ref    = self.goldendata[0].getField(field='sigma', cycle=5, region= self.whole_mesh)
        self.epsilon_tr_ref  = self.goldendata[0].getField(field='epsilon', cycle=5, region= self.whole_mesh)
        self.epsplas_tr_ref  = self.goldendata[0].getField(field='e_plastic', cycle=5, region= self.whole_mesh)

    def testElasticStress(self):
         "Compare stresses at initial temperature"

         n = 0
         
         fail = 0
         for j in range(6):
             i = j + 1
             error = max(abs(self.sigma_ini.data[:,j]-self.sigma_ini_ref.data[:,j]))

             self.logger.write("  Cycle %2d, Sigma %2d: maximum absolute error = %10.4e\n" %(n,j,error))
             if error > self.tol['Stress_abs']:
                 self.logger.write("  Unexpected change in the stress %2d for cycle %2d: FAIL\n" %(i,n))
                 fail = fail + 1

         if fail:
             self.fail("Incorrect initial stress fields: %2d failed tests\n" %(fail))
         else:
             self.logger.write(": PASS \n")

    def testElasticStrain(self):
         "Compare strains at initial temperature"

         n = 0
         fail = 0

         for j in range(6):
             i = j + 1
             error = max(abs(self.epsilon_ini.data[:,j]-self.epsilon_ini_ref.data[:,j]))
    
             self.logger.write("  Cycle %2d, Epsilon %2d: maximum absolute error = %10.4e\n" %(n,j,error))
             if error > self.tol['Strain_abs']:
                 self.logger.write("  Incorrect change in initial strain %2d for cycle %2d: FAIL\n" %(i,n))
                 fail = fail + 1

         if fail:
             self.fail("Incorrect transient strain fields: %2d failed tests\n" %(fail))
         else:
             self.logger.write(": PASS \n")



    def testInitialTraction(self):
         "Compare interface normal traction at initial temperature"

         n = 0
         nnodes = 150
         error = 0.0

         error = max(abs(self.ntrac_ini.data - self.ntrac_ini_ref.data))
         self.logger.write("  Cycle %2d: maximum absolute normal traction error error = %10.4e\n" %(n,error))

         if error > self.tol['Stress_abs']:
             self.fail("Incorrect initial interface traction")
         else:
             self.logger.write(": PASS \n")

    def testFinalStress(self):
         "Compare stresses at final temperature"

         n = 8
         
         fail = 0
         for j in range(6):
             i = j + 1
             error = max(abs(self.sigma_tr.data[:,j]-self.sigma_tr_ref.data[:,j]))

             self.logger.write("  Cycle %2d, Sigma %2d: maximum absolute error = %10.4e\n" %(n,j,error))
             if error > self.tol['Stress_abs_loose']:
                 self.logger.write("  Unexpected change in the stress %2d for cycle %2d: FAIL\n" %(i,n))
                 fail = fail + 1

         if fail:
             self.fail("Incorrect final stress fields: %2d failed tests\n" %(fail))
         else:
             self.logger.write(": PASS \n")

    def testFinalStrain(self):
         "Compare strains at final temperature"

         n = 8
         fail = 0

         for j in range(6):
             i = j + 1
             error = max(abs(self.epsilon_tr.data[:,j]-self.epsilon_tr_ref.data[:,j]))
    
             self.logger.write("  Cycle %2d, Epsilon %2d: maximum absolute error = %10.4e\n" %(n,j,error))
             if error > self.tol['Strain_abs_loose']:
                 self.logger.write("  Incorrect change in final strain %2d for cycle %2d: FAIL\n" %(i,n))
                 fail = fail + 1

         if fail:
             self.fail("Incorrect final strain fields: %2d failed tests\n" %(fail))
         else:
             self.logger.write(": PASS \n")

    def testFinalPlasticStrain(self):
         "Compare effective plastic strains at final temperature"

         n = 8
         fail = 0

         epseff = Numeric.array([],'d')
         epseff_ref = Numeric.array([],'d')

         exx =  self.epsplas_tr.data[:,0]
         eyy =  self.epsplas_tr.data[:,1]
         ezz =  self.epsplas_tr.data[:,2]
         exy =  self.epsplas_tr.data[:,3]
         exz =  self.epsplas_tr.data[:,4]
         eyz =  self.epsplas_tr.data[:,5]

         epseff = (2.0/9.0 * ((exx - eyy)**2+(eyy - ezz)**2+(ezz - exx)**2) + 4.0/3.0 *(exy**2+exz**2+eyz**2))**0.5 

         exx =  self.epsplas_tr_ref.data[:,0]
         eyy =  self.epsplas_tr_ref.data[:,1]
         ezz =  self.epsplas_tr_ref.data[:,2]
         exy =  self.epsplas_tr_ref.data[:,3]
         exz =  self.epsplas_tr_ref.data[:,4]
         eyz =  self.epsplas_tr_ref.data[:,5]

         epseff_ref = (2.0/9.0 * ((exx - eyy)**2+(eyy - ezz)**2+(ezz - exx)**2) + 4.0/3.0 *(exy**2+exz**2+eyz**2))**0.5 

         error = max(abs(epseff-epseff_ref))
    
         self.logger.write("  Cycle %2d, Epsplas: maximum absolute error = %10.4e\n" %(n,error))
         if error > self.tol['Strain_abs_loose']:
             self.logger.write("  Incorrect change in final plasticstrain for cycle %2d: FAIL\n" %(n))
             fail = fail + 1

         if fail:
             self.fail("Incorrect final plastic strain fields: %2d failed tests\n" %(fail))
         else:
             self.logger.write(": PASS \n")


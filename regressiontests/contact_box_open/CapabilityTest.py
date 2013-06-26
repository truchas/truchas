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
    "Contact box close capability test"
    
    def setDataStores(self):
        "defines testdata and goldendata output directories needed by this Capability TestCase"

        self.testdirs    = ['contact_box_open_output']
        self.goldendirs  = ['contact_box_open_golden']

    def setDefinitions(self):
        self.whole_mesh = self.testdata[0].getRegion(['Default'])
        #self.gap = self.testdata[0].getRegion(['MB',3])
        self.gap = self.testdata[0].getRegion(['VAR', 'NTRAC_02', [-1.0e99, -1.0e-10]])
    
    def testOpenStress(self):
        "Verify stresses just after gap opening"
        self.logger.write("Comparing stresses to reference just after gap opening ...\n")
        n = 21
        tol = 1.0e3 # a little sloppy
        sigma = self.testdata[0].getField(field='sigma', cycle=n, region=self.whole_mesh)
        sigma_ref = self.goldendata[0].getField(field='sigma', cycle=n, region=self.whole_mesh)
        fail = 0
        for j in range(6):
            error = max(abs(sigma.data[:,j]-sigma_ref.data[:,j]))
            self.logger.write("  Cycle %2d, Sigma %2d: max abs error = %10.4e" %(n,j+1,error))
            if error > tol:
                self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
                fail = fail + 1
            else:
                self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
        if fail:
            self.fail("  Error exceeds the tolerance in %1d stress components" %(fail))

    def testOpenStrain(self):
        "Verify strains just after gap opening"
        self.logger.write("Comparing strains to reference just after gap opening ...\n")
        n = 21
        tol = 2.0e-8
        eps = self.testdata[0].getField(field='epsilon', cycle=n, region=self.whole_mesh)
        eps_ref = self.goldendata[0].getField(field='epsilon', cycle=n, region=self.whole_mesh)
        fail = 0
        for j in range(6):
            i = j + 1
            error = max(abs(eps.data[:,j]-eps_ref.data[:,j]))
            self.logger.write("  Cycle %2d, Epsilon %2d: max abs error = %10.4e" %(n,j+1,error))
            if error > tol:
                self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
                fail = fail + 1
            else:
                self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
        if fail:
            self.fail("  Error exceeds the tolerance in %1d stress components" %(fail))

    def testOpenNormalTraction(self):
        "Verifying normal traction just after gap opening"
        self.logger.write("Comparing normal traction just after gap opening to exact ...\n")
        n = 21
        tol = 1.0e-4
        exact = 0.0
        ntrac = self.testdata[0].getField(field='NTRAC_02', cycle=n, region=self.gap)
        error = max(abs(ntrac.data-exact))
        self.logger.write("  Cycle %2d: max abs error = %10.4e" %(n,error))
        if error > tol:
            self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
            self.fail("  Error exceeds the tolerance")
        else:
            self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))

    
    def testClosedStress(self):
        "Verify stresses just before gap opening"
        self.logger.write("Comparing stresses to reference just before gap opening ...\n")
        n = 18
        tol = 1.0e4 # sloppy
        sigma = self.testdata[0].getField(field='sigma', cycle=n, region=self.whole_mesh)
        sigma_ref = self.goldendata[0].getField(field='sigma', cycle=n, region=self.whole_mesh)
        fail = 0
        for j in range(6):
            error = max(abs(sigma.data[:,j]-sigma_ref.data[:,j]))
            self.logger.write("  Cycle %2d, Sigma %2d: max abs error = %10.4e" %(n,j+1,error))
            if error > tol:
                self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
                fail = fail + 1
            else:
                self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
        if fail:
            self.fail("  Error exceeds the tolerance in %1d stress components" %(fail))

    def testClosedStrain(self):
        "Verify strains just before gap opening"
        self.logger.write("Comparing strains to reference just before gap opening ...\n")
        n = 18
        tol = 2.0e-8
        eps = self.testdata[0].getField(field='epsilon', cycle=n, region=self.whole_mesh)
        eps_ref = self.goldendata[0].getField(field='epsilon', cycle=n, region=self.whole_mesh)
        fail = 0
        for j in range(6):
            i = j + 1
            error = max(abs(eps.data[:,j]-eps_ref.data[:,j]))
            self.logger.write("  Cycle %2d, Epsilon %2d: max abs error = %10.4e" %(n,j+1,error))
            if error > tol:
                self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
                fail = fail + 1
            else:
                self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
        if fail:
            self.fail("  Error exceeds the tolerance in %1d stress components" %(fail))

    def testClosedDisplacement(self):
        "Verify displacement magnitude just before gap opening"
        self.logger.write("Comparing displacement magnitude to reference just before gap opening ...\n")

        n = 18
        u = self.testdata[0].getField(field='disp', cycle=n, region=self.whole_mesh)
        u_ref = self.goldendata[0].getField(field='disp', cycle=n, region=self.whole_mesh)
        
        abs_tol = 1.0e-9
        rel_tol = 1.0e-4
        
        d = (u.data[:,0]**2 + u.data[:,1]**2 + u.data[:,2]**2)**0.5
        d_ref = (u_ref.data[:,0]**2 + u_ref.data[:,1]**2 + u_ref.data[:,2]**2)**0.5
        
        error = max(abs(d - d_ref)/(abs_tol + rel_tol*abs(d_ref)))
        self.logger.write("  Cycle %2d: max error = %10.4e" %(n,error))
        if error > 1.0:
            self.logger.write(": FAIL (tolerance = 1.0)\n")
            self.fail("  Error exceeds the tolerance")
        else:
            self.logger.write(": PASS (tolerance = 1.0)\n")
        
    def testClosedTemperature(self):
        "Verify temperature just before gap opening"
        self.logger.write("Comparing temperature to reference just before gap opening ...\n")
        n = 18
        tol = 1.0e-5
        T = self.testdata[0].getField(field='T', cycle=n, region= self.whole_mesh)
        T_ref = self.goldendata[0].getField(field='T', cycle=n, region= self.whole_mesh)
        error = max(abs((T.data-T_ref.data)/T_ref.data))
        self.logger.write("  Cycle %2d: max rel error = %10.4e" %(n,error))
        if error > tol:
            self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
            self.fail("  Error exceeds the tolerance")
        else:
            self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))

    def testInitialStress(self):
        "Verify initial stresses"
        self.logger.write("Comparing initial stresses to exact ...\n")
        n = 0
        tol = 1.0e2
        exact = [-2.288e7, -5.720e6, -1.716e7, 0.0, 0.0, -9.90733e6]
        sigma = self.testdata[0].getField(field='sigma', cycle=n, region= self.whole_mesh)
        fail = 0
        for j in range(6):
            error = max(abs(sigma.data[:,j]-exact[j]))
            self.logger.write("  Cycle %2d, Sigma %2d: max abs error = %10.4e" %(n,j+1,error))
            if error > tol:
                self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
                fail = fail + 1
            else:
                self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
        if fail:
            self.fail("  Error exceeds the tolerance in %1d stress components" %(fail))

    def testInitialStrain(self):
        "Verify initial strains"
        self.logger.write("Comparing initial strains to exact ...\n")
        n = 0
        tol = 1.0e-8
        exact = [0.0, 3.3e-4, 1.1e-4, 0.0, 0.0, -1.90526e-4]
        eps = self.testdata[0].getField(field='epsilon', cycle=n, region= self.whole_mesh)
        fail = 0
        for j in range(6):
            error = max(abs(eps.data[:,j]-exact[j]))
            self.logger.write("  Cycle %2d, Epsilon %2d: max abs error = %10.4e" %(n,j+1,error))
            if error > tol:
                self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
                fail = fail + 1
            else:
                self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
        if fail:
            self.fail("  Error exceeds the tolerance in %1d strain components" %(fail))

    def testFinalNormalTraction(self):
        "Verifying final normal traction"
        self.logger.write("Comparing final normal traction to exact ...\n")
        n = 49
        tol = 1.0e-4
        exact = 0.0
        ntrac = self.testdata[0].getField(field='NTRAC_02', cycle=n, region=self.gap)
        error = max(abs(ntrac.data-exact))
        self.logger.write("  Cycle %2d: max abs error = %10.4e" %(n,error))
        if error > tol:
            self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
            self.fail("  Error exceeds the tolerance")
        else:
            self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))

    def testFinalStress(self):
        "Verify final stresses"
        self.logger.write("Comparing final stresses to exact steady state ...\n")
        n = 49
        tol = 1.0e2
        exact = [7.62665e6, 1.90666e6, 5.71999e6, 3.81333e6, 6.60487e6, 3.30244e6]
        sigma = self.testdata[0].getField(field='sigma', cycle=n, region= self.whole_mesh)
        fail = 0
        for j in range(6):
            error = max(abs(sigma.data[:,j]-exact[j]))
            self.logger.write("  Cycle %2d, Sigma %2d: max abs error = %10.4e" %(n,j+1,error))
            if error > tol:
                self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
                fail = fail + 1
            else:
                self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
        if fail:
            self.fail("  Error exceeds the tolerance in %1d stress components" %(fail))

    def testFinalStrain(self):
        "Verify final strains"
        self.logger.write("Comparing final strains to exact steady state ...\n")
        n = 49
        tol = 1.0e-9
        exact = [-1.46667e-4, -2.56666e-4, -1.83333e-4, 7.33332e-5, 1.27017e-4, 6.35085e-5]
        eps = self.testdata[0].getField(field='epsilon', cycle=n, region= self.whole_mesh)
        fail = 0
        for j in range(6):
            error = max(abs(eps.data[:,j]-exact[j]))
            self.logger.write("  Cycle %2d, Epsilon %2d: max abs error = %10.4e" %(n,j+1,error))
            if error > tol:
                self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
                fail = fail + 1
            else:
                self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))
        if fail:
            self.fail("  Error exceeds the tolerance in %1d strain components" %(fail))

    def testInitialNormalTraction(self):
        "Verifying initial normal traction"
        self.logger.write("Comparing initial normal traction to exact steady state ...\n")
        n = 0
        tol = 1.0e-6
        exact = -2.288e7
        ntrac = self.testdata[0].getField(field='NTRAC_02', cycle=n, region=self.gap)
        error = max(abs((ntrac.data-exact)/exact))
        self.logger.write("  Cycle %2d: max rel error = %10.4e" %(n,error))
        if error > tol:
            self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
            self.fail("  Error exceeds the tolerance")
        else:
            self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))

    def testFinalTemperature(self):
        "Verifying final temperature"
        self.logger.write("Comparing final temperature field to exact steady state ...\n")
        n = 49
        tol = 1.0e-8
        exact = 288.0
        T = self.testdata[0].getField(field='T', cycle=n, region= self.whole_mesh)
        error = max(abs(T.data-exact)/exact)
        self.logger.write("  Cycle %2d: max rel error = %10.4e" %(n,error))
        if error > tol:
            self.logger.write(": FAIL (tolerance = %7.2e)\n" %(tol))
            self.fail("  Error exceeds the tolerance")
        else:
            self.logger.write(": PASS (tolerance = %7.2e)\n" %(tol))

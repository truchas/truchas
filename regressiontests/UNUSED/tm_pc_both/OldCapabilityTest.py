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
    "Thermo-mechanics phase change capability test"
    
    def setTolerances(self):
        "defines all tolerances to be used in this testcase"

        self.tol['Stress_rel']    = 1.0e-8
        self.tol['Stress_abs']    = 1.0e-1
        self.tol['Strain_rel']    = 1.0e-8
        self.tol['Strain_abs']    = 1.0e-10
        self.tol['Strain_loose']  = 1.0e-4
        self.tol['Temp_rel']      = 1.0e-10

    def setDataStores(self):
        "defines testdata and goldendata output directories needed by this Capability TestCase"

        self.testdirs   = ['tm_pc_both_output']
        self.goldendirs = ['tm_pc_both_golden']

    def setDefinitions(self):
        "Define the data fields"

        self.whole_mesh   = self.testdata[0].getRegion(['Default'])
        self.sigma        = self.testdata[0].getField(field='sigma', cycle=202, region= self.whole_mesh)
        self.epsilon      = self.testdata[0].getField(field='epsilon', cycle=202, region= self.whole_mesh)
        self.epspc        = self.testdata[0].getField(field='epspc', cycle=202, region= self.whole_mesh)
        self.epstherm     = self.testdata[0].getField(field='epstherm', cycle=202, region= self.whole_mesh)
        self.temp         = self.testdata[0].getField(field='T', cycle=202, region= self.whole_mesh)
        self.sigma_ref    = self.goldendata[0].getField(field='sigma', cycle=202, region= self.whole_mesh)
        self.epsilon_ref  = self.goldendata[0].getField(field='epsilon', cycle=202, region= self.whole_mesh)
        self.epspc_ref    = self.goldendata[0].getField(field='epspc', cycle=202, region= self.whole_mesh)
        self.epstherm_ref = self.goldendata[0].getField(field='epstherm', cycle=202, region= self.whole_mesh)
        self.temp_ref     = self.goldendata[0].getField(field='T', cycle=202, region= self.whole_mesh)

    def testFinalTemperature(self):
         "Compare final temperature to golden output"

         n = 202
         ncells = 18
         error = 0.0

         for j in range(ncells):
             diff = abs((self.temp.data[j] - self.temp_ref.data[j]) / self.temp_ref.data[j])
             if diff > error: error = diff
                     
         self.logger.write("  Cycle %2d: maximum relative temperature error = %10.4e\n" %(n,error))

         if error > self.tol['Temp_rel']:
             self.fail("Incorrect final temperature")
         else:
             self.logger.write(": PASS \n")

    def testFinalPhaseChangeStrain(self):
         "Compare final phase change strain to golden output"

         n = 202
         ncells = 18
         error = 0.0
         error_zero = 0.0
         
         for j in range(ncells):
             if self.epspc_ref.data[j,0] != 0.0:
                 diff = abs((self.epspc.data[j,0] - self.epspc_ref.data[j,0]) / self.epspc_ref.data[j,0])
                 if diff > error: error = diff
             else:
                 if self.epspc.data[j,0] != 0.0: error_zero = self.epspc.data[j,0]
                     
         self.logger.write("  Cycle %2d: maximum relative phase change strain error = %10.4e\n" %(n,error))
         self.logger.write("  Cycle %2d: maximum zero phase change strain error = %10.4e\n" %(n,error_zero))

         if (error > self.tol['Temp_rel']) or (error_zero != 0.0):
             self.fail("Incorrect phase change strain")
         else:
             self.logger.write(": PASS \n")

    def testFinalThermalStrain(self):
         "Compare final thermal strain to value calculated from temperature"

         n = 202
         ncells = 18
         cte1 = 2.2e-5
         cte2 = 2.1e-5
         cte3 = 2.0e-5
         T0 = 500.0
         Ti = 400.0
         Th = 375.0
         Tl = 350.0
         error = 0.0
         a = (cte2*Tl - cte3*Th)/(Tl - Th)
         b = (cte3 - cte2)/2.0/(Tl - Th)
         
#         self.epspc.str()

#   Here we compare with a calculated solution where the material has not seen
#   any non-isothermal phase change, and a golden solution where it has.
#   The non-isothermal thermal strain is a simple explicit integration of
#   a nonlinear equation.

         for j in range(ncells):
             T = self.temp.data[j]
             if T >= Ti:
                 epsref = (T-T0) * cte1
                 diff = abs((self.epstherm.data[j,0] - epsref) / epsref)
                 self.logger.write("  > Ti\n")
             elif T >= Th:
                 epsref = (Ti - T0) * cte1 + (T - Ti) * cte2
                 diff = abs((self.epstherm.data[j,0] - epsref) / epsref)
                 self.logger.write("  > Th\n")
             elif T >= Tl:
#                 epsref = (Ti - T0) * cte1 + (Th - Ti) * cte2
#                 epsref = epsref + a * (T - Th) + b * (T**2 - Th**2)
                 epsref = self.epstherm_ref.data[j,0]
                 diff = abs((self.epstherm.data[j,0] - epsref) / epsref)
                 self.logger.write("  > Tl\n")
             elif T < Tl:
#                 epsref = (Ti - T0) * cte1 + (Th - Ti) * cte2
#                 epsref = epsref + a * (Tl - Th) + b * (Tl**2 - Th**2)
#                 epsref = epsref + cte3 * (T - Tl)
                 epsref = self.epstherm_ref.data[j,0]
                 diff = abs((self.epstherm.data[j,0] - epsref) / epsref)
                 self.logger.write("  < Tl\n")

             self.logger.write("  Calc, reference %13.7e %13.7e %13.7e\n" %(self.epstherm.data[j,0], epsref, diff))
             if diff > error: error = diff
                     
         self.logger.write("  Cycle %2d: maximum relative thermal strain error = %10.4e\n" %(n,error))

         if error > self.tol['Strain_rel']:
             self.fail("Incorrect thermal strain")
         else:
             self.logger.write(": PASS \n")

    def testFinalStress(self):
         "Compare stresses at final temperature"

         n = 202

         epsdil = self.epspc.data[:,0] + self.epstherm.data[:,0]
         l1 = 5.20e+10
         l2 = 2.60e+10


         XXcomp = self.sigma.data[:,0]
         YYcomp = self.sigma.data[:,1]
         ZZcomp = self.sigma.data[:,2]
         XYcomp = self.sigma.data[:,3]
         XZcomp = self.sigma.data[:,4]
         YZcomp = self.sigma.data[:,5]

         XXref = -epsdil * (2.0*l1 + 2.0*l2 - 2.0*l1**2/(l1 + 2.0*l2))
         YYref = XXref

         errorxx = max(abs((XXcomp - XXref) / XXref))
         erroryy = max(abs((YYcomp - YYref) / YYref))
         errorzz = max(abs(ZZcomp))
         errorxy = max(abs(XYcomp))
         errorxz = max(abs(XZcomp))
         erroryz = max(abs(YZcomp))

         self.logger.write("  Cycle %2d: maximum relative sigxx error = %10.4e\n" %(n,errorxx))    
         self.logger.write("  Cycle %2d: maximum relative sigyy error = %10.4e\n" %(n,erroryy))    
         self.logger.write("  Cycle %2d: maximum absolute sigzz error = %10.4e\n" %(n,errorzz))
         self.logger.write("  Cycle %2d: maximum absolute sigxy error = %10.4e\n" %(n,errorxy))    
         self.logger.write("  Cycle %2d: maximum absolute sigxz error = %10.4e\n" %(n,errorxz))
         self.logger.write("  Cycle %2d: maximum absolute sigyz error = %10.4e\n" %(n,erroryz))

         fail = 0
         if errorxx > self.tol['Stress_rel']:
             fail = fail + 1
             self.logger.write("  Final sigxx not within %12e: FAIL\n" %(self.tol['Stress_rel']))
         if erroryy > self.tol['Stress_rel']:
             fail = fail + 1
             self.logger.write("  Final sigyy not within %12e: FAIL\n" %(self.tol['Stress_rel']))
         if errorzz > self.tol['Stress_abs']:
             fail = fail + 1
             self.logger.write("  Final sigxx not within %12e: FAIL\n" %(self.tol['Stress_abs']))
         if errorxy > self.tol['Stress_abs']:
             fail = fail + 1
             self.logger.write("  Final sigxy not within %12e: FAIL\n" %(self.tol['Stress_abs']))
         if errorxz > self.tol['Stress_abs']:
             fail = fail + 1
             self.logger.write("  Final sigxz not within %12e: FAIL\n" %(self.tol['Stress_abs']))
         if erroryz > self.tol['Stress_abs']:
             fail = fail + 1
             self.logger.write("  Final sigyz not within %12e: FAIL\n" %(self.tol['Stress_abs']))

         if fail:
             self.fail("Incorrect final stress fields: %2d failed tests\n" %(fail))
         else:
             self.logger.write(": PASS \n")

    def testFinalStrain(self):
         "Compare strains at final temperature"

         n = 202

         epsdil = self.epspc.data[:,0] + self.epstherm.data[:,0]
         l1 = 5.20e+10
         l2 = 2.60e+10


         XXcomp = self.epsilon.data[:,0]
         YYcomp = self.epsilon.data[:,1]
         ZZcomp = self.epsilon.data[:,2]
         XYcomp = self.epsilon.data[:,3]
         XZcomp = self.epsilon.data[:,4]
         YZcomp = self.epsilon.data[:,5]

         ZZref = epsdil + 2.0 * l1 * epsdil/(l1 + 2.0*l2)

         errorxx = max(abs(XXcomp))
         erroryy = max(abs(YYcomp))
         errorzz = max(abs((ZZcomp - ZZref) / ZZref))
         errorxy = max(abs(XYcomp))
         errorxz = max(abs(XZcomp))
         erroryz = max(abs(YZcomp))

         self.logger.write("  Cycle %2d: maximum absolute epsxx error = %10.4e\n" %(n,errorxx))    
         self.logger.write("  Cycle %2d: maximum absolute epsyy error = %10.4e\n" %(n,erroryy))    
         self.logger.write("  Cycle %2d: maximum relative epszz error = %10.4e\n" %(n,errorzz))
         self.logger.write("  Cycle %2d: maximum absolute epsxy error = %10.4e\n" %(n,errorxy))    
         self.logger.write("  Cycle %2d: maximum absolute epsxz error = %10.4e\n" %(n,errorxz))
         self.logger.write("  Cycle %2d: maximum absolute epsyz error = %10.4e\n" %(n,erroryz))

         fail = 0
         if errorxx > self.tol['Strain_abs']:
             fail = fail + 1
             self.logger.write("  Final epsxx not within %12e: FAIL\n" %(self.tol['Strain_abs']))
         if erroryy > self.tol['Strain_abs']:
             fail = fail + 1
             self.logger.write("  Final epsyy not within %12e: FAIL\n" %(self.tol['Strain_abs']))
         if errorzz > self.tol['Strain_rel']:
             fail = fail + 1
             self.logger.write("  Final epsxx not within %12e: FAIL\n" %(self.tol['Strain_rel']))
         if errorxy > self.tol['Strain_abs']:
             fail = fail + 1
             self.logger.write("  Final epsxy not within %12e: FAIL\n" %(self.tol['Strain_abs']))
         if errorxz > self.tol['Strain_abs']:
             fail = fail + 1
             self.logger.write("  Final epsxz not within %12e: FAIL\n" %(self.tol['Strain_abs']))
         if erroryz > self.tol['Strain_abs']:
             fail = fail + 1
             self.logger.write("  Final epsyz not within %12e: FAIL\n" %(self.tol['Strain_abs']))

         if fail:
             self.fail("Incorrect final strain fields: %2d failed tests\n" %(fail))
         else:
             self.logger.write(": PASS \n")


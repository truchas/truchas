#!/usr/bin/env python
if __name__=='__main__':

    import os, sys

    #developers : please specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '/../../')

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
  "DS1 CapabilityTest"

  def setDataStores(self):
    "defines testdata and goldendata output directories needed by this Capability TestCase"

    self.testdirs   = ['ht_bc_ortho_output', 'ht_bc_so_output']
    self.goldendirs = ['ht_bc_ortho_golden', 'ht_bc_so_golden']

  def setTolerances(self):
    "define tolerances"

    # here's how we can select for the compiler
    #if 'lahey' in self.basicrunspecs.compiler:
    #    print '\n ran with lahey \n'
    #if 'nag' in self.basicrunspecs.compiler:
    #    print '\n ran with nag \n'

    # we must choose a looser tolerance for the parallel
    # runs because the parallel runs will generate
    # different time step sequences histories, leading
    # to a different final time than in the serial runs
    # where the final time is equal to within machine epsilon
    if 'serial' in self.basicrunspecs.parallel_env:
        self.tol['temperature l infinity error'] = 5e-13
    else:
        self.tol['temperature l infinity error'] = 1e-8


  def setDefinitions(self):
    "Define the regions"
    self.reg_ortho = self.testdata[0].getRegion(['Default','all'])
    self.golden_reg_ortho = self.goldendata[0].getRegion(['Default','all'])
    
    self.reg_so = self.testdata[1].getRegion(['Default','all'])
    self.golden_reg_so = self.goldendata[1].getRegion(['Default','all'])
    

  def testHThtbc(self):
    "Verifying BCs in Heat Transfer"

    self.logger.write("Verifying BCs in Heat Transfer... ")

    temp_ortho = self.testdata[0].getField(field='T',
                                           cycle=10,
                                           region=self.reg_ortho)
    
    temp_ref_ortho = self.goldendata[0].getField(field='T',
                                                 cycle=10,
                                                 region=self.golden_reg_ortho)
    
    
    temp_so = self.testdata[1].getField(field='T',
                                        cycle=10,
                                        region=self.reg_so)
    
    temp_ref_so = self.goldendata[1].getField(field='T',
                                              cycle=10,
                                              region=self.golden_reg_so)
    
    



    TL2err_ortho = self.measure.linfError(temp_ortho.data,temp_ref_ortho.data)
    TL2err_so = self.measure.linfError(temp_so.data,temp_ref_so.data)
    

    fail = 0
    
    if TL2err_ortho > self.tol['temperature l infinity error']:
        self.logger.write(": FAIL(ortho) (tolerance = %12.7e)" %(TL2err_ortho))
        fail = fail + 1
    else:
        self.logger.write(": PASS(ortho) (tolerance = %12.7e)" %(TL2err_ortho))
        

    if TL2err_so > self.tol['temperature l infinity error']:
        self.logger.write(": FAIL(so) (tolerance = %12.7e)\n" %(TL2err_so))
        fail = fail + 1
    else:
        self.logger.write(": PASS(so) (tolerance = %12.7e)\n" %(TL2err_so))


    if fail > 0:
        self.fail("The temperature does not match the golden output.")
        

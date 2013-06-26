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

try:
    import numpy.oldnumeric as Numeric
except ImportError:
    import Numeric
except:
    raise

class CapabilityTest(TruchasCapabilityTest):

  def setTolerances(self):
      "define tolerances here"
    
      # Should see zero velocities, linear hydrostatic profile in pressure
      self.tol['spec_error']  = 1.0e-15


  def setDataStores(self):
      "Define data for hydrostatic regression test"
        
      self.testdirs     = ['species_simple_output']

  def setDefinitions(self):
      "Define the region and analytic centroid values"

      self.reg = self.testdata[0].getRegion(['Default','all'])

      self.Zhat = 0.06

  def testSpeciesConcentration(self):
      "Test species concentration along the characteristic"

      # Probe 1 at t=0.5 (5 cycles)
      e1 = self.testdata[0].getProbe(name='e1',
                                     field='PC1_1',
                                     timerange=[0.5,0.505]) 

      # Probe 2 at t=1.5 (15 cycles)
      e2 = self.testdata[0].getProbe(name='e2',
                                     field='PC1_1',
                                     timerange=[1.5,1.505]) 

      # Probe 3 at t=2.5 (25 cycles)
      e3 = self.testdata[0].getProbe(name='e3',
                                     field='PC1_1',
                                     timerange=[2.5,2.505]) 

      #PC probe data is always stored [cycle,time,PC] 
      dz1 = (e1.data[0,2] - self.Zhat)
      dz2 = (e2.data[0,2] - self.Zhat)
      dz3 = (e3.data[0,2] - self.Zhat)
      err = Numeric.sqrt(dz1*dz1 + dz2*dz2 + dz3*dz3)/self.Zhat

      # Test the error
      failmsg = 'In %s: species concentration error =  %10.4e is > %10.4e' \
                %(self.methodName, err, self.tol['spec_error'])
      self.failIf(err > self.tol['spec_error'], msg = failmsg)

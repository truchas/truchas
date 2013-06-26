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
    
      # Should see zero velocities, linear hydrostatic profile in pressure
      self.tol['velx_error']  = 2.00e-9
      self.tol['vely_error']  = 2.00e-9
      self.tol['velz_error']  = 2.00e-9
      self.tol['dpdz_error']  = 4.00e-7


  def setDataStores(self):
      "Define data for hydrostatic regression test"
        
      self.testdirs     = ['hydrostatic_steady_output']

  def setDefinitions(self):
      "Define the region and analytic centroid values"

      self.reg = self.testdata[0].getRegion(['Default','all'])

      self.dpdz = -9.81e+3

  def testSteadyHydrostaticFields(self):
      "Test velocity and pressure for the steady-hydrostatic problem"

      # Check everything at cycle=20
      test_cycle = 20
      
      vel = self.testdata[0].getField(field='Velocity',
                                      cycle=test_cycle,
                                      region=self.reg)

      velx_max = abs( self.measure.maxValue(vel.data[:,0]) )
      vely_max = abs( self.measure.maxValue(vel.data[:,1]) )
      velz_max = abs( self.measure.maxValue(vel.data[:,2]) )

      P = self.testdata[0].getField(field='P',
                                    cycle=test_cycle,
                                    region=self.reg)

      X   = self.testdata[0].getField(field='CENTROIDS', cycle=0)

      # Compute dpdz in the water using a unit length for the z-dimension
      # Use cell #55 and #7 (1-based)
      top = 54
      bot = 6
      dpdz = (P.data[top] - P.data[bot])/(X.data[top,2] - X.data[bot,2])

      # Do the trivial tests first
      # x-velocity
      failmsg = 'In %s: x-velocity error =  %10.4e is > %10.4e' \
                %(self.methodName, velx_max, self.tol['velx_error'])
      self.failIf(velx_max > self.tol['velx_error'], msg = failmsg)

      # y-velocity
      failmsg = 'In %s: y-velocity error =  %10.4e is > %10.4e' \
                %(self.methodName, vely_max, self.tol['vely_error'])
      self.failIf(vely_max > self.tol['vely_error'], msg = failmsg)

      # z-velocity
      failmsg = 'In %s: z-velocity error =  %10.4e is > %10.4e' \
                %(self.methodName, velz_max, self.tol['velz_error'])
      self.failIf(velz_max > self.tol['velz_error'], msg = failmsg)

      # Check the pressure gradient in the water
      # get % error = |P_umax - self.panalytic|/|self.panalytic|
      err = self.measure.percentError(dpdz,self.dpdz)

      logmessage = '\n Pressure gradient = %10.4e\n' % dpdz
      self.logger.write(logmessage)

      logmessage   = '\n Percentage pressure-gradient error = %10.4e \n' % err 
      self.logger.write(logmessage)

      failmessage  = 'In %s: Pressure gradient = %10.3e is > %10.3e' \
                     %(self.methodName, err, self.tol['dpdz_error'])
      self.failIf(err > self.tol['dpdz_error'], msg = failmessage)


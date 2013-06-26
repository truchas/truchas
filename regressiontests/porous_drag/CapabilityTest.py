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

class CapabilityTest(TruchasCapabilityTest):

  def setTolerances(self):
      "define tolerances here"
      
      # Should see constant x-velocity, linear pressure profile, and
      # zero velocities in the y and z-directions
      # self.tol['velx_error'] = 5.0e-10
      self.tol['velx_error'] = 5.0e-8
      self.tol['vely_error'] = 5.0e-10
      self.tol['velz_error'] = 5.0e-10
      self.tol[   'P_error'] = 5.0e-7

  def setDataStores(self):
      "Define data for the porous drag problem"

      self.testdirs     = ['porous_drag_output']

  def setDefinitions(self):
      "Define the region and analytic centroid values"

      self.reg = self.testdata[0].getRegion(['Default','all'])

      X = self.testdata[0].getField(field   = 'CENTROIDS', cycle   = 0)
        
      #define analytic pressure... 1000 at inlet, 0 at outlet
      length = 20.0
      Pin    = 20000.0
      self.Panalytic = Pin*(1.0 - X.data[:,0]/length)

      #define 1D analytic velocity -- y,z should be very small
      self.Uanalytic = 0.5


  def testVelocityandPressure(self):
      "Test the velocity and pressure fields"

      # Check everything at cycle=113
      test_cycle = 113

      #get pressure field 
      P  = self.testdata[0].getField(field  = 'P', cycle  = test_cycle)
      #P.str()
        
      #log P field
      self.logger.writeField(P, desc = 'Pressure in the entire domain')

      #get linf_error(P - Panalytic) 
      err         = self.measure.linfError(P.data,self.Panalytic)
      logmessage  = '\n linf P error = %10.3e \n' %(err)
      self.logger.write(logmessage)
	
      #fail if err < P_error
      failmessage = 'In %s: Pressure percent_error = %10.3e is > %10.3e' %(
          self.methodName, err, self.tol['P_error'])

      self.failIf(err > self.tol['P_error'], msg = failmessage)

      # Now check the velocities -- velx = 0.5 , vely = 0.0, velz = 0.0
      vel = self.testdata[0].getField(field='Velocity',
                                      cycle=test_cycle,
                                      region=self.reg)

      velx_min = abs( self.measure.minValue(vel.data[:,0]) )
      velx_max = abs( self.measure.maxValue(vel.data[:,0]) )
      vely_max = abs( self.measure.maxValue(vel.data[:,1]) )
      velz_max = abs( self.measure.maxValue(vel.data[:,2]) )

      # x-velocity
      velx_min_err = self.measure.percentError(velx_min, self.Uanalytic)
      failmsg = 'In %s: Minimum x-velocity error =  %10.4e, error is > %10.4e' \
                %(self.methodName, velx_min, velx_min_err)
      self.failIf(velx_min_err > self.tol['velx_error'], msg = failmsg)

      velx_max_err = self.measure.percentError(velx_max, self.Uanalytic)
      failmsg = 'In %s: Maximum x-velocity error =  %10.4e, error is > %10.4e' \
                %(self.methodName, velx_min, velx_max_err)
      self.failIf(velx_min_err > self.tol['velx_error'], msg = failmsg)

      # y-velocity
      failmsg = 'In %s: y-velocity error =  %10.4e is > %10.4e' \
                %(self.methodName, vely_max, self.tol['vely_error'])
      self.failIf(vely_max > self.tol['vely_error'], msg = failmsg)

      # z-velocity
      failmsg = 'In %s: z-velocity error =  %10.4e is > %10.4e' \
                %(self.methodName, velz_max, self.tol['velz_error'])
      self.failIf(velz_max > self.tol['velz_error'], msg = failmsg)

      msg = 'vel x min = %10.3e and err = %10.3e \n' % (velx_min, velx_min_err)
      self.logger.write(msg)
      msg = 'vel y min = %10.3e  \n' % (vely_max )
      self.logger.write(msg)
      msg = 'vel z min = %10.3e  \n' % (velz_max )
      self.logger.write(msg)

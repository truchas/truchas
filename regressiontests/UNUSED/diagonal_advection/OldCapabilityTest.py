#!/usr/bin/env python
if __name__=='__main__':

    import os, sys

    #specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '/../../')

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
  "Diagonal advection test"

  def setTolerances(self):
    "define tolerances"

    # Setup centroid, mass and enthalpy tolerances
    self.tol['xcent_eps'] = 1.0e-4
    self.tol['ycent_eps'] = 1.0e-11
    self.tol['zcent_eps'] = 2.0e-4
    self.tol['m_h2o_eps'] = 1.0e-11
    self.tol['m_air_eps'] = 1.0e-10
    self.tol['H_eps']     = 5.0e-5
    

  def setDataStores(self):
    "Define data for this regression test"

    self.testdirs   = ['diagonal_advection_output']


  def setDefinitions(self):
    "Define the region and analytic centroid values"

    self.reg = self.testdata[0].getRegion(['Default','all'])

    # Final centroid coordinates: area centroid and center of mass
    # This is based on (u=1,v=0,w=1) and tf=0.32
    self.xf = 0.62
    self.yf = 0.015625
    self.zf = 0.62


  def testPLICAdvection(self):
    "Checking PLIC Advection"

    self.logger.write("Checking PLIC Advection ...")

    # Check everything at cycle=40
    test_cycle = 40

    vof_h2o_i = self.testdata[0].getField(field='VOF0001',
                                          cycle=0,
                                          region=self.reg)

    vof_air_i = self.testdata[0].getField(field='VOF0002',
                                          cycle=0,
                                          region=self.reg)

    vof_h2o_f = self.testdata[0].getField(field='VOF0001',
                                          cycle=test_cycle,
                                          region=self.reg)
    
    vof_air_f = self.testdata[0].getField(field='VOF0002',
                                          cycle=test_cycle,
                                          region=self.reg)
    
    # Note: mixture enthalpy is carried in truchas
    h_i = self.testdata[0].getField(field='Enthalpy',
                                    cycle=0,
                                    region=self.reg)

    h_f = self.testdata[0].getField(field='Enthalpy',
                                    cycle=test_cycle,
                                    region=self.reg)

    X   = self.testdata[0].getField(field='CENTROIDS', cycle=0)

    # What's the cell volume called?
    # V   = self.testdata[0].getField(field='VOLUMES', cycle=0)

    # P.str() -- gives parameters
    # self.regionname.str() -- gives parameters

    # H20 density -- from input file
    rho_h2o = 1.258

    # Air density -- from input file
    rho_air = 0.1

    dx  = 1.0/32.0
    dy  = 0.03125
    dz  = dx
    vol = dx*dy*dz

    xcent_i = 0.0
    ycent_i = 0.0
    zcent_i = 0.0
    V_h2o_i = 0.0
    m_h2o_i = 0.0
    H_h2o_i = 0.0
    m_air_i = 0.0
    H_air_i = 0.0
    i       = 0
    for vf in vof_h2o_i.data:
       m_h2o_i = m_h2o_i + rho_h2o*vol*vf
       m_air_i = m_air_i + rho_air*vol*vof_air_i.data[i]
       H_h2o_i = H_h2o_i + rho_h2o*vol*vf*h_i.data[i]
       H_air_i = H_air_i + rho_air*vol*vof_air_i.data[i]*h_i.data[i]
       V_h2o_i = V_h2o_i + vol*vf
       xcent_i = xcent_i + X.data[i,0]*vol*vf
       ycent_i = ycent_i + X.data[i,1]*vol*vf
       zcent_i = zcent_i + X.data[i,2]*vol*vf
       i = i + 1

    xcent_i = xcent_i/V_h2o_i
    ycent_i = ycent_i/V_h2o_i
    zcent_i = zcent_i/V_h2o_i

    xcent_f = 0.0
    ycent_f = 0.0
    zcent_f = 0.0
    V_h2o_f = 0.0
    m_h2o_f = 0.0
    H_h2o_f = 0.0
    m_air_f = 0.0
    H_air_f = 0.0
    i       = 0
    for vf in vof_h2o_f.data:
       m_h2o_f = m_h2o_f + rho_h2o*vol*vf
       m_air_f = m_air_f + rho_air*vol*vof_air_f.data[i]
       H_h2o_f = H_h2o_f + rho_h2o*vol*vf*h_f.data[i]
       H_air_f = H_air_f + rho_air*vol*vof_air_f.data[i]*h_f.data[i]
       V_h2o_f = V_h2o_f + vol*vf
       xcent_f = xcent_f + X.data[i,0]*vol*vf
       ycent_f = ycent_f + X.data[i,1]*vol*vf
       zcent_f = zcent_f + X.data[i,2]*vol*vf
       i = i + 1

    xcent_f = xcent_f/V_h2o_f
    ycent_f = ycent_f/V_h2o_f
    zcent_f = zcent_f/V_h2o_f

    m_h2o_err = abs(1.0 - m_h2o_f/m_h2o_i)
    m_air_err = abs(1.0 - m_air_f/m_air_i)
    H_h2o_err = abs(1.0 - H_h2o_f/H_h2o_i)
    H_air_err = abs(1.0 - H_air_f/H_air_i)
    H_tot_i = H_h2o_i + H_air_i
    H_tot_f = H_h2o_f + H_air_f
    H_tot_err = abs(1.0 - H_tot_f/H_tot_i)
    xcent_err = abs(1.0 - xcent_f/self.xf)
    ycent_err = abs(1.0 - ycent_f/self.yf)
    zcent_err = abs(1.0 - zcent_f/self.zf)
    
    # Do the error tests now
    # h2o mass
    failmsg = 'In %s: h2o mass error =  %10.4e is > %10.4e' \
              %(self.methodName, m_h2o_err, self.tol['m_h2o_eps'])
    self.failIf(m_h2o_err > self.tol['m_h2o_eps'], msg = failmsg)

    # air mass
    failmsg = 'In %s: air mass error =  %10.4e is > %10.4e' \
              %(self.methodName, m_air_err, self.tol['m_air_eps'])
    self.failIf(m_air_err > self.tol['m_air_eps'], msg = failmsg)

    # x-centroid
    failmsg = 'In %s: x-centroid error =  %10.4e is > %10.4e' \
              %(self.methodName, xcent_err, self.tol['xcent_eps'])
    self.failIf(xcent_err > self.tol['xcent_eps'], msg = failmsg)

    # y-centroid
    failmsg = 'In %s: y-centroid error =  %10.4e is > %10.4e' \
              %(self.methodName, ycent_err, self.tol['ycent_eps'])
    self.failIf(ycent_err > self.tol['ycent_eps'], msg = failmsg)

    # z-centroid
    failmsg = 'In %s: z-centroid error =  %10.4e is > %10.4e' \
              %(self.methodName, zcent_err, self.tol['zcent_eps'])
    self.failIf(zcent_err > self.tol['zcent_eps'], msg = failmsg)

    # Total enthalpy
    failmsg = 'In %s: total enthalpy error =  %10.4e is > %10.4e' \
              %(self.methodName, H_tot_err, self.tol['H_eps'])
    self.failIf(H_tot_err > self.tol['H_eps'], msg = failmsg)

#!/usr/bin/env python

try:
    import numpy.oldnumeric as Numeric
    from numpy.oldnumeric import *
except ImportError:
    from Numeric import *
except:
    raise
from math import *

def ScriptSolutionTime(npts):

    time  = zeros(npts,Float64)
    input = open('script_data.dat','r')
    for k in range(npts):
        S       = input.readline()
        cols    = S.split()
        time[k] = float(cols[2])

    return time

def ScriptSolutionCell1(npts):

    cell1 = zeros(npts,Float64)
    input = open('script_data.dat','r')
    for k in range(npts):
        S        = input.readline()
        cols     = S.split()
        cell1[k] = float(cols[4])
        
    return cell1

def ScriptSolutionCell5(npts):

    cell5 = zeros(npts,Float64)
    input = open('script_data.dat','r')
    for k in range(npts):
        S        = input.readline()
        cols     = S.split()
        cell5[k] = float(cols[6])
        
    return cell5

def ScriptSolutionCell10(npts):

    cell10 = zeros(npts,Float64)
    input  = open('script_data.dat','r')
    for k in range(npts):
        S         = input.readline()
        cols      = S.split()
        cell10[k] = float(cols[8])
        
    return cell10


if __name__=='__main__':

    import os, sys

    #developers : please specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '/../..')

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
  "Time Dependent Thermal BC CapabilityTest"

  def setDataStores(self):
    "defines testdata  directory needed by this Capability TestCase"

    self.testdirs   = ['t_dep_temp_bc_output']

  def setTolerances(self):
    "define tolerances"

    self.tol['Time_Error']        = 1e-8
    self.tol['Temperature_Error'] = 3e-3

  def setDefinitions(self):
    "Define the regions"
    self.reg = self.testdata[0].getRegion(['Default'])

  def testTime_Points(self):
    "Verifying time values for probe data"

    self.logger.write("Verifying time values for probe data")

    time_truchas  = self.testdata[0].getProbe(name='cell1',field='TEMP') #,timerange=[0.0,2.0])
    npts          = len(time_truchas.data)
    time_script   = ScriptSolutionTime(npts)
    truchas_times = zeros(npts,Float64)
    script_times  = zeros(npts,Float64)
    
    #  Check Time values from Probe and Script Results
    truchas_times[:] = time_truchas.data[:,1]
    script_times[:]  = time_script[:]

    err              = self.measure.linfError(truchas_times,script_times)
    logmessage       = '\n Max Error in time values = %10.3e \n' %(err)
    self.logger.write(logmessage)
    
    #fail if err < self.tol['Time_Error']
    failmessage  = 'In %s: Max Error in time values = %10.3e is > %10.3e' %(self.methodName, err, self.tol['Time_Error'])

    self.failIf(err > self.tol['Time_Error'],
                msg = failmessage)
    

  def testCell_1_Temps(self):
    "Verifying Cell 1 Temperatures"

    self.logger.write("Verifying Cell 1 Temperatures")
    
    temp_truchas  = self.testdata[0].getProbe(name='cell1',field='TEMP') #,timerange=[0.0,2.0])
    npts          = len(temp_truchas.data)
    temp_script   = ScriptSolutionCell1(npts)
    truchas_temps = zeros(npts,Float64)
    script_temps  = zeros(npts,Float64)
    
    #  Check Time values from Probe and Script Results
    truchas_temps[:] = temp_truchas.data[:,2]
    script_temps[:]  = temp_script[:]

    err          = self.measure.linfError(truchas_temps,script_temps)
    logmessage   = '\n Max Error in temperature = %10.3e \n' %(err)
    self.logger.write(logmessage)
    
    #fail if err < self.tol['Temperature_Error']
    failmessage  = 'In %s: Max Error in temperature = %10.3e is > %10.3e' %(self.methodName, err, self.tol['Temperature_Error'])

    self.failIf(err > self.tol['Temperature_Error'],
                msg = failmessage)
    

  def testCell_5_Temps(self):
    "Verifying Cell 5 Temperatures"

    self.logger.write("Verifying Cell 5 Temperatures")
    
    temp_truchas  = self.testdata[0].getProbe(name='cell5',field='TEMP') #,timerange=[0.0,2.0])
    npts          = len(temp_truchas.data)
    temp_script   = ScriptSolutionCell5(npts)
    truchas_temps = zeros(npts,Float64)
    script_temps  = zeros(npts,Float64)
    
    #  Check Time values from Probe and Script Results
    truchas_temps[:] = temp_truchas.data[:,2]
    script_temps[:]  = temp_script[:]

    err          = self.measure.linfError(truchas_temps,script_temps)
    logmessage   = '\n Max Error in temperature = %10.3e \n' %(err)
    self.logger.write(logmessage)
    
    #fail if err < self.tol['Temperature_Error']
    failmessage  = 'In %s: Max Error in temperature = %10.3e is > %10.3e' %(self.methodName, err, self.tol['Temperature_Error'])

    self.failIf(err > self.tol['Temperature_Error'],
                msg = failmessage)
    

  def testCell_10_Temps(self):
    "Verifying Cell 10 Temperatures"

    self.logger.write("Verifying Cell 10 Temperatures")
    
    temp_truchas  = self.testdata[0].getProbe(name='cell10',field='TEMP') #,timerange=[0.0,2.0])
    npts          = len(temp_truchas.data)
    temp_script   = ScriptSolutionCell10(npts)
    truchas_temps = zeros(npts,Float64)
    script_temps  = zeros(npts,Float64)
    
    #  Check Time values from Probe and Script Results
    truchas_temps[:] = temp_truchas.data[:,2]
    script_temps[:]  = temp_script[:]

    err          = self.measure.linfError(truchas_temps,script_temps)
    logmessage   = '\n Max Error in temperature = %10.3e \n' %(err)
    self.logger.write(logmessage)
    
    #fail if err < self.tol['Temperature_Error']
    failmessage  = 'In %s: Max Error in temperature = %10.3e is > %10.3e' %(self.methodName, err, self.tol['Temperature_Error'])

    self.failIf(err > self.tol['Temperature_Error'],
                msg = failmessage)
"""
THE FOLLOWING TEXT IS A DOCSTRING IT IS IGNORED WITHIN THIS PYTHON SOURCE FILE

THE FOLLOWING IS THE SCRIPT FOR GENERATING THE DATA FILE USED BY THIS CAPABILITY TEST
TO EXECUTE IT USE THE FOLLOWING COMMAND:

python Newton.py > script_data.dat


from Numeric import *
from math import *
import LinearAlgebra

class TDepConductionData:

    def __init__(self):
	self.k   = 10.
	self.rho = 100.
	self.Cp  = 10.
        self.sigma = 5.670e-8  # stefan-boltzmann constant
        self.epsilon = 1.0     # surface emmissivity

	self.Tinit = 50.
	self.TbcL  = 2.0
        self.TbcR  = 50.
        
	self.L     = 1.0   # domain length, not latent heat
	self.dx    = 0.
	self.time  = 10.
	self.dt    = 0.01

	self.tol = 1e-11

def EvalFunc(data,Ncell,h_inc,h,residual):

    F_h   = zeros((Ncell),Float64)
    H_np1 = h_inc + h

    # evaluate F(H)
    a = data.k/data.rho/data.dx/data.dx
    Nusselt = data.epsilon*data.sigma*data.dx/2.0/data.k
    T0  = EvalTemp(data,H_np1[0])
    T1  = EvalTemp(data,H_np1[1])
    T2  = EvalTemp(data,H_np1[2])
    # determine the interface temperature at the radiation cell face
    T_interface = T1+Nusselt*T0**4
    itmax=125
    for i in range(itmax):
        rinner = T_interface*(1.0+Nusselt*T_interface**3)-(T1+Nusselt*T0**4)
        deltaT = -rinner/(1.0+4.0*Nusselt*T_interface**3)
        T_interface = T_interface+deltaT
        Qcond = -2.0*data.k*(T1-T_interface)/data.dx
        Qrad = data.epsilon*data.sigma*(T0**4-T_interface**4)
        dQ = Qrad - Qcond
        TestQ = data.tol*max(abs(Qcond)+abs(Qrad),1.0)
        if abs(dQ)<TestQ: break
    else:
        print T_interface,T0,T1,Qrad,Qcond,dQ
        raise 'EvalFunc ERROR: Exiting interface temperature solve'
    
    F_h[1]     = h_inc[1]/data.dt + data.epsilon*data.sigma/data.rho/data.dx*(T_interface**4-T0**4) - a*(T2-T1)
    residual[1]=-F_h[1]
    
    T_i = T1
    T_ip1 = T2
    for i in range(2,Ncell-1):
	T_im1  = T_i
	T_i    = T_ip1
	T_ip1  = EvalTemp(data,H_np1[i+1])
	F_h[i] = (h_inc[i])/data.dt - a*(T_ip1 - 2*T_i + T_im1)
            
	residual[i] = -F_h[i]

    # BC residuals
    residual[0] = 0
    residual[Ncell-1]=0

def TangFunc(data,NCell,h_inc,h,jacobian):

    for i in xrange(len(jacobian)):
	for j in xrange(len(jacobian)):
	    jacobian[i,j] = 0.0

    a = data.k/data.rho/data.dx/data.dx
    Nusselt = data.epsilon*data.sigma*data.dx/2.0/data.k
    RadK = 4.0*data.epsilon*data.sigma/data.rho/data.Cp/data.dx
    H_np1 = h_inc + h

    T0  = EvalTemp(data,H_np1[0])
    T1  = EvalTemp(data,H_np1[1])
    # determine the interface temperature at the radiation cell face
    T_interface = T1+Nusselt*T0**4
    itmax=25
    for i in range(itmax):
        rinner = T_interface*(1.0+Nusselt*T_interface**3)-(T1+Nusselt*T0**4)
        deltaT = -rinner/(1.0+4.0*Nusselt*T_interface**3)
        T_interface = T_interface+deltaT
        Qcond = -2.0*data.k*(T1-T_interface)/data.dx
        Qrad = data.epsilon*data.sigma*(T0**4-T_interface**4)
        dQ = Qrad - Qcond
        TestQ = data.tol*max(abs(Qcond)+abs(Qrad),1.0)
        if abs(dQ)<TestQ: break
    else:
        print T_interface,T0,T1,Qrad,Qcond,dQ
        raise 'TangFunc ERROR: Exiting interface temperature solve'
    
    jacobian[1,0] = 0.
    jacobian[1,2] = -a*dTemp(data,H_np1[2])
    jacobian[1,1] = 1/data.dt - RadK*T_interface**3/(1.0+4.0*Nusselt*T_interface**3) + a*dTemp(data,H_np1[2])
    for i in range(2,Ncell-1):
        jacobian[i,i-1] = -a*dTemp(data,H_np1[i-1])
        jacobian[i,i+1] = -a*dTemp(data,H_np1[i+1])
        jacobian[i,i] = 1/data.dt + 2.*a*dTemp(data,H_np1[i])
    
    # BC equations    
    jacobian[0,:] = 0.
    jacobian[:,0] = 0.
    jacobian[0,0] = 1.
    jacobian[Ncell-1,:]=0.
    jacobian[:,Ncell-1]=0.
    jacobian[Ncell-1,Ncell-1]=1.

def EvalTemp(data,h):

    Temp = h/data.Cp
    return Temp

def dTemp(data,h):

    dT = 1./data.Cp
    return dT

def Solve(H,Ncell,tol,data):

    # initiailize
    del_H_itr = zeros((Ncell),Float64)
    del_H_cum = zeros((Ncell),Float64)
    
    func_val  = zeros((Ncell),Float64)
    jacobian  = zeros((Ncell,Ncell),Float64)
    iter      = 0
    converge  = False
    iter      = 0
    MAX_ITER  = 25 

    # evaluate the nonlinear function
    EvalFunc(data,Ncell,del_H_cum,H,func_val)

    # iterate over step solution
    while converge == False and iter < MAX_ITER:

        iter += 1

	# solve to find newton step
	TangFunc(data,Ncell,del_H_cum,H,jacobian)

	del_H_itr = LinearAlgebra.solve_linear_equations(jacobian,func_val)
	del_H_cum += del_H_itr

	EvalFunc(data,Ncell,del_H_cum,H,func_val)

	norm = 0.0
	for i in xrange(Ncell): norm += func_val[i]**2
	norm = sqrt(norm)

	if norm <= tol: converge = True

        # update the conduction boundary condition
        H[Ncell-1]=2.*data.TbcR*data.Cp-(H[Ncell-2]+del_H_cum[Ncell-2])
        T[Ncell-1]=H[Ncell-1]/data.Cp
        TNgbr = (H[Ncell-2]+del_H_cum[Ncell-2])/data.Cp
        
    if iter > MAX_ITER-1:
            print '    iter: ',iter,' norm: ',norm
            raise ' Iteration failure in SOLVE'

    # new vector
    H += del_H_cum
    return H

if __name__=='__main__':


    model = TDepConductionData()
    # account for ghost cells to do bcs
    Ncell = 10
    model.dx = model.L/Ncell
    model.L += 2*model.dx
    Ncell += 2

    steps = int(model.time/model.dt + 1)
    
    # initialize
    H = ones((Ncell),Float64)*(model.Tinit*model.Cp)
    T = ones((Ncell),Float64)*model.Tinit

    print 'time =  0.0 ',' T[1]  ',T[1],' T[5]  ',T[5],' T[10]  ',T[10]
    
    count = 0
    # start solution: loop over time steps
    for t in range(0,steps):

	# set BC's	
        # Left (radiation) boundary
        model.TbcL = 50.0+(t*model.dt-0.0)*(300.0-50.0)/(10.0-0.0)
        H[0] = model.TbcL*model.Cp
        T[0] = H[0]/model.Cp
        # Right (conduction) boundary
        if t*model.dt<1.0:
            model.TbcR = 50.0+(t*model.dt-0.0)*(1.0-50.0)/(1.0-0.0)
        elif t*model.dt<2.0:
            model.TbcR = 1.0+(t*model.dt-1.0)*(100.0-1.0)/(2.0-1.0)
        elif t*model.dt<5.0:
            model.TbcR = 100.0+(t*model.dt-2.0)*(25.0-100.0)/(5.0-2.0)
        elif t*model.dt<7.5:
            model.TbcR = 25.0+(t*model.dt-5.0)*(222.0-25.0)/(7.5-5.0)
        else:
            model.TbcR = 222.0+(t*model.dt-7.5)*(0.5-222.0)/(10.0-7.5)
            
        H[Ncell-1]=2.*model.TbcR*model.Cp-H[Ncell-2]
        T[Ncell-1]=H[Ncell-1]/model.Cp
        
	# solver for next H
	H = Solve(H,Ncell,model.tol,model)

        # copy and update state
	for i in xrange(Ncell):
	    T[i]     = EvalTemp(model,H[i])

	# print output
        count += 1
	if count == 10:
	
	    print 'time = ',(t+1)*model.dt,' T[1]  ',T[1],' T[5]  ',T[5],' T[10]  ',T[10]
	    count =0
THIS IS THE END OF THE DOCSTRING
"""    

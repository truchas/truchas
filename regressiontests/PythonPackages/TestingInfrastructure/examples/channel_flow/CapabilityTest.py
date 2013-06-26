try:
    import numpy.oldnumeric as Numeric
except ImportError:
    import Numeric
except:
    raise

def AnalyticSolution(Ncellx,Ncelly):
        

    # Analytic solution to balance between viscous stress and pressure gradient
    #      Assumes equal size cells along each axis
    Pinlet = 1000.
    Poutlet = 0.
    Length = 1.
    GradP = (Poutlet-Pinlet)/Length
    Rho = 1.
    Mu  = 1.
    Height   = 1.
    DeltaY = Height/Ncelly
    RHS = DeltaY*DeltaY*GradP/Mu
    Beta=-3.
    BCoeff=-2.
    Ux=Numeric.array([RHS/Beta],'d')
    Gamma=Numeric.array([],'d')

    #Forward Elimination Loop
    for j in range(1,Ncelly):
        if j==Ncelly-1:
            BCoeff=-3.
        Gamma.resize([j+1])
        Gamma[j]=1./Beta
        Beta=BCoeff-Gamma[j]
        Ux.resize([j+1])
        Ux[j]=(RHS-Ux[j-1])/Beta

    #Back Substitution Loop
    for i in range(1,Ncelly):
        j=Ncelly-i-1
        Ux[j]=Ux[j]-Gamma[j+1]*Ux[j+1]
    return Ux

if __name__=='__main__':

    import os, sys
    #developers : please specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '../../../../../../')

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


    def setTolerances(self):
        "define tolerances here"
        
        self.tol['VL1_error']  = 1.0e-04
        self.tol['P_error']    = 1.0e-04


    def setDataStores(self):
        "defines all datastores to be used throughout this TestCase"

        #define testdata output directories..these result from running a BasicRun TestCase
        self.testdirs     = ['channel_flow_restart_output']

    def setDefinitions(self):
        "defines analytic expressions and regions needed by this Capability TestCase"

        #get region representing the entire mesh
        self.meshregion  = self.testdata[0].getRegion(selector=['Default'],
                                                      desc    = 'Region representing the entire Truchas mesh')

    def testXVelocities(self):
	"tests x direction velocity generation in entire domain"

        # Define cell counts
        Ncellx = 2
        Ncelly = 7
        
 	#get velocities on the entire domain
	V     = self.testdata[0].getField(field  = 'Velocity',
                                          region = self.meshregion,
                                          cycle  = 51025)

	Ucomp = V.data[:,0]
        # Skip cells in data because the mesh is Ncellx wide
        i=0
        Ucomp2=Numeric.array([],'d')
        for j in range(Ncelly):
            Ucomp2.resize([j+1])
            Ucomp2[j]=Ucomp[i]
            i=i+Ncellx

        #log V field
	self.logger.writeField(V,
                               region = self.meshregion,
                               desc   = 'Velocity at cycle 1')

        # get Analytic Solution
        U_Analytic=AnalyticSolution(Ncellx,Ncelly)
        
	#get L1 error
	L1err    = self.measure.l1Error(Ucomp2,U_Analytic)

        #fail if L1err > VL1_error
        failmessage = 'In %s velocity L1_error = %10.3e is > %10.3e' %(self.methodName, L1err, self.tol['VL1_error'])

	self.failIf(L1err > self.tol['VL1_error'],
                    msg = failmessage)

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

try:
   import numpy.oldnumeric as Numeric
except ImportError:
   import Numeric
except:
   raise

import string

from TestCases import TruchasCapabilityTest

class CapabilityTest(TruchasCapabilityTest):
    "composite conduction CapabilityTest" 

    def setTolerances(self):
        "define tolerances here"
        
        self.tol['T_error']   = 4.0
        self.tol['Prb_error'] = 2.0e-2

    def setDataStores(self):
        "defines testdata and goldendata output directories needed by this Capability TestCase"

        #define testdata output directories..these are created from running a BasicRun TestCase
        self.testdirs     = ['composite_conduction_output']

        #define restartdata output directories..these are created from running a RestartRun TestCase
        self.restartdirs  = ['composite_conduction_restart_output']

    def setDefinitions(self):
	"developer defined regions to be used throughout this Capability TestCase"

	#get region representing plane y=0.5 from testdata store
	self.myplane      = self.testdata[0].getRegion(selector = ['S','[0.0,0.5,0.0->1.0,0.5,1.0]'],
                                                        desc     = 'Region:y=0.5 plane')

        #get probe at center of circle (0.5,0.5,0.05)
        self.centerprobe  = self.restartdata[0].getProbe(name  = 'Center',
                                                          field = 'TEMP')
                                                       
    def testTemperaturesOnPlane(self):
	"tests temperatures along y=0.5 at various times 0.2<=t<=1.0"

	#get list of timesteps from Truchas testdata output that satisfy 0.2<=t<=1.0 
	these_timesteps = self.testdata[0].getTimeSteps(timerange = [0.2,1.0],
                                                         dtstep   = 0.2)

        #get specialised 'Comsol' golden output of temperature along y=0.5
        mygolden        = self.__getMyGoldenOutput(file = 'Tvsx.txt')

        #get cell-center X coordinates along y=0.5
        Xplane          = self.testdata[0].getField(field    = 'CENTROIDS',
                                                    region   = self.myplane)
	for timestep in these_timesteps:

	    thistime    = timestep.time
            index       = these_timesteps.index(timestep)

            #get Temperature field from testdata store on self.myplane region
	    Temp        = self.testdata[0].getField(field    = 'T',
                                                    time     = thistime,
                                                    region   = self.myplane)
            
	    #log Temperature field
	    description = 'Test Temperature at time t=%5.2f' %(thistime)
	    self.logger.writeField(Temp,
                                   region = self.myplane,
                                   desc   = description)

            Temp_golden = self.__interpolateData(Ydata   = mygolden[:,index+1],
                                                 Xdata   = mygolden[:,0],
                                                 Xinterp = Xplane.data[:,0])

	    #calculate the % error 
	    err         = self.measure.l1Error(Temp.data,Temp_golden)

	    #fail if err > T_error
            failmessage = 'In %s: Temperature L1_error = %10.3e is > %10.3e' %(self.methodName, err, self.tol['T_error'])

	    self.failIf(err>self.tol['T_error'],msg = failmessage)


    def testTemperaturesAtCenter(self):
        "tests temperatures at the circle center (0.5,0.5,0.05)"

        #get specialised 'Comsol' golden output of temperature at center
        mygolden    = self.__getMyGoldenOutput(file     = 'Tcentervst.txt',
                                               rowstart = 0,
                                               columns  = 2)

        Temp_golden = self.__interpolateData(Ydata   = mygolden[:,1],
                                             Xdata   = mygolden[:,0],
                                             Xinterp = self.centerprobe.data[:,1])


        #calculate the % error 
        err         = self.measure.l1Error(self.centerprobe.data[:,2],Temp_golden)

        #fail if err > Prb_error
        failmessage = 'In %s: Center probe temperature L1_error = %10.3e is > %10.3e' %(self.methodName, err, self.tol['Prb_error'])
        
        self.failIf(err>self.tol['Prb_error'],msg = failmessage)


    def __getMyGoldenOutput(self,rowstart=2,columns=6,file='Tvsx.txt'):
        "reads in golden data files generated from Comsol"
        
        fp     = open(file,'r')
        mydata = fp.readlines()
        fp.close()

        T     = Numeric.array(0.0,'d')
        shape = [len(mydata[rowstart:]),columns]
        T     = Numeric.resize(T,shape)

        cnti  = 0
        for i in mydata[rowstart:]:
            tmp   = string.split(i,'\n')
            ii    = cnti
            cnti += 1
            cntj  = 0
            try:
                tmp2  = string.split(tmp[0],'\t')
                for a in tmp2:
                    if not len(a):tmp2.remove(a)
                for j in tmp2:
                    T[ii,cntj] = float(j)
                    cntj += 1
            except:
                sep   = ' '*10
                tmp2  = string.split(tmp[0],sep)
                tmp3  = []
                for a in tmp2:
                    if '.' in a:
                        tmp3.append(a)
                for j in tmp3:
                    T[ii,cntj] = float(j)
                    cntj      += 1

        return T


    def __interpolateData(self,Ydata,Xdata,Xinterp):
        "given X,Y data and Xinterp positions to interpolate to, get interpolated field"

        Yinterp = Numeric.array(0.0,'d')
        Yinterp = Numeric.resize(Yinterp,Numeric.shape(Xinterp))
            
        cnt     = 0
        for x in Xinterp:
            xa    = Numeric.choose(Numeric.less_equal(Xdata,x),(Xdata,0.0))
            if len(Numeric.nonzero(xa)):
                na   = Numeric.nonzero(xa)[0]
                inbounds = (x >= Xdata[na-1] and x <= Xdata[na])
                err      = 'Unable to linearly interpolate data for x=%d' %(x)
                assert inbounds, err
                Yinterp[cnt] = Ydata[na-1] + (x-Xdata[na-1])*(Ydata[na]-Ydata[na-1])/(Xdata[na]-Xdata[na-1]) 
            else:
                na           = len(Xdata)-1
                Yinterp[cnt] = Ydata[na]

            cnt += 1

        return Yinterp

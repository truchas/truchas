if __name__=='__main__':

    import os, sys
    #developers : please specify the location of your truchas checkout directory 'truchasdir'
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    truchasdir  = os.path.abspath(thisdir + '../../../../../')

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
    RunThisCapability(testdir,
                      loadername = 'MyAbortLoader',
                      runnername = 'RunHierarchicalTestSuite',
                      clean      = 0,
                      debug      = 0)

from TestCases import TruchasCapabilityTest
import os

class GetAbortedResidualTest(TruchasCapabilityTest):

    def setDataStores(self):
        "defines all datastores to be used throughout this TestCase"

        #define testdata output directories..these result from running a BasicRunFail TestCase
        self.testdirs     = ['ht_bc_fail_output']

    def testPlotNLRAborts(self):
        "tests plotting non-linear residual aborts with GMV"

	mylogdir  = self.basicrunspecs.logdir

        self.nlrs = self.testdata[0].getAborts(type = 'NONLINEAR')
        for nlr in self.nlrs:
            idx      = self.nlrs.index(nlr)
            pltfile  = 'gmvnlr%d' %(idx)
	    gmvfile  = os.path.join(mylogdir,pltfile)
            self.plotter.plotFields(vis      = 'GMV',
                                    plotfile = gmvfile, 
                                    mesh     = self.testdata[0].meshes[0],
                                    fields   = nlr.vlist)
            
    def testPlotLRAborts(self):
        "tests plotting linear residual aborts with GMV"

	mylogdir  = self.basicrunspecs.logdir

        self.lrs  = self.testdata[0].getAborts(type = 'LINEAR')
        for lr in self.lrs:
            idx      = self.lrs.index(lr)
            pltfile  = 'gmvlr%d' %(idx)
	    gmvfile  = os.path.join(mylogdir,pltfile)
            self.plotter.plotFields(vis      = 'GMV',
                                    plotfile = gmvfile,
                                    mesh     = self.testdata[0].meshes[0],
                                    fields   = lr.vlist)

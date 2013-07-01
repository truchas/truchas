#!/usr/bin/env python
"""
 DataPlotter

-----------------------------------------------------------------------------
   Purpose:
  
      Provides convenient methods for ...
  
   Public Interface:
  
      T = DataPlotter(filename,mode,buffering)
      T.plotFields(field,region,desc)
        
   Contains:
      class DataLogger
        __init__(filename,mode,buffering)
        plotFields(field,region,desc)
        
   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
           Erin Iesulauro Barker (eibarker@lanl.gov)
-----------------------------------------------------------------------------
"""

import os, sys, string
try:
   import numpy.oldnumeric as Numeric
except ImportError:
   import Numeric
except:
   raise

if __name__ == '__main__':
    print "\n for component test in %s \n" %(__file__)
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from GMVwriteStorageObject     import GMVwriteStorageObject
from TECPLOTwriteStorateObject import TECPLOTwriteStorageObject
from VTKwriteStorageObject     import VTKwriteStorageObject

class DataPlotter:
    "class that allows developers to write fields in their TestCase to a viz file"

    def __init__(self,fpwatch=sys.stdout,debug=0):

        self.fpwatch = fpwatch
        self.debug   = debug

    def plotFields(self,vis,plotfile,mesh,fields,time=0.0,format='ascii'):

	if vis == 'GMV':
	    if format == 'binary':
		    file = plotfile + '.bin'
	    else:
		    file = plotfile + '.ascii'
	    writer = GMVwriteStorageObject(file, format, mesh, time, 0, fields,
                                           fpwatch=self.fpwatch, debug=self.debug)
	elif vis == 'TecPlot':
	    if format == 'binary':
		#write error message saying only ascii available
		errstr = 'TecPlot Binary is not available\n' + 'Getting TecPlot Ascii\n'
		self.fpwatch.write(errstr)
		format = 'ascii'
	    writer = TECPLOTwriteStorageObject(plotfile+'.ascii', format, mesh, time, 0, fields,
                                               fpwatch=self.fpwatch, debug=self.debug)
        elif vis == 'VTK':
	    file = plotfile +'.' + format + '.vtk'
	    writer = VTKwriteStorageObject(file, format, mesh, time, 0, fields,
                                           fpwatch=self.fpwatch, debug=self.debug)

if __name__=='__main__':
  
    from DataCreator import getDataCreator
    sys.path.append(thisdir + '/../')
    from TestCases   import CheckPPGMV, CheckPPTecPlot, CheckPPVTK
    import unittest
    from Runners     import TestRunner, RunTimeSpecs
    from PYTHONutils import uTestOpts
    import Clean

    testdir  = '../TestCases/ht_bc_fail_output'
    basedir  = '/../TestCases'
    test     = 'ht_bc_fail_output'
    testname = 'ht_bc_fail'
            
    # file, debug, output prefix, binary, ascii, clean
    opts = uTestOpts('fdco', 
                     defaults = {'d' : False,
                                 'c' : False,
				 'o' : 'outputs_',
				 'f' : testname+'.TBrook.xml'},
                     actions  = {'d' : 'store_true',
                                 'c' : 'store_true'},
                     dir      = testdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    try:

	if opt.c:
	    print "Cleaning up from a past test run...\n"
	    cmd = 'rm -r ' + opt.o + '*'
	    os.system(cmd)
	    Clean.Clean()

	else:
	    if opt.d: debug=1
	    else:     debug=0
    
            TestRunner = unittest.TextTestRunner

            errstr   = "\n The directory " + testdir + " does not exsist \n"
            errstr  += "Please prodcue the directory or alter the component test \n"
            assert test in os.listdir(os.path.join(os.getcwd()+basedir)), errstr
            
            testdata = getDataCreator(testdir)
            testdata.getStorage()
        
            #create dummy log dir and log file for purpose of running CheckPPGMV Testcase...
            logdir   = opt.o + __file__[0:-3]
            if not os.path.isdir(logdir):
               os.mkdir(logdir)
            logfile  = os.path.join(logdir,'ht_bc_fail_logs','CheckPPGMV.log')   
        
            fpwatch  = open('plotFields.log','w')
            plotter  = DataPlotter(fpwatch = fpwatch)
            nlrs     = testdata.getAborts(type = 'NONLINEAR')
            for nlr in nlrs:
        	idx     = nlrs.index(nlr)
                gmvfile = 'gmvnlr%d' % (idx)
                pltfile = os.path.join(logdir,gmvfile) 
                plotter.plotFields(vis       = 'GMV',
                                   plotfile  = pltfile,
                                   mesh      = testdata.meshes[0],
                                   fields    = nlr.vlist,
                                   time      = 0.0)
        
                tecfile = 'tecplotnlr%d' % (idx)
                pltfile = os.path.join(logdir,tecfile) 
                plotter.plotFields(vis       = 'TecPlot',
                                   plotfile  = pltfile,
                                   mesh      = testdata.meshes[0],
                                   fields    = nlr.vlist,
                                   time      = 0.0)
        
                vtkfile = 'vtknlr%d' % (idx)
                pltfile = os.path.join(logdir,vtkfile) 
                plotter.plotFields(vis       = 'VTK',
                                   plotfile  = pltfile,
                                   mesh      = testdata.meshes[0],
                                   fields    = nlr.vlist,
                                   time      = 0.0)
            fpwatch.close()
        
        
            LastRunSpecs = RunTimeSpecs(logdir,logfile) 
            suite        = unittest.TestSuite()
            T1           = CheckPPGMV(os.getcwd(),LastRunSpecs) 
            suite.addTest(T1)
            T2           = CheckPPTecPlot(os.getcwd(),LastRunSpecs) 
            suite.addTest(T2)
            T3           = CheckPPVTK(os.getcwd(),LastRunSpecs) 
            suite.addTest(T3)
            runner       = TestRunner(verbosity=2)
            result       = runner.run( suite )    
            print T1.str()
            print T2.str()
            print T3.str()
        
    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise

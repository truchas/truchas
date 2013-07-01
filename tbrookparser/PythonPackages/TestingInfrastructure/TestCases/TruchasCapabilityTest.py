#!/usr/bin/env python
"""
TruchasCapabilityTest

-----------------------------------------------------------------------------
   Purpose:
  
      TestCase to run a developer-written Truchas Capability TestCase
  
   Public Interface:
  
      T = TruchasCapabilityTest(unittest.TestCase,TruchasBaseTest)
      T.setDataStores()
      T.setDefinitions()
      T.setTolerances()
      T.setUp()
      T.shortDescription()
      T.runTest()
      T.tearDown()
      T.str()
  
   Contains:  
      class TruchasCapabilityTest
            setDataStores()
            setDefinitions()
            setTolerances()
            setUp()
            shortDescription()
            runTest()
            tearDown()
            str()
            __specstr(column)
            __tolstr(column1,column2)
  
   Unit Test Block
  
   Author(s): Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""
import unittest
import trutest
import os, sys, string, re, platform

thisdir     = os.path.abspath(os.path.dirname(__file__))
if __name__ == '__main__':
    print "\n for component test in %s \n" %(__file__)
    testingdir  = thisdir + '/../'
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)
    
from DataHandlers       import DataMeasure, DataLogger, getDataCreator, DataOperator, DataPlotter
from Runners            import TestRunner, RunTimeSpecs
from TruchasBaseTest    import TruchasBaseTest

class TruchasCapabilityTest(trutest.TestCase,TruchasBaseTest):

    def setDataStores(self):
        "defines datastores to be used throughout the Capability TestCase"

    def setDefinitions(self):
	"defines specific common inputs (such as regions) to be used throughout the Capability TestCase"

    def setTolerances(self):
	"developer inputs tolerances here"

    def setUp(self):
	"setups infrastructure needed for the test"
        
	#first establish if any input file(s), restart input file(s) and a golden output directory exists 
        self.tol             = {}
        self.testdirs        = []
        self.goldendirs      = []
        self.restartdirs     = []
	self.testdata        = [] # i.e list of testdata storages
	self.goldendata      = [] # i.e list of goldendata storages
	self.restartdata     = [] # i.e list of restartdata storages

	logdir               = self.basicrunspecs.logdir
	srcdir               = os.path.abspath(thisdir + '/../../../../regressiontests')
	self.dir             = os.path.abspath(logdir + '/../')
        L                    = string.split(self.dir,'/')
        self.caseName        = L[len(L)-1]        
        
	#get data store(s) 'testdata' created from the BasicRun TestCase(s)

        try:
            os.mkdir(logdir)
        except:
            "log dir already exists"

        self.setTolerances()
        self.setDataStores()

        self.testdirwarning      = "\n\n    WARNING!! No 'testdirs' list specified in setDefinitions!!   \n"
        self.testdirwarning     += "          You will not be able to utilize testdata in your script!! \n\n"

        self.goldendirwarning    = "\n\n    WARNING!! No 'goldendirs' list specified in setDefinitions!!   \n"
        self.goldendirwarning   += "          You will not be able to utilize goldendata in your script!!   \n\n"

        self.restartdirwarning   = "\n\n    WARNING!! No 'restartdirs' list specified in setDefinitions!!   \n"
        self.restartdirwarning  += "          You will not be able to utilize restartdata in your script!!   \n\n"
        
        #if len(self.testdirs)    == 0: print self.testdirwarning
        #if len(self.goldendirs)  == 0: print self.goldendirwarning
        #if len(self.restartdirs) == 0: print self.restartdirwarning

	#get data store(s) 'datastore' created from BasicRun TestCase (if truchas inpfile(s) exists)
        if len(self.datastores):
            for ds in self.datastores:
                for testdir in self.testdirs:
                    outdir     = os.path.join(logdir,testdir)
                    if os.path.abspath(ds.dir) == os.path.abspath(outdir):
                        self.testdata.append(ds)
        else:
            for testdir in self.testdirs:
                outdir  = os.path.abspath(os.path.join(logdir,testdir))
                T       = getDataCreator(dir=outdir)
                T.getStorage()
                self.testdata.append(T)            

	#get data store(s) 'restartdata' created from RestartRun TestCase (if restartinpfile(s) exists)
        if len(self.restartdatastores):
            for ds in self.restartdatastores:
                for restartdir in self.restartdirs:
                    outdir     = os.path.join(logdir,restartdir)                    
                    if os.path.abspath(ds.dir) == os.path.abspath(outdir):
                        self.restartdata.append(ds)
        else:
            for restartdir in self.restartdirs:
                outdir  = os.path.abspath(os.path.join(logdir,restartdir))
                T       = getDataCreator(dir=outdir)
                T.getStorage()
                self.restartdata.append(T)

        if len(self.goldendatastores):
            for gs in self.goldendatastores:
                gs.getStorage()
                self.goldendata.append(gs)

        else:
            for goldir in self.goldendirs:
                outdir      = os.path.abspath(os.path.join(srcdir,goldir))
                G           = getDataCreator(dir=outdir)
                G.getStorage()
                self.goldendata.append(G)

        self.setDefinitions()

        self.measure      = DataMeasure()
        self.operator     = DataOperator()

        self.plotfile     = '/dev/null'
	self.plotter      = DataPlotter(fpwatch=open(self.plotfile,'w'))

        self.logname      = 'Capability.log'
        self.logfile      = os.path.join(logdir,self.logname)

	#for formatting log file
	self.column1      = 0
	self.column2      = 5

        #if we have multiple testMethods within a TestCase then ensure we append to the log file....
        if self.testNumber == 1:
            self.logger    = DataLogger(filename=self.logfile,mode='w',buffering=1)
            self.logger.write('\n')
            self.logger.write('In Capability TestCase')
            self.logger.write('\n')
            self.logger.write('\n')
                
            tmp     = '*Setting up'
            str     = self.rjustln(self.column1, tmp)
            self.logger.write(str)
                
            tmp     = 'Current dir     : ' + self.dir 
            str     = self.rjustln(self.column2, tmp)
            self.logger.write(str)
            self.logger.write(self.__specstr(self.column2))
            self.logger.write('\n')
            self.logger.write(self.__tolstr(self.column2))
        else:
            self.logger    = DataLogger(filename=self.logfile,mode='a',buffering=1)

        self.logger.write('\n')
        tmp      = '*testMethod = %s' %(self.methodName)
        str      = self.rjustln(self.column1, tmp)
        self.logger.write(str)
        self.logger.write('\n')


    def fail(self, msg=None):
        "automatic failure; log a message to Capability.log file"

        mes  = "\n\nFAILURE in 'fail' criteria for testCase '%s' in testMethod '%s'!!! \n\n" %(self.caseName,self.methodName)
        tmp  = self.rjustln(self.column1, mes)
        if msg != None:
            mes += self.rjustln(self.column1, str(msg))
                
        self.logger.write(mes)
        self.logger.write('\n')

        raise self.failureException, mes


    def failIf(self, expr, msg=None):
        "simple check for failure if expr is true; log the result to Capability.log file"

        if expr:
            mes  = "\n\nFAILURE in 'failIf' criteria for testCase '%s' in testMethod '%s'!!! \n\n" %(self.caseName,self.methodName)
            tmp  = self.rjustln(self.column1, mes)
            if msg != None:
                mes += self.rjustln(self.column1, str(msg))

            self.logger.write(mes)
            self.logger.write('\n')

            raise self.failureException, mes

    def runTest(self):
	"basic capability testMethod - to be overloaded by developer"

	self.logger.write('\n')
	tmp     = '*testMethod = runTest'
	str     = self.rjustln(self.column1, tmp)
	self.logger.write(str)

	self.failIf(7>8)

    def shortDescription(self):
	"modified unittest short description to allow class name and input file name"

        classname = self.__class__.__name__
        classname = re.sub('Test[ers]*$', '', classname)
        docname   = '['+classname + '] ' + '(' + self.methodName + ') '
        return  docname.ljust(40) + ': ' + str(trutest.TestCase.shortDescription(self)) 


    def tearDown(self):
	"ensures working directory returned to original directory before test started"

	os.chdir(self.dir)

        if self.testNumber == self.numberOfTests:
            tmp     = '*Tearing down'
            str     = self.rjustln(self.column1, tmp)
            self.logger.write(str)

            tmp     = 'Returning to dir: ' + self.dir
            str     = self.rjustln(self.column2, tmp)
            self.logger.write(str)
            self.logger.write('\n')

            self.logger.close()

    def __specstr(self,column=15):
	"provides formatted info about postprocessor runtime specs for a given macro file"

        info = ''
        if (len(self.testdirs) > 0):
            tmp   = 'Testdata dirs   : ' + str(self.testdirs)
            info += self.rjustln(column, tmp)
        if (len(self.restartdirs) > 0):
            tmp   = 'Restartdata dirs: ' + str(self.restartdirs)
            info += self.rjustln(column, tmp)
        if (len(self.goldendirs) > 0):
            tmp   = 'Goldendata dirs : ' + str(self.goldendirs)
            info += self.rjustln(column, tmp)
        L  = string.split(self.logfile,'/')
        i  = 0
        for word in L:
            if 'outputs.t' in word: i=L.index(word)
        thelogfile = string.join(L[i:],'/')
        tmp   = 'Log file        : ' + thelogfile
        info += self.rjustln(column, tmp)
        info += '\n'
	tmp   = 'Basic RunSpecs: ' + self.basicrunspecs.str()
        info += self.rjustln(column, tmp)
	if self.restartrunspecs != None:
	    tmp   = 'Restart RunSpecs: ' + self.restartrunspecs.str()
            info += self.rjustln(column, tmp)

        return info

    def __tolstr(self,column=10,column2=38):
        "provides formatted info about developer chosen tolerances"

        tmp   = 'Defined Tolerances:'
        info  = self.rjustln(column, tmp)

        length = -1
        for tolerance in self.tol:
            length = max(length,len(str(tolerance)))
        for tolerance in self.tol:
            value  = '%10.5e' %(self.tol[tolerance])
            extra  = length - len(str(tolerance)) + 1
            tmp    = str(tolerance) + ' '*extra + value
            info  += self.rjustln(column, tmp)

	return info


if __name__=='__main__':
    currdir     = os.getcwd()
    TestRunner  = unittest.TextTestRunner

    from PYTHONutils import uTestOpts
    import Clean
    
    # file, debug, output prefix, binary, ascii, clean
    opts = uTestOpts('dc', 
                     defaults = {'d' : False,
                                 'c' : False},
                     actions  = {'d' : 'store_true',
                                 'c' : 'store_true'},
                     dir      = currdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,None)

    try:

        if opt.c: #clean

	    print "Cleaning up from a past test run...\n"
	    Clean.Clean()

	else:
    #determine runtime         directories
            output_dirs = []
            for dir in os.listdir(currdir):
        	if 'outputs.t-' in dir and os.path.isdir(dir): output_dirs.append(dir)
        	    
            # loop over available output directories
            for dir in output_dirs:
                logdir       = dir
                logfile      = logdir+'/static_drop_logs/BasicRun.log'
        
                # parse directory name to get runtime options
        	vals = string.splitfields(dir,'.')
        	if vals[1][2:] == string.lower(platform.system()):
                    parallel_env = vals[4] #'serial'
                    np           = 1
                    compiler     = vals[3] #'lahey'
        	    compile_mode = vals[5][:-2] #'opt'
                    LastRunSpecs = RunTimeSpecs(logdir,logfile,parallel_env,np,compiler,compile_mode)
        
                    suite        = unittest.TestSuite()
                    T            = TruchasCapabilityTest(datastores        = [],
                                                         restartdatastores = [],
                                                         goldendatastores  = [],
                                                         numberoftests     = 1,
                                                         testnumber        = 1,
                                                         BasicRunSpecs     = LastRunSpecs,
                                                         methodName        = 'runTest')  
                    suite.addTest(T)
                    runner       = TestRunner(stream=sys.stdout,verbosity=2)
                    result       = runner.run( suite )    
        
    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise
        
        
        
        
        
        

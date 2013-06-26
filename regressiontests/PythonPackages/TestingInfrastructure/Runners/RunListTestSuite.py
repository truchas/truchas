#!/usr/bin/env python
"""
 RunListTestSuite

-----------------------------------------------------------------------------
   Purpose:
  
      Creates a TestSuite from an inputted list.
      Each entry in the list invokes the hierarchical runner (RunHierarchicalTestSuite)
  
   Public Interface:

      T = RunListTestSuite(self,argv,currdir,outputdir,testlist)
      T.options(self)
      T.createOutputDir(self)
      T.listLoadAndRun(self)
      
   Contains:

      class RunListTestSuite(self,argv,currdir,outputdir,testlist)
            __init__(self,argv,currdir,outputdir,testlist)
            options(self)
            createOutputDir(self)
            listLoadAndRun(self)

   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import os, sys, string, fnmatch, shutil

thisdir     = os.path.abspath(os.path.dirname(__file__))
testingdir  = thisdir + '/../'

if __name__=='__main__':
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from RunBaseTestSuite         import RunBaseTestSuite
from RunHierarchicalTestSuite import RunHierarchicalTestSuite
from Loaders                  import TruchasCapabilityLoader,TruchasPostProcessorLoader

class RunListTestSuite(RunBaseTestSuite):
    """
    Creates a TestSuite from an inputted list.
    Each entry in the list invokes the hierarchical runner (RunHierarchicalTestSuite).
    This class interfaces with the Truchas GNUmakefile via the file
    tdir/runRegressionTests.py
    """

    def __init__(self,
                 argv,
                 bindir,
                 currdir,
                 outputdir,
                 testlist = ['postprocessor'],
                 debug    = 0):

        self.list         = testlist
        self.argv         = argv
        self.currdir      = currdir
        self.cwd          = self.currdir
        self.outputdir    = outputdir
        self.bindir       = bindir
        self.parallel_env = None
        self.np           = 1
        self.compile_mode = None
        self.compiler     = None
	self.version      = '2.3'
        self.executables  = []
        self.exevalue     = None
        self.no_rerun     = 0
        #for formatting output
        self.column       = 5
        self.column1      = 10
        self.column2      = 15
        self.debug        = debug

        self.setLogDir(outputdir = self.outputdir)

    def options(self):
        "additional -l (list) and -o (output dir) options over and above basic set of options"

        self.defineUsage()
        self.usage        = self.baseusage.replace('RunTestSuite.py','RunListTestSuite.py')        

        self.usage       += """
        python RunListTestSuite.py -l list
             - runs test suite contained in 'list'

        python RunListTestSuite.py -o outputdir
             - places all output from a testsuite run in 'outputdir'              
             
        """

        self.getOptions()
        self.mycommands.add_option("-l","--list", action="store",type="string",
                              help="LIST containing all TestCases to be run")
        self.mycommands.add_option("-o","--outputdir", action="store",type="string",
                                   help="Places all output from a testsuite run in OUTPUTDIR")
        
        self.parseOptions()
        self.checkOptions()

        if self.options.list != None:
            self.list = string.split(self.options.list,' ')

        #remove spurious spaces that may result from makefile list...
        for a in self.list:
            if not len(a):self.list.remove(a)

    def createOutputDir(self):
        "output directory to place testcase logs in"

        if self.options.outputdir != None:
            self.outputdir = self.options.outputdir
            
        otptdr = os.path.abspath(self.outputdir)
        try:
            os.chdir(otptdr)
        except:
            os.mkdir(otptdr)
            os.chdir(otptdr)

        #cmd = 'cp -r %s/* %s/.' %(self.currdir,otptdr)
        #os.system(cmd)

        self.__recursiveCopy(self.currdir,otptdr)

        self.outptdr = os.getcwd()
        
    def loadAndRun(self):

        fpsummary = sys.stdout

        fpsummary.write('\n')
        fpsummary.write('RESULTS DIR : %s' %(str(self.outptdr)))
        fpsummary.write('\n')

        for executable in self.executables:

            fpsummary.write('\n')
            desc = 'Executable : %s' %(str(executable))
            if self.np > 1 and 'parallel' in executable:
                desc = 'Executable : %s    np : %i' %(str(executable),self.np)
            fpsummary.write(desc)
            fpsummary.write('\n')
            fpsummary.write('-'*70)
            fpsummary.write('\n')

            totaltime = 0.0
            for dir in self.list:

                thisdir = os.path.abspath(os.path.join(self.outptdr,dir))
                srcdir  = os.path.abspath(os.path.join(self.currdir,dir))

                X       = RunHierarchicalTestSuite(self.argv,
                                                   srcdir  = srcdir,
                                                   currdir = thisdir,
                                                   debug   = self.debug)
                                                   
            
                X.setRunTimeEnvironment(self.no_rerun,
                                        self.np,
                                        self.compiler,
                                        self.parallel_env,
                                        self.compile_mode,
					self.version,
                                        self.bindir,
                                        self.executables)

                X.loadAndRun(loader     = None,
                             executable = executable,
                             fpsummary  = fpsummary,
                             verbtofile = 1)

                totaltime += X.timeTaken

        t = '%7.3f' %(totaltime)
        fpsummary.write('\nTime taken : %s secs \n\n' %t)
                
        fpsummary.close()

    def __recursiveCopy(self,currdir,otptdir): 
        
        print currdir
        print 'os.path.isdir=', os.path.isdir(currdir)
        if os.path.isdir(currdir):
            self.__myDirCopy(os.path.abspath(currdir),
                             os.path.abspath(otptdir))
            thedirs = os.listdir(currdir)
            for f in thedirs:
                self.__recursiveCopy(os.path.abspath(os.path.join(currdir,f)),
                                     os.path.abspath(os.path.join(otptdir,f)))
        else:
            return


    def __myDirCopy(self,currdir,otptdir,exclude_patterns=['*_golden','CVS*','testdoc']):

        files = os.listdir(currdir)
        for file in files:
            thismember = os.path.join(currdir,file)
            newmember  = os.path.join(otptdir,file)
            match      = 0
            for pattern in exclude_patterns:
                if fnmatch.fnmatch(file,pattern):
                    match += 1
            if not match and os.path.isdir(thismember) and not os.path.isdir(newmember):
                os.mkdir(newmember)
            if not match and os.path.isfile(thismember) and os.path.isdir(otptdir):
                shutil.copy(thismember,newmember)

if __name__=="__main__":
    argv     = sys.argv[1:]
    """
    testlist = ['postprocessor', 'em1', 'em2', 'em3', 'natural_convection', 'turbulence']
    testlist.extend(['void_collapse', 'porous_drag', 'channel_flow', 'static_drop','ht_bc'])
    testlist.extend(['species_advection', 'ds1', 'diverging_duct', 'conduction_45_so', 'conduction_90'])
    testlist.extend(['../tools/PythonPackages/TestingInfrastructure/examples'])
    """
    #testlist = ['../tools/PythonPackages/TestingInfrastructure/examples']
    #testlist  = ['postprocessor/restart/mapping', 'grain_growth']
    testlist  = ['postprocessor/vis']
    #testlist  = ['tm_pc_both']
    #testlist.extend(['viscoplastic_ring', 'contact_box_close'])
    #testlist.extend(['postprocessor', 'hierarchical_example', 'em1', 'em2', 'em3', 'natural_convection', 'turbulence'])
    #testlist.extend(['broken_dam', 'void_collapse', 'porous_drag', 'channel_flow', 'static_drop','ht_bc'])
    #testlist.extend(['species_advection', 'ds1', 'diverging_duct', 'conduction_45_so', 'conduction_90'])
    #testlist.extend(['phasechange_eutectic', 'phasechange_mixed', 'phasechange_CnK', 'sensitivity_stefan'])
    #testlist.sort()

    currdir  = os.getcwd() + '/../../../../regressiontests/'
    outputdir= os.getcwd() + '/../../../../regressiontests_output/'

    Y        = RunListTestSuite(argv,currdir,outputdir,testlist,debug=1)
    Y.options()
    Y.createOutputDir()
    Y.loadAndRun()


    



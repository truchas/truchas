#!/usr/bin/env python
"""
 TruchasTestLogger

-----------------------------------------------------------------------------
   Purpose:
  
      To provide ability to log status of TestSuite for debugging purposes.  
  
   Public Interface:
  
      T = TestLoggers(basename)
      T.getThem()
  
   Contains:
      class TestLoggers
             __init__(basename)
             getThem()
             
   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
           Erin Iesulauro Barker (eibarker@lanl.gov)
-----------------------------------------------------------------------------
"""
import os, sys 

if __name__=='__main__':
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    testingdir  = thisdir + '/../'
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

import logging

class TestLoggers:

    def __init__(self,basename,debugging=1,getinfo=0):

	self.basename  = basename
        self.debugging = debugging
        self.getinfo   = getinfo

    def getThem(self):

	#set up logging files...(.debug, .log)
	logger = logging.getLogger(self.basename)

        #create formatter
	formatter = logging.Formatter("Problem - %(name)s - %(levelname)s - %(message)s")

        if self.debugging:
            logger.setLevel(logging.DEBUG)
            filename = self.basename+'.debug'
            fh_debug = logging.FileHandler(filename,'w')
            fh_debug.setLevel(logging.DEBUG)
            #add formatter to fh_debug
            fh_debug.setFormatter(formatter)
            logger.addHandler(fh_debug)	

        if self.getinfo:
            filename = self.basename+'.log'
            fh_info  = logging.FileHandler(filename,'w')
            fh_info.setLevel(logging.INFO)
            logger.addHandler(fh_info)	
	
	return logger

if __name__ == '__main__':
    " TruchasTestLogger Component test"

    from PYTHONutils import uTestOpts
    
    # file, debug, output prefix, binary, ascii, clean
    dfltdir  = '../TestCases/'
    opts = uTestOpts('dc', 
                     defaults = {'d' : False,
                                 'c' : False},
                     actions  = {'d' : 'store_true',
                                 'c' : 'store_true'},
                     dir      = dfltdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,None)

    try:

        if opt.c: #clean

	    print "Cleaning up from a past test run...\n"
	    cmd = 'rm MyFile*'
	    os.system(cmd)

	else:
            X = TestLoggers('MyFile',debugging=1).getThem()
            X.debug('Testing logging levels - debug')
            X.info('Testing logging levels - info')

    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise

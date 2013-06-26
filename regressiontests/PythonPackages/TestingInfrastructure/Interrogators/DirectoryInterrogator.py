#!/usr/bin/env python
"""
 DirectoryInterrogator

-----------------------------------------------------------------------------
   Purpose:
  
      Provides methods to extract most previous:
      1) Truchas BasicRunSpecs
      2) Truchas RestartRunSpecs
  
      by searching through the outputs.* directories within the current directory
  
   Public Interface:
  
      T = DirectoryInterrogator(currdir)
      T.getBasicRunSpecs()
      T.getRestartRunSpecs()
      T.getLatest(logdir,pattern)
      T.extractSpecs(logdir)
  
   Contains:
      class DirectoryInterrogator
             __init__(currdir)
             getBasicRunSpecs()
             getRestartRunSpecs()
             getLatest(logdir,pattern)
             extractSpecs(logdir)

   Unit Test Block
   
   Author: Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import os, sys, string, fnmatch

thisdir     = os.path.abspath(os.path.dirname(__file__))
testingdir  = thisdir + '/../'

if __name__ == "__main__":
    sys.path.append(testingdir)
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from Runners     import RunTimeSpecs
from PYTHONutils import getAbortString

class DirectoryInterrogator:

    def __init__(self,currdir):

        self.currdir = currdir

        self.files         = os.listdir(self.currdir)
        self.logDirPat     = 'outputs.t*'
        self.basicRunPat   = 'BasicRun*.log'
        self.restartRunPat = 'RestartRun*.log'
        
        self.abortstring  = getAbortString()

        self.nologdirerr  = '\n \n                        No outputs.* directory exists in \n '
        self.nologdirerr += '     %s \n ' %(self.currdir)  
        self.nologdirerr += '              Cannot use the no-rerun option when running your Capability TestCase \n'
        self.nologdirerr += self.abortstring
        

        self.nobasicrunerr = '\n \n                No BasicRun.log file exists in any outputs.* directory in %s \n ' %(self.currdir)  
        self.nobasicrunerr += '                    Cannot use the no-rerun option when running your Capability TestCase \n'
        self.nobasicrunerr += self.abortstring

        self.norestartrunerr = '\n \n                No RestartRun.log file exists in any outputs.* directory in %s \n ' %(self.currdir)  
        self.norestartrunerr += '                    Cannot use the no-rerun option when running your Capability TestCase \n'
        self.norestartrunerr += self.abortstring

        
    def getBasicRunSpecs(self):

        latestbasicrundate = '-2000'
        logdir             = None
        runfile            = None
        rundate            = '-1000'
        
        for f in self.files:
            if fnmatch.fnmatch(f, self.logDirPat):
                rundate, runfile  = self.getLatest(f, self.basicRunPat)
                logdir            = f
                if runfile != None:
                    if  rundate > latestbasicrundate:
                        latestbasicrundate = basicrundate
                        logdir             = f
                    
        assert logdir  != None, self.nologdirerr
        assert runfile != None, self.nobasicrunerr
            
        BasicRunSpecs   = self.extractSpecs(self.currdir+'/'+logdir)

	return BasicRunSpecs


    def getRestartRunSpecs(self):

        latestrestartrundate = '-2000'
        logdir             = None
        runfile            = None
        rundate            = '-1000'
        
        for f in self.files:
            if fnmatch.fnmatch(f, self.logDirPat):
                rundate, runfile  = self.getLatest(f, self.restartRunPat)
                logdir            = f
                if runfile != None:
                    if  rundate > latestrestartrundate:
                        latestrestartrundate = restartrundate
                        logdir             = f
                    
        assert logdir  != None, self.nologdirerr
        assert runfile != None, self.norestartrunerr
            
        RestartRunSpecs = self.extractSpecs(logdir)

	return RestartRunSpecs


    def getLatest(self,logdir,pattern):

        thisdir = os.path.abspath(self.currdir+'/'+logdir)
        
        datemax = -10000
        filemax = None
        for logsfile in os.listdir(thisdir):
            if fnmatch.fnmatch(logsfile,'*_logs'):
                logsdir = thisdir + '/' + logsfile
                for file in os.listdir(logsdir):
                    if fnmatch.fnmatch(file,pattern):
                        date = os.stat(logsdir + '/' + file)[7]
                        if date > datemax:
                            datemax = date
                            filemax = file

        return datemax, filemax


    def extractSpecs(self,logdir):

        specs   = string.split(logdir,'.')
        npspecs = string.split(logdir,'_')

        np     = 0
        if len(npspecs) > 1:
            np = npspecs[1][2:]
            
        mode   = specs[5][0:3]

        runSpecs = RunTimeSpecs(logdir,
                                logdir+'/Capability.log',
                                specs[4],
                                np,
                                specs[3],
                                mode,
                                )

        return runSpecs


if __name__ == "__main__":

    currdir = os.getcwd()
    X               = DirectoryInterrogator(currdir)
    BasicRunSpecs   = X.getBasicRunSpecs()
    print BasicRunSpecs.str()
    RestartRunSpecs = X.getRestartRunSpecs()
    print RestartRunSpecs.str()    

        

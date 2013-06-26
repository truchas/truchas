#!/usr/bin/env python
"""
TruchasBaseTest

-----------------------------------------------------------------------------
   Purpose:
  
      Provides base utilities for the following TestCases
      TruchasBasicRunTest
      TruchasRestartRunTest
      TruchasPostProcessorRunTest      
      TruchasCapabilityTest
  
   Public Interface:
  
      T = TruchasBaseTest
      T.constructTruchasExe()
      T.deconstructTruchasExe()
      T.constructLogDirName()
      T.constructCompilerErrorMessage()
      T.constructPlatformErrorMessage()
      T.rjustln(column)
      T.createStatistics()

   Contains:  
      class TruchasBaseTest
            constructTruchasExe()
            deconstructTruchasExe()
            constructLogDirName()
            constructCompilerErrorMessage()
            constructPlatformErrorMessage()
            rjustln(column)
            createStatistics()
  
   Unit Test Block
  
   Author(s): Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import string, os, sys, commands

thisdir     = os.path.abspath(os.path.dirname(__file__))
testingdir  = thisdir + '/../'
sys.path.append(testingdir)
parserdir   = thisdir + '/../../TBrookParser'
sys.path.append(parserdir)

from writeStorageObject import writeStorageObject
from DataHandlers       import getDataCreator

class TruchasBaseTest:

    def deconstructTruchasExe(self):
        "given executable..gets parallel_env, compile_mode, compiler"

        tmp               = string.split(self.truchas_exe,'.')
        self.compiler     = tmp[2]
        self.parallel_env = tmp[3]
        tmp2              = string.split(tmp[4],'-')
        self.compile_mode = tmp2[0]
	#BROKENself.version      = tmp2[1] + '.' + tmp[5]
	self.version      = '2' + '.' + '8'


    def constructTruchasExe(self):
	"constructs the Truchas executable name"

	pltfrm = sys.platform
	if 'linux' in pltfrm:
	    self.platform  = 'linux'
	    self.arch      = commands.getstatusoutput('uname -i')[1]
	    comperr        = self.constructCompilerErrorMessage()
	    assert self.compiler != 'compaq', comperr
	elif 'osf1' in pltfrm:
	    self.platform  = 'osf1'
            tmp            = commands.getstatusoutput('uname -a')[1]
            self.arch      = string.split(tmp,' ')[4]
	    comperr        = self.constructCompilerErrorMessage()
	    assert self.compiler != 'nag', comperr
	    assert self.compiler != 'lahey', comperr
	else:
	    pltfrmerr      = self.constructPlatformErrorMessage()
	    raise pltfrmerr

	truchas_exe  = 't-'+self.platform+'.'+self.arch+'.'+self.compiler+'.'
	truchas_exe += self.parallel_env+'.'+self.compile_mode+'-' + self.version + '.dev'

	return truchas_exe

    def constructLogDirName(self):
	"constructs directory name to place log files in"

        logdir   = 'outputs.'
        logdir  += self.truchas_exe
	if self.parallel_env == 'parallel':
	    logdir +='_np'+str(self.np)

	return logdir

    def constructCompilerErrorMessage(self):
	"provides error message if wrong compiler chosen for given platform"

        Dcompilers = {'linux':'nag or lahey','osf1':'compaq'}

	comperr  ='\n \n Error! Chosen %s compiler on %s platform!!' %(self.compiler,self.platform)
	comperr +='\n Try compiler = %s.' %(Dcompilers[self.platform]) 
	comperr +='\n ABORTING %s. \n' %(self.__class__.__name__)

	return comperr

    def constructPlatformErrorMessage(self):
	"provides error message if running tests on unsupported platform"

	pltfrmerr  = '\n \n Error! Running on unknown platform!!'
	pltfrmerr +='\n Platforms must be either %s or %s \n' %('osf1','linux')
	pltfrmerr +='\n ABORTING %s. \n' %(self.__class__.__name__)

	return pltfrmerr


    def rjustln(self,column,s):

        tmp = s.rjust(column + len(s))
        tmp += '\n'

        return tmp

    def createStatistics(self,column=15):
        "creates cycle, time, variable statistics associated with this Truchas Run" 

        self.datastore.getStorage()
        stats = writeStorageObject().writeStatistics(self.datastore.storage,column)

        return stats




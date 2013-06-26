#!/usr/bin/env python

"""
 runRegressionTests

-------------------------------------------------------------------------------
   Purpose: Interface between GNUmakefile and Python testing infrastructure

   Author: Sharen Cummins (scummins@lanl.gov)

-------------------------------------------------------------------------------
""" 

import os, sys, string

# add directories to sys.path - directories are searched for importing modules

#rootdir    = os.path.abspath(os.path.dirname(__file__))
#testingdir = rootdir + '/tools/PythonPackages/TestingInfrastructure'
#runnersdir = rootdir + '/tools/PythonPackages/TestingInfrastructure/Runners'
#parserdir  = rootdir + '/tools/PythonPackages/TBrookParser'
#bindir     = rootdir + '/bin'

# Variable replacement from CMake
rootdir    = '@Truchas_SOURCE_DIR@'
testingdir = '@TruchasRegTests_SOURCE_DIR@' + '/PythonPackages/TestingInfrastructure'
runnersdir = '@TruchasRegTests_SOURCE_DIR@' + '/PythonPackages/TestingInfrastructure/Runners'
parserdir  = '@TruchasRegTests_SOURCE_DIR@' + '/PythonPackages/TBrookParser'
bindir     = '@TruchasExe_BINARY_DIR@'

sys.path.append(testingdir)
sys.path.append(runnersdir)
sys.path.append(parserdir)
sys.path.append(bindir)

from RunListTestSuite import RunListTestSuite

defcurrdir     = rootdir + '/regressiontests/'
defoutputdir   = rootdir + '/regressiontests_output/'
#defcurrdir     = '@TruchasRegTests_BINARY_DIR@'
#defoutputdir   = defcurrdir + '/regressiontests_output/'

argv        = sys.argv[1:]

#begin interfacing Makefile PARALLEL, DEBUG strings with our options...
parallel    = 'PARALLEL=yes'
serial      = 'PARALLEL=no'
debug       = 'DEBUG=yes'
opt         = 'DEBUG=no'

if parallel in argv:
    argv[argv.index(parallel)]  = '-m'
    if debug in argv:
        argv[argv.index(debug)] = 'parallel-dbg'
    if opt in argv:
        argv[argv.index(opt)]   = 'parallel-opt'
    
if serial in argv:
    argv[argv.index(serial)]    = '-m'
    if debug in argv:
        argv[argv.index(debug)] = 'serial-dbg'
    if opt in argv:
        argv[argv.index(opt)]   = 'serial-opt'

np = 0
#check for number of processors
for i in range(len(argv)-1,0,-1):
    if 'NUMPROCS=' in argv[i]:
        np = string.split(argv[i],'=')[-1]
        argv.remove(argv[i])
if np > 0:
    argv.append('--np')
    argv.append(np)

sys.argv[1:] = argv
Y            = RunListTestSuite(argv,bindir,defcurrdir,defoutputdir)
Y.options()
Y.createOutputDir()
Y.loadAndRun()



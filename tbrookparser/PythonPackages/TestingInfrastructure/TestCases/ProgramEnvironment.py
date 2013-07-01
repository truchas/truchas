#!/usr/bin/env python
"""
ProgramEnvironment

-----------------------------------------------------------------------------
   Purpose:
  
      To create the serial or parallel executable commands for the TestCases
      TruchasBasicRun
      TruchasRestartRun
      TruchasPostProcessorRun
  
   Public Interface:
  
      T = createRunTimeProgram(progname,args,outfile)
      T.run()
      T.str()
  
   Contains:  
      class __Program
            __init__(self, name, args, outfile)
  
      class __ProgramEnvironment
            run(self)
            str(self,column_indent)

      class __SerialEnvironment(__ProgramEnvironment)
            run(self)
      class __MPIEnvironment(__ProgramEnvironment)
            run(self)
  
   Unit Test Block
  
   Author(s): Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------
"""

import os, sys, string, subprocess ,popen2

class __Program:
    "specifies the name, arguments and outfile when a particular program is run"

    def __init__(self, name, args, outfile):

	self.name     = name
	self.args     = args
	self.outfile  = outfile

class __ProgramEnvironment:
    "defines the run time environment for any program (i.e Truchas, Chaparral etc)"

    def run(self):
	pass

    def str(self, column_indent=10):
	n    = string.rfind(self.cmd,'/')
	tmp  = self.cmd[n+1:]
	info ='Running %s' %(self.cmd)
	str  = info.rjust(column_indent + len(info))
	str +='\n'
        info ='Creating %s' %(self.program.outfile)
	str += info.rjust(column_indent + len(info))
	str +='\n'
        str +='\n'
	return str

    def teardown(self):
        pass

class __SerialEnvironment(__ProgramEnvironment):

    def __init__(self,program):

	self.program = program
	self.cmd     = ''
	self.bpsh    = ''
	if 'BEOWULF_JOB_MAP' in os.environ.keys():
	    nodes = string.splitfields(os.environ['NODES'],',')
	    self.bpsh = 'bpsh ' + nodes[0] + ' '

    def run(self):
	"call sub program and return exit status"

	self.cmd = self.bpsh + self.program.name + ' '
	for a in self.program.args:
	    self.cmd += (a+' ')
        #P = popen2.Popen4(self.cmd)
        #P.tochild.close()
        #file(self.program.outfile, 'w').write(P.fromchild.read())

	P = subprocess.Popen(self.cmd, stdout=subprocess.PIPE,
			    stderr=subprocess.PIPE, close_fds=True,
		            shell=True)
	(stdout,stderr) = P.communicate()
	file(self.program.outfile, 'w').write(stdout)
	
	return P.wait()

class __MPIEnvironment(__ProgramEnvironment):    
 
    def __init__(self, program, MPI_command):
        
	self.program    = program
	self.cmd        = ''
        self.mpicommand = MPI_command
	
    def run(self):
	"call sub program and return exit status"

	# could be replaced by subprocess for python >= 2.4
	self.cmd = self.mpicommand + self.program.args[0] + self.program.name + ' '
        for a in self.program.args[1:]:
	    self.cmd += (a+' ')

        #P = popen2.Popen4(self.cmd)
        #P.tochild.close()
        #file(self.program.outfile, 'w').write(P.fromchild.read())

	P = subprocess.Popen(self.cmd, stdout=subprocess.PIPE,
			    stderr=subprocess.PIPE, close_fds=True,
		            shell=True)
	(stdout,stderr) = P.communicate()
	file(self.program.outfile, 'w').write(stdout)

	return P.wait() 

    
    def teardown(self):
	pass


def createRunTimeProgram(program_name, args, outfile, maxnodes = None):

    ThisProg = __Program(program_name, args, outfile)

    if '-np' in args[0]:
        MPI_command = 'mpirun '
        try:
            if 'qsc' in os.environ['HOSTNAME']:
                MPI_command      = 'prun '
                tmp              = args[0].replace('np','n')
                ThisProg.args[0] = tmp
        except:
            'HOSTNAME check not relevant'
            
        return __MPIEnvironment(ThisProg, MPI_command)
    else:
        return __SerialEnvironment(ThisProg)
    

if __name__=="__main__":
    tdir           = os.path.abspath(os.path.dirname(__file__)+'../../../../bin')
    truchas_exe    = os.path.join(tdir,'truchas')
    S              = createRunTimeProgram(truchas_exe, ['-np 2 ', '-o:static_drop_out', 'static_drop.inp'], 'static_drop.out')
    retcode        = S.run()
    info           = S.str()
    S.teardown()
    
    print
    print info
    print '          retcode : ' + str(retcode)
    print


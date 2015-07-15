#!/usr/bin/env python

import sys
import os

import unittest

import Truchas
import TruchasTest

class RestartTest(TruchasTest.BaseTestCase):


  
  '''
  Runs before each test_* need to call the setUpClass
  to initialize the truchas instance correctly.
  '''
  def setUp(self):
    # Call the class setup function setUpClass
    if self._is_initialized is False:
      self.setUpClass()
    # Input file used in all tests  
    self.truchas.input = self.get_input_rootdir('restart') + \
	                 os.path.sep + 'ds11.inp'
    # Default number of procs is 4			 
    if self.truchas_is_parallel:
      self.num_procs = 4
      self.truchas.nprocs=self.num_procs
    self.test_name='restart'  

  def test_read_restart(self):
    '''Test creating and reading a restart file'''
    
    # Output location
    outdir='read_test'
    full_outdir=self.build_output_directory(self.test_name,outdir)
    self.clean_output_directory(self.test_name,outdir)

    # Define the outfiles
    self.truchas.outdir = full_outdir
    self.truchas.h5file = self.build_output_filename(self.test_name,outdir,'ds11.h5')
    self.truchas.stdout = self.build_output_filename(self.test_name,outdir,'ds11-first.tty')
    self.truchas.stderr = self.build_output_filename(self.test_name,outdir,'ds11-first.err')

    # First run 
    # Depending on the order of the tests restart != None
    # which will cause a fail. 
    self.truchas.restart=None 
    self.truchas.run()

    # Output location
    outdir='read_test_restart'
    full_outdir=self.build_output_directory(self.test_name,outdir)
    self.clean_output_directory(self.test_name,outdir)

    # Write a restart file
    restart_cycle=39
    restart_file=self.build_output_filename(self.test_name,outdir,'restart.%i'%restart_cycle)
    self.truchas.write_restart(restart_cycle,output=restart_file)

    # Define new outfiles in the restart directory
    self.truchas.outdir = full_outdir
    self.truchas.h5file = self.build_output_filename(self.test_name,outdir,'ds11-next.h5')
    self.truchas.stdout = self.build_output_filename(self.test_name,outdir,'ds11-next.tty')
    self.truchas.stderr = self.build_output_filename(self.test_name,outdir,'ds11-next.err')

    # Now run with the restart file
    self.truchas.run()

  def test_fields(self):
    '''Test comparing final field datasets when restarting'''

    # First run
    outdir1='out1'
    full_outdir1=self.build_output_directory(self.test_name,outdir1)
    self.clean_output_directory(self.test_name,outdir1)

    # Define the output files
    self.truchas.outdir=full_outdir1
    h5file1=self.build_output_filename(self.test_name,outdir1,'ds11.h5')
    self.truchas.h5file = h5file1
    self.truchas.stdout = self.build_output_filename(self.test_name,outdir1,'ds11.tty')
    self.truchas.stderr = self.build_output_filename(self.test_name,outdir1,'ds11.err')

    # First run
    # No restart
    self.truchas.restart=None
    self.truchas.run()

    # Second run, direct output to another directory 
    outdir2='out2'
    full_outdir=self.build_output_directory(self.test_name,outdir2)
    self.clean_output_directory(self.test_name,outdir2)

    # Write a restart file
    restart_cycle=39
    restart_file=self.build_output_filename(self.test_name,outdir2,'restart.%i'%restart_cycle)
    try:
      self.truchas.write_restart(restart_cycle,output=restart_file)
    except:
      raise TruchasError('Failed to create a restart file')

    # Define the output files
    self.truchas.outdir=full_outdir
    h5file2=self.build_output_filename(self.test_name,outdir2,'ds11.h5')
    self.truchas.h5file = h5file2
    self.truchas.stdout = self.build_output_filename(self.test_name,outdir2,'ds11.tty')
    self.truchas.stderr = self.build_output_filename(self.test_name,outdir2,'ds11.err')

    # Second run
    self.truchas.run()

    # Compare the final output 
    h51 = Truchas.TruchasOutput(h5file1)
    sim1 = h51.get_simulation()
    h52 = Truchas.TruchasOutput(h5file2)
    sim2 = h52.get_simulation()

    num_series1 = len(sim1.series_names)
    num_series2 = len(sim2.series_names)

    series_last1 = sim1.get_series(id=num_series1)
    series_last2 = sim2.get_series(id=num_series2)

    # Check all the dataset fields
    for ds in series_last1.datasets:
      if ds not in series_last2.datasets:
	self.fail('Missing dataset')
      d1 = series_last1.get_data(ds)
      d2 = series_last2.get_data(ds)
      if d1.all() != d2.all():
	self.fail('Mismatched dataset values')
    

  def test_repartition(self):
    '''Test restarting on a different number of processors'''

    # Use half the numprocs defined in other tests
    try:
      np1 = max(self.num_procs/2,1)
    except:
      np1 = 1
    self.truchas.nprocs=np1
    
    # Output location
    outdir='repartition_np%i'%np1
    full_outdir=self.build_output_directory(self.test_name,outdir)
    self.clean_output_directory(self.test_name,outdir)

    # Define the outfiles
    self.truchas.outdir = full_outdir
    self.truchas.h5file = self.build_output_filename(self.test_name,outdir,'ds11.h5')
    self.truchas.stdout = self.build_output_filename(self.test_name,outdir,'ds11.tty')
    self.truchas.stderr = self.build_output_filename(self.test_name,outdir,'ds11.err')

    # First run 
    # Depending on the order of the tests restart != None
    # which will cause a fail. 
    self.truchas.restart=None 
    self.truchas.run()

    # Use half the numprocs defined in other tests
    try:
      np2 = self.num_procs
    except:
      np2 = 1
    self.truchas.nprocs=np2

    # Output location
    outdir='repartition_restart_np%i'%np2
    full_outdir=self.build_output_directory(self.test_name,outdir)
    self.clean_output_directory(self.test_name,outdir)

    # Write a restart file
    restart_cycle=39
    restart_file=self.build_output_filename(self.test_name,outdir,'restart.%i'%restart_cycle)
    self.truchas.write_restart(restart_cycle,output=restart_file)

    # Define new outfiles in the restart directory
    self.truchas.outdir = full_outdir
    self.truchas.h5file = self.build_output_filename(self.test_name,outdir,'ds11.h5')
    self.truchas.stdout = self.build_output_filename(self.test_name,outdir,'ds11.tty')
    self.truchas.stderr = self.build_output_filename(self.test_name,outdir,'ds11.err')

    # Now run with the restart file
    self.truchas.run()


if __name__ == '__main__':
  unittest.main()
 
  

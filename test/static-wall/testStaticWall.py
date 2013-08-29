#!/usr/bin/env python

import sys
import os

import unittest
import numpy

import Truchas
import TruchasTest

class StaticWallTest(TruchasTest.BaseTestCase):

  def runTest(self):
    '''Test static wall zero velocity'''

    test_name='static-wall'

    tol=1.0e-11

    # Same input file, notice pointing to the restart directory
    self.truchas.input=os.path.join(self.get_input_rootdir(test_name),'static-wall.inp')

    # Root directory for output
    output_rootdir=self.get_output_rootdir('static-wall')
    self.truchas.outdir=output_rootdir

    # Setup for parallel test
    if self.truchas_is_parallel:
      self.truchas.nprocs=4

    # Run
    self.truchas.run()

    # Grab the velocity field
    field='Z_VC'
    h5_file=os.path.join(output_rootdir,'static-wall.h5')
    o=Truchas.TruchasOutput(h5_file)

    v=o.get_simulation().get_series_data(field,cycle=1)
    (nc,d)=v.shape
    v_norm=numpy.ones(nc,dtype=numpy.float64)
    for c in range(nc):
      v_norm[c]=numpy.linalg.norm(v[c,:],2)
    v_max=numpy.linalg.norm(v_norm,numpy.inf)  

    self.assertTrue(v_max<tol)

    

if __name__ == '__main__':
  unittest.main()




    





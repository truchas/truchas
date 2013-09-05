################################################################################

# --- Python NumPy module
import numpy

# --- Truchas 
import Truchas
import TruchasTest


class GapRadTest(TruchasTest.GoldenTestCase):

  test_name='gap-rad'

  num_procs=4


  probe_field='TEMP'
  probe_names=['left end', 'right end', 'gap left', 'gap right']
  probe={}
  probe['left end']   = [1.4999998474, 1.371828885]
  probe['right end']  = [0.500001526, 0.628171115 ]
  probe['gap left']   = [1.498533630, 1.370837833]
  probe['gap right']  = [0.5014663700, 0.629162167]

  def setUp(self):
    if self._is_initialized is False:
      self.setUpClass()
    self.test_output=Truchas.TruchasOutput(self.get_output_file())  

  def runTest(self):
    '''Test probe temperatures at left,right ends and gap'''

    fails=0

    # Find the last cycle
    probe=self.test_output.get_simulation().get_probe(self.probe_names[0],field=self.probe_field)
    data=probe.get_data()
    (n,d)=data.shape
    last_cyc=data[n-1,probe.CYCLE_INDEX]
    del data
    

    # Grab the initial probe data and compare
    tol=1.0e-9
    cyc=0
    for pname in self.probe_names:
      probe=self.test_output.get_simulation().get_probe(name=pname,field=self.probe_field)
      probe_data=probe.find_cycle_data(cycle=cyc)[0] # Returns an array slice, need a number to compare
      error=abs(self.probe[pname][0]-probe_data)
      if error > tol:
	print '(test) probe=%1.9e'%(self.probe[pname][0])
	print '(output) probe_data=%1.9e'%(probe_data)
	print 'Error exceeds tolerance %s cyc=%d %e %e\n'%(pname,cyc,error,tol) 
	fails+=1

    # Grab the final time
    tol=5.0e-5
    for pname in self.probe_names:
      probe=self.test_output.get_simulation().get_probe(name=pname,field=self.probe_field)
      probe_data=probe.find_cycle_data(cycle=last_cyc)[0] # Returns an array slice
      error=abs(self.probe[pname][1]-probe_data)
      if error > tol:
	print '(test) probe=%1.9e'%(self.probe[pname][1])
	print '(output) probe_data=%1.9e'%(probe_data)
	print 'Error exceeds tolerance %s cyc=%d %e %e\n'%(pname,cyc,error,tol) 
	fails+=1

    self.assertTrue(fails == 0)


if __name__ == '__main__':
  import unittest
  unittest.main()

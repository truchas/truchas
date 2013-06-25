# ############################################################################ #
#
# Python UnitTest for Outputs
#
# Really this is unit test of the H5Object, Can not create and  H5Object
# without generating an HDF5 file.
# 
# ############################################################################ #

import os, sys
from types import *

import unittest

from Danu import Output

class TestOutputAccess(unittest.TestCase):

  def setUp(self):
    self.filename         = 'test-Output.h5'

    # Create the file
    fh = Output(self.filename)
    fh.close()

  def tearDown(self):
    if os.path.exists(self.filename):
      os.remove(self.filename)

  def test_open_dflt(self):
    file = Output(self.filename)
    del file

  def test_open_rdonly(self):
    from Danu import FILE_ACC_RDONLY
    bad_file='does-not-exist-'+self.filename
    try:
      file = Output(bad_file,'r')
    except:
      print 'Caught the bad file error'
    else:
      print 'Failed to catch the error'
      raise
    file = Output(self.filename,'r')
    access = file.access
    self.assertEqual(access,FILE_ACC_RDONLY)
    del file

  def test_open_rdwr(self):
    from Danu import FILE_ACC_RDWR
    file = Output(self.filename,'w')
    access = file.access
    self.assertEqual(access,FILE_ACC_RDWR)
    del file

  def test_open_append(self):
    from Danu import FILE_ACC_APPEND
    file = Output(self.filename,'a')
    access = file.access
    self.assertEqual(access,FILE_ACC_APPEND)
    del file

  def test_flush(self):
    file = Output(self.filename,'a')
    file.flush()
    file.flush(1)
    del file

class TestOutputWriteAttribute(unittest.TestCase):

  def setUp(self):
    self.filename         = 'test-Output.h5'
    self.int_attr_name    = 'Dummy Int Attribute'
    self.double_attr_name = 'Dummy Double Attribute'
    self.str_attr_name    = 'Dummy String Attribute'
    self.fh               = Output(self.filename,'a')

  def tearDown(self):
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def test_string_attribute(self):
    import random
    str_len = random.randint(16,1024)
    self.fh.set_attribute(self.str_attr_name,os.urandom(str_len))
  
  def test_double_attribute(self):
    import random
    self.fh.set_attribute(self.double_attr_name,random.random())
  
  def test_int_attribute(self):
    import random
    self.fh.set_attribute(self.int_attr_name,random.randint(0,100000))

class TestOutputReadAttribute(unittest.TestCase):

  def setUp(self):
    import random
    import string

    self.filename         = 'test-Output.h5'

    self.fh               = Output(self.filename,'w')

    self.attributes       = []

    self.int_attr_name    = 'Dummy Int Attribute'
    self.attributes.append(self.int_attr_name)
    self.int_attr         = random.randint(0,1000000)
    self.fh.set_attribute(self.int_attr_name,self.int_attr)

    self.double_attr_name = 'Dummy Double Attribute'
    self.attributes.append(self.double_attr_name)
    self.double_attr      = random.random()
    self.fh.set_attribute(self.double_attr_name,self.double_attr)

    self.str_attr_name    = 'Dummy String Attribute'
    self.attributes.append(self.str_attr_name)
    str_len = random.randint(10,1024)
    self.str_attr         = ''.join(random.choice(string.letters) for i in xrange(str_len))
    self.fh.set_attribute(self.str_attr_name,self.str_attr)

  def tearDown(self):
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def test_string_attribute(self):
    test_str = self.fh.get_attribute(self.str_attr_name)
    self.assertEqual(test_str,self.str_attr)
   
  def test_double_attribute(self):
    test_double = self.fh.get_attribute(self.double_attr_name)
    self.assertEqual(test_double,self.double_attr)
  
  def test_int_attribute(self):
    test_int = self.fh.get_attribute(self.int_attr_name)
    self.assertEqual(test_int,self.int_attr)

  def test_attributes(self):
    attributes = self.fh.attributes()
    for attr in self.attributes:
      if attr not in attributes: TestCase.fail('Failed to read attributes correctly')


class TestOutputMesh(unittest.TestCase):
  
  def setUp(self):
    from Danu import LINE_ELEM, TRI_ELEM, QUAD_ELEM, TET_ELEM, HEX_ELEM

    self.filename = 'test-Output.h5'
    self.fh       = Output(self.filename,'w')

    self.mesh_count = 0
    self.mesh_names = []

    self.valid_elems = [LINE_ELEM, TRI_ELEM, QUAD_ELEM, TET_ELEM, HEX_ELEM]

  def tearDown(self):
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def get_random_elem(self):
    from random import choice
    return choice(self.valid_elems)

  def test_add_mesh_u(self):
    # Require and element type
    try:
      self.fh.add_unstruct_mesh('Test Mesh') 
    except RuntimeError:
      print 'Caught the invalid mesh value'
    # Test all avail mesh types  
    for elem in self.valid_elems:
      mesh_name = 'Test Mesh ' + str(elem)
      self.fh.add_unstruct_mesh(mesh_name,elem)
    # Try to add an existing mesh
    elem = self.valid_elems[0]
    mesh_name = 'Test Mesh ' + str(elem)
    try:
      self.fh.add_unstruct_mesh(mesh_name,elem)
    except:
      print 'Caught the exception when mesh already exists'

  def test_mesh_exists(self):
    import random
    import string
    mesh_name = ''.join(random.choice(string.letters) for i in xrange(16))
    self.assertEqual(self.fh.mesh_exists(mesh_name),0,'Failed to query DNE mesh correctly')
    mesh_name = 'Test Mesh'
    elem = self.get_random_elem()
    mesh = self.fh.add_unstruct_mesh(mesh_name,elem)
    self.assertFalse(mesh is NoneType, 'Failed to add random mesh type')
    self.assertTrue(self.fh.mesh_exists(mesh_name), 'Failed to locate existing mesh')

  def test_mesh_count(self):
    import random
    self.assertEqual(self.mesh_count,self.fh.mesh_count(), 'Failed to return correct mesh count')
    num_mesh = random.randint(1,32)
    m = 1
    while m <= num_mesh:
      elem = self.get_random_elem()
      mesh_name = 'Test Mesh ' + str(m)
      self.fh.add_unstruct_mesh(mesh_name,elem)
      self.mesh_count = self.mesh_count + 1
      m = m + 1
    self.assertEqual(self.mesh_count,self.fh.mesh_count(), 'Failed to return correct mesh count')

  def test_mesh_list(self):
    import random
    mesh_list = self.fh.mesh_list()
    self.assertEqual(len(mesh_list),self.mesh_count, 'Incorrect mesh list length')
    num_mesh = random.randint(1,32)
    m = 1
    while m <= num_mesh:
      elem = self.get_random_elem()
      mesh_name = 'Test Mesh ' + str(m)
      self.fh.add_unstruct_mesh(mesh_name,elem)
      self.mesh_count = self.mesh_count + 1
      self.mesh_names.append(mesh_name)
      m = m + 1
    mesh_list = self.fh.mesh_list()
    for name in mesh_list:
      self.assertTrue(name in self.mesh_names)

class TestOutputSimulations(unittest.TestCase):

  def setUp(self):

    import random

    self.filename = 'test-Output.h5'
    self.fh       = Output(self.filename,'w')

  def tearDown(self):
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def test_sim_list(self):
    import random
    list=self.fh.simulation_list()
    self.assertEqual(len(list),0)
    n=random.randint(10,20)
    i=0
    save_names=[]
    while i < n:
      sim_name='Simulation ' + str(i)
      save_names.insert(i,sim_name)
      self.fh.add_simulation(sim_name)
      i=i+1
    list=self.fh.simulation_list()
    self.assertEqual(len(list),n)
    for name in save_names:
      if not name in list:
	raise RuntimeError, 'Name not found in list'

  def test_add_simulation(self):
    sim_name='New Simulation'
    self.fh.add_simulation(sim_name)

  def test_sim_count(self):
    import random
    sim_count=self.fh.simulation_count()
    self.assertEqual(sim_count,0)
    n=random.randint(10,20)
    i=0
    while i < n:
      sim_name='Simulation ' + str(i)
      self.fh.add_simulation(sim_name)
      i=i+1
    cnt=self.fh.simulation_count()
    self.assertEqual(cnt,n)

  def test_sim_exists(self):
    sim_name='Sim DNE'
    self.assertEqual(self.fh.simulation_exists(sim_name),0)
    sim_name='Smin Does Exist'
    self.fh.add_simulation(sim_name)
    self.assertEqual(self.fh.simulation_exists(sim_name),1)

  def test_get_simulation(self):
    import Danu
    sim_name='Fetchme'
    try:
      sim=self.get_simulation(sim_name)
    except:
      print 'Caught get_simulation DNE test'
    else:
      raise RuntimeError, 'Failed to raise an error'
    self.fh.add_simulation(sim_name)
    sim=self.fh.get_simulation(sim_name)
    if not isinstance(sim,Danu.Simulation):
      raise RuntimeError, 'Failed to return simulation object'


if __name__ == '__main__':
   #unittest.main()
   suite=unittest.TestSuite()
   suite.addTest(TestOutputAccess('test_open_rdonly'))
   unittest.TextTestRunner().run(suite)



    


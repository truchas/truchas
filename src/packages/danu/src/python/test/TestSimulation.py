# ############################################################################ #
#
# Python UnitTest for Danu Simulation objects
#
# ############################################################################ #

import os, sys
import unittest
from unittest import TestCase

from Danu import Simulation

class TestUtils(TestCase):

  def setUp(self):
    import os
    from Danu import Output

    self.filename = 'test-Simulation.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Simulation')

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def test_data_count(self):
    cnt=self.sim.data_count()
    self.assertEqual(cnt,0)

  def test_data_exists(self):
    data_name='Dataset'
    self.assertEqual(self.sim.data_exists(data_name),0)

  def test_data_list(self):
    list=self.sim.data_list()
    self.assertEqual(len(list),0)

class TestAttributes(TestCase):

  def setUp(self):
    import os
    from Danu import Output

    self.filename = 'test-Simulation.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
    import random
    sim_name='Simulation For Testing'
    sim=Simulation(self.fh,sim_name)

    # integers
    int_name='Integer Attribute'
    try:
      sim.get_attribute(int_name)
    except:
      print 'Caught DNE attribute'
    else:
      raise RuntimeError, 'Failed to catch DNE error'
    int_attr=random.randint(1,102400)
    sim.set_attribute(int_name,int_attr)
    read_attr=sim.get_attribute(int_name)
    self.assertEqual(int_attr,read_attr)

    # double
    dbl_name='Double Attribute'
    dbl_attr=random.random()
    sim.set_attribute(dbl_name,dbl_attr)
    read_dbl=sim.get_attribute(dbl_name)
    self.assertEqual(read_dbl,dbl_attr)

    # string
    str_name='Chacracter Name Attribute'
    str_attr='ajksdweiouamnv 9iajemn  oiwiwe'
    sim.set_attribute(str_name,str_attr)
    string=sim.get_attribute(str_name)
    self.assertEqual(string,str_attr)


class TestWriteData(TestCase):

  def setUp(self):
    import os
    from Danu import Output

    self.filename = 'test-Simulation.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Simulation')

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()

  def test_write_int(self):
    import random
    import numpy
    data_name='Data Set 1D Integer'
    nsize=random.randint(1000,100000)
    data=numpy.zeros((nsize),numpy.int32)
    i=0
    while i < nsize:
      data[i]=random.randint(0,100000)
      i=i+1
    self.sim.data_write(data_name,data)  
    del data

    data_name='Data Set 2D Integer'
    n1=random.randint(100,1000)
    n2=random.randint(100,1000)
    data=numpy.zeros((n1,n2),numpy.int32)
    i=0
    j=0
    while i < n1:
      while j < n2:
        data[i][j]=random.randint(0,100000)
        j=j+1
      i=i+1
    self.sim.data_write(data_name,data)  
    del data

    data_name='Data Set 3D Integer'
    n1=random.randint(100,1000)
    n2=random.randint(100,1000)
    n3=random.randint(100,1000)
    data=numpy.zeros((n1,n2,n3),numpy.int32)
    i=0
    j=0
    k=0
    while i < n1:
      while j < n2:
	while k < n3:
	  data[i][j][k]=random.randint(0,100000)
	  k=k+1
        j=j+1
      i=i+1
    self.sim.data_write(data_name,data)  
    del data

  def test_write_double(self):
    import random
    import numpy
    data_name='Data Set 1D Double'
    nsize=random.randint(1000,100000)
    data=numpy.random.random_sample((nsize))
    self.sim.data_write(data_name,data)  
    del data

    data_name='Data Set 2D Double'
    n1=random.randint(100,1000)
    n2=random.randint(100,1000)
    data=numpy.random.random_sample((n1,n2))
    self.sim.data_write(data_name,data)  
    del data

    data_name='Data Set 3D Double'
    n1=random.randint(100,1000)
    n2=random.randint(100,1000)
    n3=random.randint(100,1000)
    data=numpy.random.random_sample((n1,n2,n3))
    self.sim.data_write(data_name,data)  
    del data


  def test_write_float(self):
    import random
    import numpy
    data_name='Data Set 1D Float'
    nsize=random.randint(1000,100000)
    data=numpy.float32(numpy.random.random_sample((nsize)))
    self.sim.data_write(data_name,data)  
    del data

    data_name='Data Set 2D Float'
    n1=random.randint(100,1000)
    n2=random.randint(100,1000)
    data=numpy.float32(numpy.random.random_sample((n1,n2)))
    self.sim.data_write(data_name,data)  
    del data

    data_name='Data Set 3D Float'
    n1=random.randint(10,100)
    n2=random.randint(10,100)
    n3=random.randint(10,100)
    data=numpy.float32(numpy.random.random_sample((n1,n2,n3)))
    self.sim.data_write(data_name,data)  
    del data


class TestDataReadInteger(TestCase):

  def setUp(self):
    import os
    import numpy
    import random
    from Danu import Output

    self.filename = 'test-Simulation.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Simulation')

    self.data_name_1d='Integer Data 1D'
    self.n=random.randint(1000,10000)
    self.data_1d=numpy.zeros((self.n),numpy.int32)
    i=0
    while i < self.n:
      self.data_1d[i]=random.randint(1,100000)
      i=i+1
    self.sim.data_write(self.data_name_1d,self.data_1d) 

    self.data_name_2d='Integer Data 2D'
    self.n1=random.randint(100,1000)
    self.n2=random.randint(100,1000)
    self.data_2d=numpy.zeros((self.n1,self.n2),numpy.int32)
    i=0
    while i < self.n1:
      j=0
      while j < self.n2:
	self.data_2d[i][j]=random.randint(1,100000)
	j=j+1
      i=i+1
    self.sim.data_write(self.data_name_2d,self.data_2d) 
    
    self.data_name_3d='Integer Data 3D'
    self.n1=random.randint(10,100)
    self.n2=random.randint(10,100)
    self.n3=random.randint(10,100)
    self.data_3d=numpy.zeros((self.n1,self.n2,self.n3),numpy.int32)
    i=0
    while i < self.n1:
      j=0
      while j < self.n2:
        k=0
	while k < self.n3:
	  self.data_3d[i][j][k]=random.randint(1,100000)
	  k=k+1
	j=j+1
      i=i+1
    self.sim.data_write(self.data_name_3d,self.data_3d) 

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()

  def runTest(self):
    try:
      data=self.sim.data_read('Data DNE')
    except:
      print 'Caught bad non-series data name'
    else:
      raise RuntimeError, 'Failed to catch DNE dataset name'

    data1=self.sim.data_read(self.data_name_1d)
    self.assertEqual(data1.all(),self.data_1d.all())

    data2=self.sim.data_read(self.data_name_2d)
    self.assertEqual(data2.all(),self.data_2d.all())

    data3=self.sim.data_read(self.data_name_3d)
    self.assertEqual(data3.all(),self.data_3d.all())


class TestDataReadFloat(TestCase):

  def setUp(self):
    import os
    import numpy
    import random
    from Danu import Output

    self.filename = 'test-Simulation.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Simulation')

    self.data_name_1d='Float Data 1D'
    self.n=random.randint(1000,10000)
    self.data_1d=numpy.float32(numpy.random.random_sample((self.n)))
    self.sim.data_write(self.data_name_1d,self.data_1d) 

    self.data_name_2d='Float Data 2D'
    self.n1=random.randint(100,1000)
    self.n2=random.randint(100,1000)
    self.data_2d=numpy.float32(numpy.random.random_sample((self.n1,self.n2)))
    self.sim.data_write(self.data_name_2d,self.data_2d) 
    
    self.data_name_3d='Float Data 3D'
    self.n1=random.randint(10,100)
    self.n2=random.randint(10,100)
    self.n3=random.randint(10,100)
    self.data_3d=numpy.float32(numpy.random.random_sample((self.n1,self.n2,self.n3)))
    self.sim.data_write(self.data_name_3d,self.data_3d) 

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()

  def runTest(self):
    data1=self.sim.data_read(self.data_name_1d)
    self.assertEqual(data1.all(),self.data_1d.all())

    data2=self.sim.data_read(self.data_name_2d)
    self.assertEqual(data2.all(),self.data_2d.all())

    data3=self.sim.data_read(self.data_name_3d)
    self.assertEqual(data3.all(),self.data_3d.all())

class TestDataReadDouble(TestCase):

  def setUp(self):
    import os
    import numpy
    import random
    from Danu import Output

    self.filename = 'test-Simulation.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Simulation')

    self.data_name_1d='Double Data 1D'
    self.n=random.randint(1000,10000)
    self.data_1d=numpy.random.random_sample((self.n))
    self.sim.data_write(self.data_name_1d,self.data_1d) 

    self.data_name_2d='Double Data 2D'
    self.n1=random.randint(100,1000)
    self.n2=random.randint(100,1000)
    self.data_2d=numpy.random.random_sample((self.n1,self.n2))
    self.sim.data_write(self.data_name_2d,self.data_2d) 
    
    self.data_name_3d='Double Data 3D'
    self.n1=random.randint(10,100)
    self.n2=random.randint(10,100)
    self.n3=random.randint(10,100)
    self.data_3d=numpy.random.random_sample((self.n1,self.n2,self.n3))
    self.sim.data_write(self.data_name_3d,self.data_3d) 

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()

  def runTest(self):
    data1=self.sim.data_read(self.data_name_1d)
    self.assertEqual(data1.all(),self.data_1d.all())

    data2=self.sim.data_read(self.data_name_2d)
    self.assertEqual(data2.all(),self.data_2d.all())

    data3=self.sim.data_read(self.data_name_3d)
    self.assertEqual(data3.all(),self.data_3d.all())

class TestSequence(TestCase):
  
  def setUp(self):
    import os
    import numpy
    import random
    from Danu import Output

    self.filename = 'test-Simulation.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Simulation')

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()

  def test_basic(self):
    from Danu import Sequence
    import numpy
    import random
    cycle=random.randint(10,100)
    time=numpy.random.random_sample((1))[0]
    seq=self.sim.get_nextSequence(cycle,time);
    print seq.id
    print seq.time
    print seq.cycle

  def test_sequence_count(self):
    import random
    import numpy
    cnt=self.sim.sequence_count()
    self.assertEqual(cnt,0,'Failed to return zero count')
    gen_cnt=random.randint(1,100)
    cycle=1
    i=1
    time=numpy.random.random_sample((1))[0]
    while i <= gen_cnt:
      time=time+numpy.random.random_sample((1))[0]
      seq=self.sim.get_nextSequence(cycle,time)
      del seq
      i=i+1

    cnt=self.sim.sequence_count()
    self.assertEqual(cnt,gen_cnt)

  def test_sequence_exists(self):
    import random
    import numpy
    cycle=1
    time=numpy.random.random_sample((1))[0]

    flag=self.sim.sequence_exists('Series DNE')
    self.assertEqual(flag,0,'Failed to return false status')
    seq=self.sim.get_nextSequence(cycle,time)
    id=seq.id
    del seq
    seq_name=self.sim.get_sequence_name(id)
    flag=self.sim.sequence_exists(seq_name)
    self.assertEqual(flag,1,'Failed to return true status')

  def test_sequence_list(self):
    import random
    import numpy
    list=self.sim.sequence_list()
    self.assertEqual(len(list),0,'Failed to return empty list')
    gen_cnt=random.randint(1,100)
    cycle=1
    i=1
    time=numpy.random.random_sample((1))[0]
    names=[]
    while i <= gen_cnt:
      time=time+numpy.random.random_sample((1))[0]
      seq=self.sim.get_nextSequence(cycle,time)
      seq_name=self.sim.get_sequence_name(seq.id)
      names.append(seq_name)
      del seq
      i=i+1

    list=self.sim.sequence_list()
    self.assertEqual(len(list),gen_cnt,'Failed to return correct list size')
    for name in names:
      self.assertTrue(name in list)

      



if __name__ == '__main__':
   unittest.main()
   #suite=unittest.TestSuite()
   #tests=[TestSequence.test_sequence_count]
   #suite.addTest(TestSequence('test_sequence_list'))
   #unittest.TextTestRunner().run(suite)





    


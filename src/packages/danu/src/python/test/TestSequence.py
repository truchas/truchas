# ############################################################################ #
#
# Python UnitTest for Danu Sequence objects
#
# ############################################################################ #

import os, sys
import unittest
from unittest import TestCase

from Danu import Sequence

class TestBasic(TestCase):

  def setUp(self):
    import os
    from Danu import Output

    self.filename = 'test-Sequence.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Sequence')

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
    import random
    import numpy

    print 'Begin basic test'
    try:
      seq=Sequence()
    except:
      print 'Caught invalid number of arguments'

    try:
      seq=Sequence(self.sim)
    except:
      print 'Caught invalid name or id'

    try:
      seq=Sequence(self.sim,'Series DNE')
    except:
      print 'Caught id group that does not exist'

    try:
      seq=Sequence(self.sim,seriesname='Series 10')
    except:
      print 'Caught group name that does not exist'

    cycle=random.randint(0,10)
    time=numpy.float64(random.randint(0,1000))
    seq=self.sim.get_nextSequence(cycle,time)
    id=seq.id
    del seq
    seq_name=self.sim.get_sequence_name(id)
    new_seq=Sequence(self.sim,seq_name)
    del new_seq


class TestWriteData(TestCase):

  def setUp(self):
    import os
    from Danu import Output

    self.filename = 'test-Sequence.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Sequence')

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def test_write_1D(self):
    import numpy
    import random

    seq=self.sim.get_nextSequence(0,0.0)

    size=random.randint(128,512)
    idata=numpy.zeros((size))
    i=0
    while i < size:
      idata[i]=random.randint(0,100000)
      i=i+1
    seq.data_write('Integer Data',idata)

    size=random.randint(128,512)
    data=numpy.float32(numpy.random.random_sample((size)))
    seq.data_write('Float Data',data)

    size=random.randint(128,512)
    data=numpy.random.random_sample((size))
    seq.data_write('Double Data',data)

  def test_write_2D(self):
    import numpy
    import random

    seq=self.sim.get_nextSequence(0,0.0)

    n1=random.randint(10,128)
    n2=random.randint(10,128)
    idata=numpy.zeros((n1,n2))
    i=0
    while i < n1:
      j=0
      while j < n2:
        idata[i][j]=random.randint(1,1000)
        j=j+1
      i=i+1
    seq.data_write('Integer Data',idata)

    data=numpy.float32(numpy.random.random_sample((n1,n2)))
    seq.data_write('Float Data',data)

    data=numpy.random.random_sample((n1,n2))
    seq.data_write('Double Data',data)

  def test_write_3D(self):
    import numpy
    import random

    seq=self.sim.get_nextSequence(0,0.0)

    n1=random.randint(1,20)
    n2=random.randint(1,20)
    n3=random.randint(1,20)
    idata=numpy.zeros((n1,n2,n3))
    i=0
    while i < n1:
      j=0
      while j < n2:
        k=0
        while k < n3:
          idata[i][j][k]=random.randint(1,1000)
          k=k+1
        j=j+1
      i=i+1
    seq.data_write('Integer Data',idata)

    data=numpy.float32(numpy.random.random_sample((n1,n2,n3)))
    seq.data_write('Float Data',data)

    data=numpy.random.random_sample((n1,n2,n3))
    seq.data_write('Double Data',data)

class TestReadInteger1D(TestCase):

  def setUp(self):
    import os
    from Danu import Output
    import random
    import numpy

    self.filename = 'test-Sequence.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Sequence')

    seq=self.sim.get_nextSequence(0,0.0)
    n=random.randint(128,512)
    self.data=numpy.zeros((n))
    i=0
    while i < n:
      self.data[i]=random.randint(0,10000)
      i=i+1

    self.data_name='Integer Data'
    seq.data_write(self.data_name,self.data)

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
  
    seq=self.sim.get_sequence('Series 1')
    try:
      read_data=seq.data_read('Data DNE')
    except:
      print 'Caught Data set that did not exist'

    read_data=seq.data_read(self.data_name)

    self.assertEqual(read_data.all(), self.data.all())
                     

class TestReadInteger2D(TestCase):

  def setUp(self):
    import os
    from Danu import Output
    import random
    import numpy

    self.filename = 'test-Sequence.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Sequence')

    seq=self.sim.get_nextSequence(0,0.0)
    n1=random.randint(128,512)
    n2=random.randint(128,512)
    self.data=numpy.zeros((n1,n2))
    i=0
    while i < n1:
      j=0
      while j < n2:
        self.data[i][j]=random.randint(0,10000)
        j=j+1
      i=i+1

    self.data_name='Integer Data'
    seq.data_write(self.data_name,self.data)

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
  
    seq=self.sim.get_sequence('Series 1')
    read_data=seq.data_read(self.data_name)

    self.assertEqual(read_data.all(), self.data.all())
                     

class TestReadInteger3D(TestCase):

  def setUp(self):
    import os
    from Danu import Output
    import random
    import numpy

    self.filename = 'test-Sequence.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Sequence')

    seq=self.sim.get_nextSequence(0,0.0)
    n1=random.randint(10,20)
    n2=random.randint(10,20)
    n3=random.randint(10,20)
    self.data=numpy.zeros((n1,n2,n3))
    i=0
    while i < n1:
      j=0
      while j < n2:
        k=0
        while k < n3:
          self.data[i][j]=random.randint(0,10000)
          k=k+1
        j=j+1
      i=i+1

    self.data_name='Integer Data'
    seq.data_write(self.data_name,self.data)

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
  
    seq=self.sim.get_sequence('Series 1')
    read_data=seq.data_read(self.data_name)

    self.assertEqual(read_data.all(), self.data.all())
     
class TestReadFloat1D(TestCase):

  def setUp(self):
    import os
    from Danu import Output
    import random
    import numpy

    self.filename = 'test-Sequence.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Sequence')

    seq=self.sim.get_nextSequence(0,0.0)
    n=random.randint(128,512)
    self.data=numpy.float32(numpy.random.random_sample((n)))

    self.data_name='Float Data'
    seq.data_write(self.data_name,self.data)

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
  
    seq=self.sim.get_sequence('Series 1')
    read_data=seq.data_read(self.data_name)
    self.assertEqual(read_data.all(), self.data.all())

class TestReadFloat2D(TestCase):

  def setUp(self):
    import os
    from Danu import Output
    import random
    import numpy

    self.filename = 'test-Sequence.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Sequence')

    seq=self.sim.get_nextSequence(0,0.0)
    n1=random.randint(128,512)
    n2=random.randint(128,512)
    self.data=numpy.float32(numpy.random.random_sample((n1,n2)))

    self.data_name='Float Data'
    seq.data_write(self.data_name,self.data)

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
  
    seq=self.sim.get_sequence('Series 1')
    read_data=seq.data_read(self.data_name)
    self.assertEqual(read_data.all(), self.data.all())


class TestReadFloat3D(TestCase):

  def setUp(self):
    import os
    from Danu import Output
    import random
    import numpy

    self.filename = 'test-Sequence.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Sequence')

    seq=self.sim.get_nextSequence(0,0.0)
    n1=random.randint(12,50)
    n2=random.randint(10,60)
    n3=random.randint(2,40)
    self.data=numpy.float32(numpy.random.random_sample((n1,n2,n3)))

    self.data_name='Float Data'
    seq.data_write(self.data_name,self.data)

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
  
    seq=self.sim.get_sequence('Series 1')
    read_data=seq.data_read(self.data_name)
    self.assertEqual(read_data.all(), self.data.all())

class TestReadDouble1D(TestCase):

  def setUp(self):
    import os
    from Danu import Output
    import random
    import numpy

    self.filename = 'test-Sequence.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Sequence')

    seq=self.sim.get_nextSequence(0,0.0)
    n=random.randint(128,512)
    self.data=numpy.random.random_sample((n))

    self.data_name='Double Data'
    seq.data_write(self.data_name,self.data)

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
  
    seq=self.sim.get_sequence('Series 1')
    read_data=seq.data_read(self.data_name)
    self.assertEqual(read_data.all(), self.data.all())

class TestReadDouble2D(TestCase):

  def setUp(self):
    import os
    from Danu import Output
    import random
    import numpy

    self.filename = 'test-Sequence.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Sequence')

    seq=self.sim.get_nextSequence(0,0.0)
    n1=random.randint(128,512)
    n2=random.randint(128,512)
    self.data=numpy.random.random_sample((n1,n2))

    self.data_name='Double Data'
    seq.data_write(self.data_name,self.data)

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
  
    seq=self.sim.get_sequence('Series 1')
    read_data=seq.data_read(self.data_name)
    self.assertEqual(read_data.all(), self.data.all())


class TestReadDouble3D(TestCase):

  def setUp(self):
    import os
    from Danu import Output
    import random
    import numpy

    self.filename = 'test-Sequence.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Sequence')

    seq=self.sim.get_nextSequence(0,0.0)
    n1=random.randint(12,50)
    n2=random.randint(10,60)
    n3=random.randint(2,40)
    self.data=numpy.random.random_sample((n1,n2,n3))

    self.data_name='Double Data'
    seq.data_write(self.data_name,self.data)

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
  
    seq=self.sim.get_sequence('Series 1')
    read_data=seq.data_read(self.data_name)
    self.assertEqual(read_data.all(), self.data.all())


class TestAttributes(TestCase):
 
  def setUp(self):
    import os
    from Danu import Output
    import random
    import numpy

    self.filename = 'test-Sequence.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Sequence')
    self.seq=self.sim.get_nextSequence(0,0.0)

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
    import random
    
    seq=self.seq

    # integers
    int_name='Integer Attribute'
    try:
      seq.get_attribute(int_name)
    except:
      print 'Caught DNE attribute'
    else:
      raise RuntimeError, 'Failed to catch DNE error'
    int_attr=random.randint(1,102400)
    seq.set_attribute(int_name,int_attr)
    read_attr=seq.get_attribute(int_name)
    self.assertEqual(int_attr,read_attr)

    # double
    dbl_name='Double Attribute'
    dbl_attr=random.random()
    seq.set_attribute(dbl_name,dbl_attr)
    read_dbl=seq.get_attribute(dbl_name)
    self.assertEqual(read_dbl,dbl_attr)

    # string
    str_name='Chacracter Name Attribute'
    str_attr='ajksdweiouamnv 9iajemn  oiwiwe'
    seq.set_attribute(str_name,str_attr)
    string=seq.get_attribute(str_name)
    self.assertEqual(string,str_attr)


if __name__ == '__main__':
   unittest.main()
   #suite=unittest.TestSuite()
   #suite.addTest(TestDataReadDouble('runTest'))
   #unittest.TextTestRunner().run(suite)





    


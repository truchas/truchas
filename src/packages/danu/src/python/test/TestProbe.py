# ############################################################################ #
#
# Python UnitTest for Danu Probe objects
#
# ############################################################################ #

import os, sys
import unittest
from unittest import TestCase

from Danu import Probe

class TestBasic(TestCase):

  def setUp(self):
    import os
    from Danu import Output

    self.filename = 'test-Probe.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Probe')

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
      probe=Probe()
    except:
      print 'Caught invalid number of arguments'

    try:
      probe=Probe(self.sim)
    except:
      print 'Caught invalid number of arguments'

    try:
      probe=Probe(self.sim,'Dummy Probe')
    except:
      print 'Caught invalid number of arguments'

    try:
      data=numpy.ones((10))
      probe=Probe(self.sim,'Dummy Probe',data)
    except:
      print 'Caught incorrect argument type'

    # Simple constructor and delete
    len=2
    num=20
    data=numpy.ones((num,len))
    probe= Probe(self.sim,'Dummy Probe',data)
    del probe

    probe=self.sim.probe_open('Dummy Probe')
    self.assertEqual(probe.length(),len)
    self.assertEqual(probe.number(),num)




class TestAttributes(TestCase):
 
  def setUp(self):
    import os
    from Danu import Output
    import random
    import numpy

    self.filename = 'test-Probe.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Probe')
    self.probe_name='Dummy Probe'
    data=numpy.ones((3,10))
    Probe(self.sim,self.probe_name,data)

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
    import random
    
    probe=self.sim.probe_open(self.probe_name)

    # integers
    int_name='Integer Attribute'
    try:
      probe.get_attribute(int_name)
    except:
      print 'Caught DNE attribute'
    else:
      raise RuntimeError, 'Failed to catch DNE error'
    int_attr=random.randint(1,102400)
    probe.set_attribute(int_name,int_attr)
    read_attr=probe.get_attribute(int_name)
    self.assertEqual(int_attr,read_attr)

    # double
    dbl_name='Double Attribute'
    dbl_attr=random.random()
    probe.set_attribute(dbl_name,dbl_attr)
    read_dbl=probe.get_attribute(dbl_name)
    self.assertEqual(read_dbl,dbl_attr)

    # string
    str_name='Chacracter Name Attribute'
    str_attr='ajksdweiouamnv 9iajemn  oiwiwe'
    probe.set_attribute(str_name,str_attr)
    string=probe.get_attribute(str_name)
    self.assertEqual(string,str_attr)


class TestCreate(TestCase):

  def setUp(self):
    import os
    from Danu import Output

    self.filename = 'test-Probe.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Probe')

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def test_create_int(self):
    import numpy
    import random

    len=random.randint(1,7)
    num=random.randint(128,512)
    idata=numpy.zeros((num,len))
    i=0
    while i < num:
      l=0
      while l < len:
	idata[i][l]=random.randint(0,100000)
	l=l+1
      i=i+1

    probe=Probe(self.sim,'Int Data',idata)  

  def test_create_float(self):
    import numpy
    import random

    len=random.randint(1,7)
    num=random.randint(128,512)
    data=numpy.float32(numpy.random.random_sample((num,len)))

    probe=Probe(self.sim,'Float data',data)  


  def test_create_double(self):
    import numpy
    import random

    len=random.randint(1,7)
    num=random.randint(128,512)
    data=numpy.random.random_sample((num,len))

    probe=Probe(self.sim,'Double data',data)  
 

class TestWrite(TestCase):

  def setUp(self):
    import os
    import random
    import numpy
    from Danu import Output

    self.filename = 'test-Probe.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Probe')

    # Int data
    self.idata_name='Int data'
    len=random.randint(1,7)
    num=random.randint(128,512)
    self.idata=numpy.zeros((num,len))
    self.idata_len=len
    i=0
    while i < num:
      l=0
      while l < len:
	self.idata[i][l]=random.randint(0,100000)
	l=l+1
      i=i+1
    probe=Probe(self.sim,self.idata_name,self.idata)  

    # Float data
    self.fdata_name='Float data'
    len=random.randint(1,7)
    num=random.randint(128,512)
    self.fdata=numpy.float32(numpy.random.random_sample((num,len)))
    self.fdata_len=len
    probe=Probe(self.sim,self.fdata_name,self.fdata)  

    # Double data
    self.rdata_name='Double data'
    len=random.randint(1,7)
    num=random.randint(128,512)
    self.rdata=numpy.random.random_sample((num,len))
    self.rdata_len=len
    probe=Probe(self.sim,self.rdata_name,self.rdata)  

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def test_write_int(self):
    import numpy
    import random

    len=self.idata_len + 1
    num=random.randint(128,512)
    idata=numpy.zeros((num,len))
    probe=self.sim.probe_open(self.idata_name)
    try:
      probe.write(idata)
    except:
      print 'Caught error when trying to write mis-sized data'

    idata=numpy.zeros((num,self.idata_len))
    probe.write(idata)


  def test_write_float(self):
    import numpy
    import random

    len=self.fdata_len+1
    num=random.randint(128,512)
    data=numpy.float32(numpy.random.random_sample((num,len)))

    probe=self.sim.probe_open(self.fdata_name)

    try:
      probe.write(data)
    except:
      print 'Caught error when trying to write mis-sized data'

    data=numpy.float32(numpy.random.random_sample((num,self.fdata_len)))
    probe.write(data)


  def test_write_double(self):
    import numpy
    import random

    len=self.rdata_len+1
    num=random.randint(128,512)
    data=numpy.random.random_sample((num,len))
    probe=self.sim.probe_open(self.rdata_name)

    try:
      probe.write(data)
    except:
      print 'Caught mis-sized data array'

    len=self.rdata_len  
    data=numpy.random.random_sample((num,len))
    probe.write(data)
 


class TestReadData(TestCase):

  def setUp(self):
    import os
    from Danu import Output
    import random
    import numpy

    self.filename='test-Probe.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')
    self.sim=self.fh.add_simulation('Test Probe')

    # Int data
    self.idata_name='Int data'
    len=random.randint(1,7)
    num=random.randint(128,512)
    self.idata=numpy.zeros((num,len))
    i=0
    while i < num:
      l=0
      while l < len:
	self.idata[i][l]=random.randint(0,100000)
	l=l+1
      i=i+1
    probe=Probe(self.sim,self.idata_name,self.idata)  

    # Float data
    self.fdata_name='Float data'
    len=random.randint(1,7)
    num=random.randint(128,512)
    self.fdata=numpy.float32(numpy.random.random_sample((num,len)))
    probe=Probe(self.sim,self.fdata_name,self.fdata)  

    # Double data
    self.rdata_name='Double data'
    len=random.randint(1,7)
    num=random.randint(128,512)
    self.rdata=numpy.random.random_sample((num,len))
    probe=Probe(self.sim,self.rdata_name,self.rdata)  

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def test_read_int(self):
 
    probe=self.sim.probe_open(self.idata_name)
    read_data=probe.read()
    self.assertEqual(read_data.all(), self.idata.all())

  def test_read_float(self):
 
    probe=self.sim.probe_open(self.fdata_name)
    read_data=probe.read()
    self.assertEqual(read_data.all(), self.fdata.all())

  def test_read_double(self):

    probe=self.sim.probe_open(self.rdata_name)
    read_data=probe.read()
    self.assertEqual(read_data.all(), self.rdata.all())

if __name__ == '__main__':
   unittest.main()
   #suite=unittest.TestSuite()
   #suite.addTest(TestReadData('test_read_int'))
   #result=unittest.TextTestRunner().run(suite)
   #if result.wasSuccessful():
   #  print 'Test suite passed'
   #  sys.exit(0)
   #else:
   #  print result.printErrors()
   #  sys.exit(1)








    


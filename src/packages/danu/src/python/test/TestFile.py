# ############################################################################ #
#
# Python UnitTest for Files
#
# Really this is unit test of the H5Object, Can not create and  H5Object
# without generating an HDF5 file.
# 
# ############################################################################ #

import os, sys
import unittest

from Danu import File
  



class TestFileAccess(unittest.TestCase):

  def setUp(self):
    self.filename         = 'test-File.h5'

  def tearDown(self):
    if os.path.exists(self.filename):
      os.remove(self.filename)

  def test_open_dflt(self):
    file = File(self.filename)
    del file

  def test_open_rdonly(self):
    file = File(self.filename,'r')
    del file

  def test_open_rdwr(self):
    file = File(self.filename,'w')
    del file

  def test_open_append(self):
    file = File(self.filename,'a')
    del file

  def test_flush(self):
    file = File(self.filename,'a')
    file.flush()
    file.flush(1)
    del file

class TestFileWriteAttribute(unittest.TestCase):

  def setUp(self):
    self.filename         = 'test-File.h5'
    self.int_attr_name    = 'Dummy Int Attribute'
    self.double_attr_name = 'Dummy Double Attribute'
    self.str_attr_name    = 'Dummy String Attribute'
    self.fh               = File(self.filename,'a')

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

class TestFileReadAttribute(unittest.TestCase):

  def setUp(self):
    import random
    import string

    self.filename         = 'test-File.h5'

    self.fh               = File(self.filename,'w')

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
      if attr not in attributes: fail('Failed to read attributes correctly')

if __name__ == '__main__':
   unittest.main()



    


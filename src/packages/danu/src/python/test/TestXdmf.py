# ############################################################################ #
#
# Python UnitTest for Xdmf.py
#
# Really this is unit test suite for the Xdmf module
# 
# ############################################################################ #

import os, sys
import unittest

import xml.etree.ElementTree as ET

import Xdmf

TEST_LOG='stuff.xml'

def print_out(a):
  if os.path.exists(TEST_LOG):
    fh = open(TEST_LOG,'a')
  else:
    fh = open(TEST_LOG,'w')
  fh.write("BEGIN TEST\n")
  Xdmf.indent(a.elem)
  ET.ElementTree(a.elem).write(fh)
  fh.write("END TEST\n")
  fh.write("\n\n")
  fh.close()

class TestXdmf(unittest.TestCase):

  def setUp(self):
    #if os.path.exists(TEST_LOG):
    #  os.remove(TEST_LOG)
    self.root = Xdmf.XdmfRoot()

  def tearDown(self):
    print_out(self.root)

  def test_root(self):
    root = Xdmf.XdmfRoot()
    print_out(root)

  def test_time(self):
    time = Xdmf.XdmfTime(1.0)
    t1=1.0
    t2=2.0
    self.assertEqual(t1,time.get_value(),'Incorrect time value returned')
    time.set_value(t2)
    self.assertEqual(t2,time.get_value(),'Incorrect time value')
    self.root.append(time)

  def test_info(self):
    self.failUnlessRaises(TypeError,Xdmf.XdmfInformation)
    info=Xdmf.XdmfInformation('Dummy')
    msg='This is a test'
    info.set_value(msg)
    cmp = info.get_value()
    self.assertEqual(msg,cmp,'Wrong value information element')
    self.root.append(info)

  def test_geometry(self):
    self.failUnlessRaises(TypeError,Xdmf.XdmfGeometry)
    geom=Xdmf.XdmfGeometry('geo','X_Y_Z')
    self.root.append(geom)

  def test_domain(self):
    self.failUnlessRaises(TypeError,Xdmf.XdmfDomain)
    dom=Xdmf.XdmfDomain('Test domain')
    self.root.append(dom)

  def test_grid(self):
    self.failUnlessRaises(TypeError,Xdmf.XdmfGrid)
    g=Xdmf.XdmfGrid('somegrid')
    t1=1.0
    t2=2.0
    g.set_time(t1)
    self.assertEqual(t1,g.get_time(),'Incorrect time returned')
    self.root.append(g)


  def test_topology(self):
    self.failUnlessRaises(TypeError,Xdmf.XdmfTopology)
    self.failUnlessRaises(TypeError,Xdmf.XdmfTopology,'My Mesh Topo')
    top = Xdmf.XdmfTopology('My Other Mesh','Hexahedron')
    offset=3
    top.set_offset(offset)
    self.assertEqual(offset,top.get_offset(),'Wrong offset')
    dim=2
    top.set_dimensions(dim)
    self.assertEqual(dim,top.get_dimensions(),'Wrong dimensions')
    order=2
    top.set_order(order)
    self.assertEqual(order,top.get_order(),'Wrong order')
    num_of_nodes=8
    top.set_npelem(num_of_nodes)
    self.assertEqual(num_of_nodes,top.get_npelem(),'Wrong number of nodes/element')
    num_nodes=10000
    top.set_nelem(num_nodes)
    self.assertEqual(num_nodes,top.get_nelem(),'Wrong number of elements')
    self.root.append(top)

  def test_attribute(self):
    self.failUnlessRaises(TypeError,Xdmf.XdmfAttribute)
    a = Xdmf.XdmfAttribute('Test Attribute')
    type='Vector'
    a.set_type(type)
    self.assertEqual(type,a.get_type(),'Incorrect type returned')
    center='Cell'
    a.set_center(center)
    self.assertEqual(center,a.get_center(),'Incorrect center returned')
    self.root.append(a)

  def test_data_item(self):
    self.failUnlessRaises(TypeError,Xdmf.XdmfDataItem)
    d = Xdmf.XdmfDataItem('Test Data Item')
    type='tree'
    d.set_type(type)
    self.assertEqual(type,d.get_type(),'Incorrect type returned')
    dims=[2, 4]
    d.set_dimensions(dims)
    self.assertEqual(dims,d.get_dimensions(),'Incorrect dimensions (1)')
    dims=100
    d.set_dimensions(dims)
    self.assertEqual(dims,d.get_dimensions(),'Incorrect dimensions (2)')
    num_type='Int'
    d.set_num_type(num_type)
    self.assertEqual(num_type,d.get_num_type(),'Incorrect num type returned')
    p=8
    d.set_precision(p)
    self.assertEqual(p,d.get_precision(),'Incorrect precision returned')
    format='HDF'
    d.set_format(format)
    self.assertEqual(format,d.get_format(),'Incorrect precision returned')
    self.root.append(d)

class TestXdmfModel(unittest.TestCase):

  def setUp(self):
    self.model = Xdmf.XdmfModel()
    self.fh = open('deleteme.xml','w')

  def test_write(self):
    self.model.write(self.fh)

  def test_domain(self):
    self.model.add_domain('My Simulation')




if __name__ == '__main__':
   unittest.main()


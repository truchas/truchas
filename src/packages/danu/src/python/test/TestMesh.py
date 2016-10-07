# ############################################################################ #
#
# Python UnitTest for Danu Mesh objects
#
# ############################################################################ #

import os, sys
import unittest
from unittest import TestCase

from Danu import Mesh

class TestMeshUtils(unittest.TestCase):
  
  def setUp(self):
    from Danu import Output
    from Danu import UNSTRUCTURED_MESH, STRUCTURED_MESH
    from Danu import LINE_ELEM, TRI_ELEM, QUAD_ELEM, TET_ELEM, HEX_ELEM 

    self.filename = 'test-Mesh.h5'
    if os.path.exists(self.filename):
      os.remove(self.filename)
    self.fh       = Output(self.filename,'w')

    self.mesh_count = 0
    self.mesh_names = []

    self.valid_mesh_types = [UNSTRUCTURED_MESH, STRUCTURED_MESH]
    self.valid_mesh_elems1 = [LINE_ELEM]
    self.valid_mesh_elems2 = [TRI_ELEM, QUAD_ELEM]
    self.valid_mesh_elems3 = [TET_ELEM, HEX_ELEM]

  def tearDown(self):
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def test_mesh_open(self):
    import Danu 
    mesh_name='Test Mesh DNE'
    self.failUnlessRaises(RuntimeError,self.fh.get_mesh,mesh_name) 
    Mesh(self.fh,mesh_name,Danu.UNSTRUCTURED_MESH,Danu.LINE_ELEM)
    mesh=self.fh.get_mesh(mesh_name)

  def test_mesh_destory(self):
    import Danu
    mesh_name='Test Mesh'
    mesh=Mesh(self.fh,mesh_name,Danu.UNSTRUCTURED_MESH,Danu.LINE_ELEM)
    del mesh

  def test_mesh_create_1D(self):
    import random
    import string
    for mesh_type in self.valid_mesh_types:
      for elem in self.valid_mesh_elems1:
        mesh_name = ''.join(random.choice(string.letters) for i in xrange(16))
        mesh = Mesh(self.fh,mesh_name,mesh_type,elem)
        self.assertEqual(mesh.dim,1, 'Incorrect dimension')
        self.assertEqual(mesh.mesh_type,mesh_type,'Incorrect mesh type')
        self.assertEqual(mesh.elem_type,elem,'Incorrect mesh element type')

  def test_mesh_create_2D(self):
    import random
    import string
    for mesh_type in self.valid_mesh_types:
      for elem in self.valid_mesh_elems2:
        mesh_name = ''.join(random.choice(string.letters) for i in xrange(16))
        mesh = Mesh(self.fh,mesh_name,mesh_type,elem)
        self.assertEqual(mesh.dim,2, 'Incorrect dimension')
        self.assertEqual(mesh.mesh_type,mesh_type,'Incorrect mesh type')
        self.assertEqual(mesh.elem_type,elem,'Incorrect mesh element type')

  def test_mesh_create_3D(self):
    import random
    import string
    for mesh_type in self.valid_mesh_types:
      for elem in self.valid_mesh_elems3:
        mesh_name = ''.join(random.choice(string.letters) for i in xrange(16))
        mesh = Mesh(self.fh,mesh_name,mesh_type,elem)
        self.assertEqual(mesh.dim,3, 'Incorrect dimension')
        self.assertEqual(mesh.mesh_type,mesh_type,'Incorrect mesh type')
        self.assertEqual(mesh.elem_type,elem,'Incorrect mesh element type')


class TestMeshWriteCoordinates(unittest.TestCase):

  def setUp(self):
    import os
    from Danu import Output
    from Danu import UNSTRUCTURED_MESH, STRUCTURED_MESH
    from Danu import LINE_ELEM, TRI_ELEM, QUAD_ELEM, TET_ELEM, HEX_ELEM 

    self.filename = 'test-Mesh.h5'
    self.fh       = Output(self.filename,'w')

    self.mesh_count = 0
    self.mesh_names = []

    self.valid_mesh_types = [UNSTRUCTURED_MESH, STRUCTURED_MESH]
    self.valid_mesh_elems1 = [LINE_ELEM]
    self.valid_mesh_elems2 = [TRI_ELEM, QUAD_ELEM]
    self.valid_mesh_elems3 = [TET_ELEM, HEX_ELEM]

  def tearDown(self):
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def test_write_coordinates_1D(self):
    import numpy
    import random
    from Danu import LINE_ELEM

    mesh_name = '1D Mesh'
    nnodes=random.randint(128,1024)
    mesh = self.fh.add_unstruct_mesh(mesh_name,LINE_ELEM)
    xcoordinates = numpy.random.random_sample((nnodes))
    try:
      mesh.write_coordinates()
      raise RuntimeError, "Failed to raise exception"
    except:
      print "Caught null argument exception"

    try:
      mesh.write_coordinates(1)
      raise RuntimeError, "Failed to raise exception"
    except:
      print "Caught invalid data type"

    mesh.write_coordinates(xcoordinates)


class TestMeshReadCoordinates(unittest.TestCase):

  def _writeUnstructMesh(self,name,ndim,elem):
    import numpy
    import random
    nnodes=random.randint(10,512)
    mesh = self.fh.add_unstruct_mesh(name,elem)
    x=numpy.random.random_sample((nnodes))
    if ndim == 1:
      mesh.write_coordinates(x) 
      return x
    if ndim == 2:
      y=numpy.random.random_sample((nnodes))
      mesh.write_coordinates(x,y)
      return (x,y)
    else:
      y=numpy.random.random_sample((nnodes))
      z=numpy.random.random_sample((nnodes))
      mesh.write_coordinates(x,y,z)
      return (x,y,z)

  def setUp(self):
    import os
    from Danu import Output
    from Danu import UNSTRUCTURED_MESH, STRUCTURED_MESH
    from Danu import LINE_ELEM, TRI_ELEM, QUAD_ELEM, TET_ELEM, HEX_ELEM 

    self.filename = 'test-Mesh.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh       = Output(self.filename,'w')

    self.valid_mesh_types = [UNSTRUCTURED_MESH, STRUCTURED_MESH]
    self.valid_mesh_elems1 = [LINE_ELEM]
    self.valid_mesh_elems2 = [TRI_ELEM, QUAD_ELEM]
    self.valid_mesh_elems3 = [TET_ELEM, HEX_ELEM]

    self.valid_mesh_1D = []
    self.valid_mesh_2D = []
    self.valid_mesh_3D = []
     
    for elem in self.valid_mesh_elems1:
      mesh_name='Mesh 1D' + str(elem)
      self._writeUnstructMesh(mesh_name,1,elem)
      self.valid_mesh_1D.insert(0,mesh_name)

    for elem in self.valid_mesh_elems2:
      mesh_name='Mesh 2D' + str(elem)
      self._writeUnstructMesh(mesh_name,2,elem)
      self.valid_mesh_2D.insert(0,mesh_name)

    for elem in self.valid_mesh_elems3:
      mesh_name='Mesh 3D' + str(elem)
      self._writeUnstructMesh(mesh_name,3,elem)
      self.valid_mesh_3D.insert(0,mesh_name)


  def tearDown(self):
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def test_read_coord_1D(self):
    import numpy
    for mesh_name in self.valid_mesh_1D:
      mesh=self.fh.get_mesh(mesh_name)
      x=mesh.read_coordinates(0,1024)
    mesh=self.fh.get_mesh(self.valid_mesh_1D[0])
    try:
      mesh.read_coordiantes(1,1024)
    except:
      print 'Caught invalid dimension request'
  
  def test_read_coord_2D(self):
    import numpy
    for mesh_name in self.valid_mesh_2D:
      mesh=self.fh.get_mesh(mesh_name)
      x=mesh.read_coordinates(0,1024)
      y=mesh.read_coordinates(1,1024)
  
  def test_read_coord_3D(self):
    import numpy
    for mesh_name in self.valid_mesh_3D:
      mesh=self.fh.get_mesh(mesh_name)
      x=mesh.read_coordinates(0,1024)
      y=mesh.read_coordinates(1,1024)
      z=mesh.read_coordinates(2,1024)

  def test_read_fails(self):
    import numpy
    mesh=self.fh.get_mesh(self.valid_mesh_2D[0])
    try:
      mesh.read_coordinates(3,1024)
    except:
      print 'Caught invalid dimension request'
    try:
      mesh.read_coordinates(0,5)
    except:
      print 'Caught invalid nnodes size'

  def test_coordinates_1D(self):
    import numpy
    mesh=self.fh.get_mesh(self.valid_mesh_1D[0])
    coordinates=mesh.coordinates()
    self.assertEqual(len(coordinates),1)
    nnodes=coordinates[0].size
    idx=0
    for array in coordinates:
      read_array=mesh.read_coordinates(idx,nnodes)
      self.assertEqual(read_array.all(),array.all())

  def test_coordinates_2D(self):
    import numpy
    mesh=self.fh.get_mesh(self.valid_mesh_2D[1])
    coordinates=mesh.coordinates()
    self.assertEqual(len(coordinates),2)
    nnodes=coordinates[0].size
    idx=0
    for array in coordinates:
      read_array=mesh.read_coordinates(idx,nnodes)
      self.assertEqual(read_array.all(),array.all())

  def test_coordinates_3D(self):
    import numpy
    mesh=self.fh.get_mesh(self.valid_mesh_3D[1])
    coordinates=mesh.coordinates()
    self.assertEqual(len(coordinates),3)
    nnodes=coordinates[0].size
    idx=0
    for array in coordinates:
      read_array=mesh.read_coordinates(idx,nnodes)
      self.assertEqual(read_array.all(),array.all())

class TestConnectivity(unittest.TestCase):

  def setUp(self):
    import os
    import numpy
    import random
    from Danu import Output
    from Danu import UNSTRUCTURED_MESH
    from Danu import HEX_ELEM, HEX_ELEM_ORDER 

    self.filename = 'test-Mesh.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')

    self.mesh_name='Test Mesh 3D HEX'
    self.mesh=self.fh.add_unstruct_mesh(self.mesh_name,HEX_ELEM)
    
    self.data_name='Data to Read'
    self.nelem=random.randint(10,2048)
    self.data=numpy.zeros((self.nelem,HEX_ELEM_ORDER),dtype=numpy.int32)
    nc=0
    while nc < self.nelem:
      i=0
      while i < HEX_ELEM_ORDER:
        self.data[nc][i]=random.randint(0,100000)
        i=i+1
      nc=nc+1	

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      #os.remove(self.filename)

  def test_basic(self):
    import numpy
    import random

    # try to read before a write
    try:
      r_data=self.mesh.read_connectivity()
    except:
      print 'Caught the read before write error'

    self.mesh.write_connectivity(self.data)  
    read_data=self.mesh.read_connectivity()
    self.assertEqual(read_data.all(),self.data.all())

  def test_offset(self):
    import numpy
    import random

    self.mesh.write_connectivity(self.data)  

    offset=self.mesh.connectivity_offset()
    self.assertEqual(offset,0)


class TestAttributes(unittest.TestCase):

  def setUp(self):
    import os
    import numpy
    import random
    from Danu import Output
    from Danu import UNSTRUCTURED_MESH
    from Danu import HEX_ELEM, HEX_ELEM_ORDER 

    self.filename = 'test-Mesh.h5'
    if  os.path.exists(self.filename):
      os.remove(self.filename)

    self.fh=Output(self.filename,'w')

    self.mesh_name='Test Mesh 3D HEX'
    self.mesh=self.fh.add_unstruct_mesh(self.mesh_name,HEX_ELEM)
    
    self.n=random.randint(10,2048)
    self.x=numpy.random.random_sample((self.n))
    self.y=numpy.random.random_sample((self.n))
    self.z=numpy.random.random_sample((self.n))
    self.coordinates=[self.x,self.y,self.z]
    self.mesh.write_coordinates(self.x,self.y,self.z)
    self.nelem=random.randint(1024,2048)
    self.data=numpy.zeros((self.nelem,HEX_ELEM_ORDER),dtype=numpy.int32)
    nc=0
    while nc < self.nelem:
      i=0
      while i < HEX_ELEM_ORDER:
        self.data[nc][i]=random.randint(0,100000)
        i=i+1
      nc=nc+1	
    self.mesh.write_connectivity(self.data)

  def tearDown(self):
    import os
    if os.path.exists(self.filename):
      self.fh.close()
      os.remove(self.filename)

  def runTest(self):
    import Danu
    self.assertEqual(self.x.size,self.mesh.nnodes())
    self.assertEqual(Danu.HEX_ELEM_ORDER,self.mesh.elem_order())
    print self.mesh.nelem()



      
if __name__ == '__main__':
   unittest.main()
   #suite=unittest.TestSuite()
   #suite.addTest(TestConnectivity('runTest'))
   #unittest.TextTestRunner().run(suite)





    


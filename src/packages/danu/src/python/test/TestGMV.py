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

import numpy
import GMV

class SimpleMesh:

  def __init__(self,x0,x1,y0,y1,z0,z1,nx,ny,nz):
    
    self.num_cells=nx*ny*nz
    self.num_nodes=(nx+1)*(ny+1)*(nz+1)
    self.num_faces=(nx+1)*ny*nz + nx*(ny+1)*nz + nx*ny*(nz+1)
    self.x_min=x0
    self.x_max=x1
    self.y_min=y0
    self.y_may=y1
    self.z_min=z0
    self.z_max=z1

    self.nx=nx
    self.ny=ny
    self.nz=nz
    
    dx=(x1-x0)/nx
    dy=(y1-y0)/ny
    dz=(z1-z0)/nz

    # Create the nodal coordinates
    self.x=numpy.float64(numpy.zeros(self.num_nodes))
    self.y=numpy.float64(numpy.zeros(self.num_nodes))
    self.z=numpy.float64(numpy.zeros(self.num_nodes))
    k=0
    while k <= nz:
      j=0
      while j <= ny:
        i=0
        while i <= nx:
	  inode = i + j*(nx+1) + k*(nx+1)*(ny+1)
	  self.x[inode] = x0 + i*dx
	  self.y[inode] = y0 + j*dy
	  self.z[inode] = z0 + k*dx
	  mystring="inode=%d %d,%d,%d=%f,%f,%f" % (inode,i,j,k,self.x[inode],self.y[inode],self.z[inode])
	  #print mystring
	  i=i+1
        j=j+1
      k=k+1

    # Create the connectivity -- dummy data for now
    self.connectivity=numpy.zeros((self.num_cells,8),dtype=numpy.int32)
    nc=0
    while nc < self.num_cells:
      i=0
      while i < 8:
	self.connectivity[nc][i] = nc
	i=i+1
      nc=nc+1	

  def node_index(self,i,j,k): 
    return i + j*(self.nx+1) + k*(nx+1)*(ny+1)

  def cell_index(self,i,j,k):
    return i + j*self.nx + k*self.nx*self.ny



class TestMeshFile(unittest.TestCase):

  def setUp(self):
    self.filename         = 'test-gmv-meshfile.gmv'

    self.mesh=SimpleMesh(0.0,1.0,0,1.0,0.0,1.0,4,1,1)


  def test_basic(self):
    file = GMV.MeshFile(self.filename)
    del file

  def test_active_file(self):
    file = GMV.MeshFile(self.filename)
    file.open()
    try:
      other_file=GMV.MeshFile('this-will-fail-'+self.filename)
    except:
      print 'Caught the exception with more than one file open'
    file.close()  

  def test_write(self):
    file = GMV.MeshFile(self.filename)
    file.open()
    file.write(self.mesh.x,self.mesh.y,self.mesh.z,self.mesh.connectivity)
    file.close()


class TestDataFile(unittest.TestCase):

  def setUp(self):

    self.filename        = 'test-gmv-data.gmv'

    self.meshname = 'test-gmv-mesh.gmv'
    self.mesh=SimpleMesh(0.0,1.0,0,1.0,0.0,1.0,2,2,2)
    self.meshfile=GMV.MeshFile(self.meshname)
    self.meshfile.open()
    self.meshfile.write(self.mesh.x,self.mesh.y,self.mesh.z,self.mesh.connectivity)
    self.meshfile.close()

  def test_basic(self):
    
    file=GMV.DataFile(self.filename,self.meshfile)
    del file

  def test_create(self):
 
    file=GMV.DataFile(self.filename,self.meshfile)
    file.open()
    file.close()

  def test_create_time(self):
 
    file=GMV.DataFile(self.filename,self.meshfile,0,0.0)
    file.open()
    file.close()

  def test_write_cell(self):

    file=GMV.DataFile(self.filename,self.meshfile)
    dummy_cell=numpy.float64(numpy.zeros(self.mesh.num_cells))
    file.open()
    file.write_cell_data('dummy_cell_data',dummy_cell)
    file.close()

  def test_write_node(self):

    file=GMV.DataFile(self.filename,self.meshfile)

    dummy_node=numpy.float64(numpy.zeros(self.mesh.num_nodes))
    file.open()
    file.write_node_data('dummy_node_data',dummy_node)
    file.close()
    

if __name__ == '__main__':
    
  unittest.main()
  #mesh=SimpleMesh(0.0,1.0,0,1.0,0.0,1.0,4,1,1)
  #print mesh.cell_index(1,1,1)
  #print mesh.x[1]
  #print mesh.y[1]
  #print mesh.z[1]
  #suite=unittest.TestSuite()
  #suite.addTest(TestDataFile('test_basic'))
  #result=unittest.TextTestRunner().run(suite)
  #if result.wasSuccessful():
  #  print 'Test suite passed'
  #  sys.exit(0)
  #else:
  #  sys.exit(1)



    


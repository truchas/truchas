#!/usr/bin/env python
#
#
#
#
#  Usage:
#
#    xdmf-parser.py [options] DANU_FILE
#
#    Creates light data model described in file DANU_FILE.xml
#     of the mesh associated data in DATA_FILE
#
#
# --- Standard Python Modules
import sys
import os

import numpy

from xml.etree import ElementTree

# --- Attempt to import Danu. Exit gracefully
try:
  import Danu
except ImportError:
  danu_py_install='@Danu_Python_INSTALL_DIR@'
  sys.path.append(danu_py_install)
  try:
    import Danu
  except ImportError:
    print "Attempted to add %s to Python to import Danu. Failed" % (danu_py_install)
    raise


import Xdmf

def setSequenceFilename(id,sim_name,h5file_name):
  base_name=os.path.splitext(h5file_name)[0]
  f='%s-%s-%d.xmf'%(base_name,sim_name,id)
  return f

def setCollectionFilename(sim_name,h5file_name):
  base_name=os.path.splitext(h5file_name)[0]
  f='%s-%s.Collection.xmf'%(base_name,sim_name)
  return f

def main(filename):

  # Open the file for reading
  f=Danu.Output(filename)

  # Open simulation (MAIN ONLY for now)
  simulation_group_name='/Simulations'+'/MAIN'
  sim_name='MAIN'
  sim = f.get_simulation(sim_name)

  # Xdmf model
  model = Xdmf.XdmfModel()
  grid=model.add_grid('MESH')
  grid.set_type('Uniform')

  # Open the mesh
  mesh_group_name=simulation_group_name+ '/Mesh'
  mesh = sim.open_mesh_link()
  dim = mesh.dim
  nnodes = mesh.nnodes()
  nelem  = mesh.nelem()
  elem_type = mesh.elem_type
  elem_order = mesh.elem_order()
  offset = mesh.connectivity_offset()


  # Coordinates
  coordinate_dataset_name=mesh_group_name+'/Nodal Coordinates'
  geo = Xdmf.XdmfGeometry('geo')
  geo.set_type('XYZ')
  grid.append(geo)
  data = Xdmf.XdmfDataItem('Coordinates')
  data.set_num_type('Float')
  data.set_precision(8)
  data.set_format('HDF')
  data.set_dimensions([nnodes,dim])
  loc=Xdmf.munge_hdf5_dataset(filename,coordinate_dataset_name)
  data.set_text(loc)
  geo.append(data)

  # Connectivity
  connect_dataset_name=mesh_group_name+'/Element Connectivity'
  type=Xdmf.xdmf_mesh_type(elem_type)
  top=Xdmf.XdmfTopology('topo',type)
  top.set_offset(offset)
  top.set_nelem(nelem)
  grid.append(top)
  data = Xdmf.XdmfDataItem('Connectivity')
  data.set_num_type('Int')
  data.set_precision(4)
  data.set_format('HDF')
  data.set_dimensions([nelem,elem_order])
  loc=Xdmf.munge_hdf5_dataset(filename,connect_dataset_name)
  data.set_text(loc)
  top.append(data)

  # Loop through the sequences
  time_files=[]
  series_group_name=simulation_group_name+'/Series Data'
  id=1
  while id <= sim.sequence_count():
    model_file=setSequenceFilename(id,sim_name,filename)
    print 'Creating file: %s'%(model_file)
    fh=open(model_file,'w')
    seq_name=sim.get_sequence_name(id)
    seq=sim.get_sequence(seq_name)
    grid.set_time(seq.time)
    # Now loop through each dataset
    for dname in seq.data_list():
      print '\tAdding %s to %s'%(dname,model_file)
      try:
        d=Xdmf.XdmfDanuDataset(seq,dname)
      except:
        msg='\t\tFailed to read data set %s will skip'%(dname)
        print msg
      else:
        data_loc='%s/Series %d/%s'%(series_group_name,id,dname)
        d.set_data_location(Xdmf.munge_hdf5_dataset(filename,data_loc))
        grid.append(d)
    model.write(fh)
    fh.close()
    time_files.append(model_file)
    print
    id=id+1

  # Now loop through the time step files to create the movie file
  movie_file=setCollectionFilename(sim_name,filename)
  mf=open(movie_file,'w')
  movie=Xdmf.XdmfTemporalCollection()
  xpath='/Xdmf/Domain/Grid'
  for f in time_files:
    print 'Adding %s to the movie file:%s'%(f,movie_file)
    movie.add_file_ref(f,xpath)
  movie.write(mf)
  mf.close()


if __name__ == '__main__':
  try:
    main(sys.argv[1])
  except IndexError:
    print 'Usage: <parser> HDF5 file'
    raise







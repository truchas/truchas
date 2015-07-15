#!/usr/bin/env python
# ############################################################################ #
#
#  Truchas GMV Parser
#
# ############################################################################ #

# ----- Import modules
import sys
import os
from optparse import OptionParser as OptParse

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

# --- Attempt to import GMV. Exit gracefully
try:
  import GMV
except ImportError:
  gmv_py_install='@Danu_Python_INSTALL_DIR@'
  sys.path.append(gmv_py_install)
  try:
    import GMV
  except ImportError:
    print "Attempted to add %s to Python to import GMV. Failed" % (danu_py_install)
    raise

# --- Process command line options
class CLError(Exception):
  def __init__(self,msg):
    Exception.__init__(self)
    self.msg=msg
  def __str__(self):
    return self.msg

class CL_Options:

  _usage = '''Usage: %prog [options] DANU_FILE

  Truchas GMV Parser

  '''

  def __init__(self):
    self.parser=OptParse(usage=CL_Options._usage,add_help_option=True)
    p=self.parser
    p.add_option("-s", "--simulation", dest="sim_name", action="store",
                  default="MAIN",
		  help="Danu simulation name",
		  metavar="STRING")
    p.add_option("-m", "--mesh", dest="mesh_name", action="store",
                  default="DEFAULT",
		  help="Danu mesh name",
		  metavar="STRING")
    p.add_option("-b", "--binary", dest="binary", action="store_true",
                  default=False,
		  help="Write GMV files in a binary format")
    p.add_option("--base", dest="base_name", action="store",
	          default=None,
                  help="Use this string instead of the simulation name for file names",
		  metavar="STRING")
    p.add_option("--stride", dest="stride", type="int", action="store",
                  default=1,
		  help="GMV data file suffix stride",
		  metavar="NUM")
    p.add_option("--start", dest="start", type="int", action="store",
                  default=0,
		  help="GMV data file initial suffix start value",
		  metavar="NUM")

  def parse(self):
    (options,args)=self.parser.parse_args()

    try:
      infile=args[0]
      if not os.path.isfile(infile):
	raise CLError('Can not locate Danu file %s' % infile)
    except IndexError:
      raise CLError('Missing Danu file')
    except CLError, e:
      self.parser.print_usage()
      raise

    return (options,infile)

# ---------------------------------------------------------------------------
#
# Truchas Danu Class
#
# Someday this class will  grow up and
# live as a separate module in the Truchas area
#
# ---------------------------------------------------------------------------

# From a Google search on creating nested dicts
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def build_series_name(num):
  return 'Series %d' % (num)

class TruchasDanuFile:

  def __init__(self,filename):

    # Initialize the file open
    self.filename=filename
    try:
      self.file=Danu.Output(filename,'r')
    except RuntimeError:
      raise RuntimeError, 'Failed to open %s' % (filename)
    except:
      raise

    # Create some useful arrays instead of rereading the file
    self.simulations=self.file.simulation_list()
    self.meshes=self.file.mesh_list()

    # Create useful dictionaries....del sim object to close the HDF5 object
    # self.datasets[sim_name]=array of dataset names
    # self.series[sim_name][series name]['datasets']=array of dataset names
    # self.series[sim_name][series name]['time']= problem time
    # self.series[sim_name][series name]['cycle']= problem cycle
    self.series=AutoVivification()
    self.datasets={}
    for sim_name in self.simulations:
      sim=self.file.get_simulation(sim_name)
      self.datasets[sim_name]=sim.data_list()

      for group in sim.sequence_list():
	seq=sim.get_sequence(group)
	self.series[sim_name][group]['datasets']=seq.data_list()
	self.series[sim_name][group]['time']=seq.time
	self.series[sim_name][group]['cycle']=seq.cycle
	del seq
      del sim

  # For debugging
  def pretty_print(self):
    from operator import attrgetter

    for mesh in self.meshes:
      print 'Mesh: %s' % (mesh)

    for sim in self.simulations:
      print 'Simulation: %s' % (sim)
      print
      print '\tDatasets'
      for ds in self.datasets[sim]:
	print '\t\t%s' % (ds)
      print
      print '\tSeries'
      print self.series[sim]
      num_names=range(1,1+len(self.series[sim]))
      ordered_series=map(build_series_name,num_names)
      for g in ordered_series:
	print '\t\t%s t=%1.5e cycle=%d' % (g,self.series[sim][g]['time'],self.series[sim][g]['cycle'])
	print '\t\t\tDatasets: ' + str(self.series[sim][g]['datasets'])

  def series_data_search(self,sim_name,series_name,key_name,key_value):
    return_datasets=[]
    try:
      datasets=self.series[sim_name][series_name]['datasets']
    except:
      print 'Failed to define datasets'
    else:
      seq=self.file.get_simulation(sim_name).get_sequence(series_name)
      for data in datasets:
	if seq.data_attribute_exists(data,key_name):
	  value=seq.get_data_attribute(data,key_name)
	  if value == key_value:
	    return_datasets.append(data)
      del seq
    return return_datasets

  def get_series_data_ndim(self,sim_name,series_name,data_name):
    ndim=-1
    datasets=[]
    try:
      datasets=self.series[sim_name][series_name]['datasets']
    except:
      print 'Failed to define dataset'
    else:
      if data_name in datasets:
	seq=self.file.get_simulation(sim_name).get_sequence(series_name)
	ndim=seq.get_data_ndims(data_name)
	del seq
    return ndim

  def get_cell_fields(self,sim_name,series_name):
    key_name='FIELDTYPE'
    key_value='CELL'
    return self.series_data_search(sim_name,series_name,key_name,key_value)

  def get_cell_fields_scalars(self,sim_name,series_name):
    scalar_fields=[]
    seq=self.file.get_simulation(sim_name).get_sequence(series_name)
    for field in self.get_cell_fields(sim_name,series_name):
      if seq.get_data_ndims(field) == 1:
	scalar_fields.append(field)
    del seq
    return scalar_fields

  def get_cell_fields_vectors(self,sim_name,series_name):
    vector_fields=[]
    seq=self.file.get_simulation(sim_name).get_sequence(series_name)
    for field in self.get_cell_fields(sim_name,series_name):
      if seq.get_data_ndims(field) > 1 and field != 'VOF':
	vector_fields.append(field)
    del seq
    return vector_fields


  def get_node_fields(self,sim_name,series_name):
    key_name='FIELDTYPE'
    key_value='NODE'
    return self.series_data_search(sim_name,series_name,key_name,key_value)

  def get_node_fields_scalars(self,sim_name,series_name):
    scalar_fields=[]
    seq=self.file.get_simulation(sim_name).get_sequence(series_name)
    for field in self.get_node_fields(sim_name,series_name):
      if seq.get_data_ndims(field) == 1:
	scalar_fields.append(field)
    del seq
    return scalar_fields

  def get_node_fields_vectors(self,sim_name,series_name):
    vector_fields=[]
    seq=self.file.get_simulation(sim_name).get_sequence(series_name)
    for field in self.get_node_fields(sim_name,series_name):
      if seq.get_data_ndims(field) > 1:
	vector_fields.append(field)
    del seq
    return vector_fields

  def get_field_names(self,sim_name,series_name,field,ncomp):
    names=[]
    seq=self.file.get_simulation(sim_name).get_sequence(series_name)
    field_attrs=seq.data_attributes(field)
    if field_attrs.has_key('FIELDTYPE'):
      i=1
      while i <= ncomp:
	fname_key='FIELDNAME' + str(i)
	if field_attrs.has_key(fname_key):
	  names.append(seq.get_data_attribute(field,fname_key))
	i=i+1
    return names





cl_opts=CL_Options()
(options,danu_infile)=cl_opts.parse()


# ----- Create Danu objects

# --- File
try:
  truchas_file=TruchasDanuFile(danu_infile)
except:
  print 'Failed to open: %s' % (danu_infile)
  raise

# --- Simulation
try:
  sim=truchas_file.file.get_simulation(options.sim_name)
except:
  print 'Failed to open Simualtion %s in %s' % (options.sim_name,danu_infile)
  raise

# --- Mesh
try:
  mesh=truchas_file.file.get_mesh(options.mesh_name)
except:
  print 'Failed to open mesh %s in %s' % (options.mesh_name,danu_infile)
  raise

# ----- Create the GMV mesh file
# TypeError rasied if base_name is None
try:
  gmv_meshname=options.base_name+'-mesh'+'.gmv'
except TypeError:
  gmv_meshname=options.sim_name+'-mesh'+'.gmv'
except:
  raise
gmv_mesh=GMV.MeshFile(gmv_meshname)

print 'Creating mesh file %s for Simulation %s' % (gmv_meshname,options.sim_name)
coordinates=mesh.coordinates()
connectivity=mesh.read_connectivity()
gmv_mesh.open(options.binary)
gmv_mesh.write(coordinates[0],coordinates[1],coordinates[2],connectivity)
gmv_mesh.close()

# ----- Define useful datasets from the simulation

# --- Number of processors
#     Danu returns a int32 array...want a Python int
numprocs=0
try:
  npes=sim.data_read('NUMPROCS')
except RuntmeError:
  print 'NUMPROCS not found'
except:
  raise
else:
  numprocs=int(npes[0])

i=1
proc_names=[]
while i <= numprocs:
  proc_name='%s%04d'%('PE',i)
  proc_names.append(proc_name)
  i=i+1


# --- The cell map rank-1 integer array
cellmap=[]
try:
  cellmap=sim.data_read('CELLMAP')
except RuntmeError:
  print 'CELLMAP not found'
except:
  raise

# --- The node id map rank-1 integer array
nodemap=[]
try:
  nodemap=sim.data_read('NODEMAP')
except RuntmeError:
  print 'NODEMAP not found'
except:
  raise

# --- The material ids rank-1 integer array
blockids=[]
materials=[]
nmats=0
try:
  blockids=sim.data_read('BLOCKID')
except RuntimeError:
  print 'BLOCKID not found'
except:
  raise
else:
  unique_ids = list(set(blockids))
  unique_ids.sort()
  nmats=len(unique_ids)
  id_to_index = {}
  for i in range(nmats):
    id_to_index[unique_ids[i]] = i+1  # want 1-based
  for i in range(len(blockids)):
    blockids[i] = id_to_index[blockids[i]]
  for i in unique_ids:
    mat_name='blk_%d'%(i)
    materials.append(mat_name)
  print 'Number of materials: ' + str(nmats)
  print 'Materials: ' + str(materials)

# --- The cell partition array
cellpart=[]
if numprocs > 1:
  try:
    cellpart=sim.data_read('CELLPART')
  except RuntmeError:
    print 'CELLPART not found'
  except:
    raise

# --- The node partition array
nodepart=[]
if numprocs > 1:
  try:
    nodepart=sim.data_read('NODEPART')
  except RuntmeError:
    print 'NODEPART not found'
  except:
    raise

# ----- Generate the GMV data file

# --- Define the base file name
# TypeError is raised if base_name is None
try:
  data_file_base=options.base_name+'-data'+'.gmv'
except TypeError:
  data_file_base=options.sim_name+'-data'+'.gmv'
except:
  raise

# --- Loop through each time series and generate a data file
cnt=options.start
num_names=range(1,1+len(truchas_file.series[options.sim_name]))
ordered_groups=map(build_series_name,num_names)
for group in ordered_groups:

  # Find the cell based data fields
  cell_fields=AutoVivification()
  cell_fields['scalars']=truchas_file.get_cell_fields_scalars(options.sim_name,group)
  cell_fields['vectors']=truchas_file.get_cell_fields_vectors(options.sim_name,group)

  # Find the node based data fields
  node_fields=AutoVivification()
  node_fields['scalars']=truchas_file.get_node_fields_scalars(options.sim_name,group)
  node_fields['vectors']=truchas_file.get_node_fields_vectors(options.sim_name,group)

  # Open the sequence HDF5 group
  seq=sim.get_sequence(group)

  # Create the data file
  gmv_data_filename='%s.%04d'%(data_file_base,cnt)
  print '\tWriting file ' + gmv_data_filename
  gmv_data=GMV.DataFile(gmv_data_filename,gmv_mesh,seq.cycle,seq.time)
  gmv_data.open(options.binary)

  # Write BLOCKID -- GMV material block
  if len(blockids) > 0:
    gmv_data.write_materials(nmats,materials,blockids)

  # Write CELLMAP -- GMV cellids
  if len(cellmap) > 0:
    gmv_data.write_cellids(cellmap)

  # Partition maps -- GMV flag blocks
  if ( ( len(cellpart) > 0 ) or ( len(nodepart) > 0 ) ):
    gmv_data.begin_flag_block()

  # Write CELLPART
  if len(cellpart) > 0:
    gmv_data.write_cell_flag('cellpart',proc_names,numprocs,cellpart)

  # Write NODEPART
  if len(nodepart) > 0:
    gmv_data.write_node_flag('nodepart',proc_names,numprocs,nodepart)

  if ( ( len(cellpart) > 0 ) or ( len(nodepart) > 0 ) ):
    gmv_data.end_flag_block()

  # Write the cell and node based scalars
  gmv_data.begin_variable_block()
  for field in cell_fields['scalars']:
    #print '\tWriting %s to %s'%(field,gmv_data_filename)
    data_name=seq.get_data_attribute(field,'FIELDNAME')
    data=seq.data_read(field)
    gmv_data.write_cell_data(data_name,data)

  # VOF Field special cell field case that is stored as a rank 2 array, but should
  # be dumped out as scalar data for each component
  if seq.data_exists('VOF'):
    vof_data=seq.data_read('VOF').transpose()
    (ncomp,ncells)=vof_data.shape
    vof_names=truchas_file.get_field_names(options.sim_name,group,'VOF',ncomp)
    i=0
    while i < ncomp:
      gmv_data.write_cell_data(vof_names[i],vof_data[i,:])
      i=i+1

  for field in node_fields['scalars']:
    #print '\tWriting %s to %s'%(field,gmv_data_filename)
    data_name=seq.get_data_attribute(field,'FIELDNAME')
    data=seq.data_read(field)
    gmv_data.write_node_data(data_name,data)

  gmv_data.end_variable_block()

  # Write the cell and node based vectors
  gmv_data.begin_vector_block()
  for field in cell_fields['vectors']:
    #print '\tWriting %s to %s'%(field,gmv_data_filename)
    data=seq.data_read(field).transpose()
    (ncomp,ncells)=data.shape
    data_names=truchas_file.get_field_names(options.sim_name,group,field,ncomp)
    gmv_data.write_cell_data_vectors(field,data_names,data)

  for field in node_fields['vectors']:
    #print '\tWriting %s to %s'%(field,gmv_data_filename)
    data=seq.data_read(field).transpose()
    (ncomp,nnodes)=data.shape
    print 'ncomp=%d nnodes=%d\n'%(ncomp,nnodes)
    data_names=truchas_file.get_field_names(options.sim_name,group,field,ncomp)
    gmv_data.write_node_data_vectors(field,data_names,data)
  gmv_data.end_vector_block()

  # Close the GMV data file
  gmv_data.close()

  # Close the HDF5 sequence group
  del seq

  # Update the counter
  cnt=cnt+options.stride


#
# And exit .....
sys.exit()

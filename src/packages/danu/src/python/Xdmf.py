# ##############################################################################
#
#  XDMF Module
#    http://wwww.xdmf.org
#
# ##############################################################################

import sys
import os

import numpy

import xml.etree.ElementTree as ElementTree

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

# --- Global defines
def isPyFloat(num):

  if isinstance(num,float):
    return True
  else:
    return False

def isPyInt(num):
  if isinstance(num,int):
    return True
  else:
    return False

def xdmf_precision(inp):
  in_str=''
  if isinstance(inp,numpy.dtype):
    in_str=inp.dtype.name
  elif isinstance(inp,str):
    in_str=inp
  else:
    raise 'Unknown input type in xdmf_precision'

  if in_str == 'int64':
    return  8
  elif in_str == 'int32':
    return 4
  elif in_str == 'float' or in_str == 'float64':
    return 8
  elif in_str == 'float32':
    return 4
  else:
    raise 'Invliad numpy data type in xdmf_precision'

def xdmf_num_type(inp):

  in_str=''
  if isinstance(inp,numpy.dtype):
    in_str=inp.dtype.name
  elif isinstance(inp,str):
    in_str=inp
  else:
    raise 'Unknown input type in xdmf_precision'

  if in_str == 'int64' or in_str == 'int32':
    return  'Int'
  elif in_str == 'float' or 'in_str' == 'float64' or in_str == 'float32':
    return 'Float'
  else:
    raise 'Invliad numpy data type in xdmf_num_type'

def xdmf_mesh_type(elem_type):
  r_elem_type=''
  if elem_type == Danu.HEX_ELEM:
    r_elem_type='Hexahedron'
  elif elem_type == Danu.QUAD_ELEM:
    r_elem_type='Quadrilateral'
  elif elem_type == Danu.TET_ELEM:
    r_elem_type='Tetrahedron'
  elif elem_type == Danu.TRI_ELEM:
    r_elem_type='Triangle'
  elif elem_type == Danu.LINE_ELEM:
    r_elem_type='Polyline'
  else:
    raise 'Invalid elemant type'
  return r_elem_type

def munge_hdf5_dataset(file,loc):
  return file + ':' + loc

def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


'''
Basic Xdmf XML Element container
'''
class XdmfElement(object):

  # XML tag
  tag = None

  def __init__(self,attrib={}):

    # Element
    for k in attrib:
      attrib[k]=str(attrib[k])
    self.elem = ElementTree.Element(tag=self.tag,attrib=attrib)

  def tag(self):
    return self.elem.tag

  def attributes(self):
    return self.attrib

  def get_elem_object(self):
    return self.elem

  def set_text(self,data):
    self.elem.text = str(data)
    return self.elem.text

  def get_text(self):
    return self.elem.text

  def set_attribute(self,key,value):
    self.elem.set(key,str(value))
    return self.elem.attrib[key]

  def get_attribute(self,key):
    return self.elem.get(key)

  def append(self,element):
    return self.elem.append(element.elem)

  def insert(self,index,element):
    return self.elem.insert(index,element)

  def remove(self,element):
    return self.elem.remove(element)


'''
Base class for XDMF elements that have a
'Name' attribute
'''
class XdmfElementName(XdmfElement):

  def __init__(self,name,attr={}):
    attr['Name']=name
    XdmfElement.__init__(self,attr)

  def set_name(self,name):
    return self.set_attribute('Name',name)

  def get_name(self):
    return self.get_attribute('Name')

'''
XDMF Attribute Element
'''
class XdmfAttribute(XdmfElementName):

  tag = 'Attribute'

  def __init__(self,name,type='Scalar',center='Node'):
    attr={'Type' : type, 'Center' : center}
    XdmfElementName.__init__(self,name,attr)

  def set_type(self,type):
    return self.set_attribute('Type',type)

  def get_type(self):
    return self.get_attribute('Type')

  def set_center(self,center):
    return self.set_attribute('Center',center)

  def get_center(self):
    return self.get_attribute('Center')

'''
XDMF DataItem
'''
class XdmfDataItem(XdmfElementName):

    tag = 'DataItem'

    def __init__(self,name,num_type='Float',precision=4,format='XML'):
      attr={'DataType' : num_type, 'Precision' : precision, 'Format' : format}
      XdmfElementName.__init__(self,name,attr)

    def set_type(self,type):
      return self.set_attribute('ItemType',type)

    def get_type(self):
      return self.get_attribute('ItemType')

    def set_dimensions(self,dims):
      try:
        l=len(dims)
      except TypeError:
        string=str(dims)
      except:
        raise
      else:
        string=''
        i=0
        while i < l-1:
          string+='%s '%(str(dims[i]))
          i=i+1
        string+='%s'%(str(dims[-1]))
      return self.set_attribute('Dimensions',string)

    def get_dimensions(self):
      l=self.get_attribute('Dimensions').split()
      i=0
      while i < len(l):
        l[i] = int(l[i])
        i=i+1
      if len(l) > 1:
        return l
      elif len(l) == 1:
        return l[0]
      else:
        return []
      return list

    def set_num_type(self,num_type):
      return self.set_attribute('DataType',num_type)

    def get_num_type(self):
      return self.get_attribute('DataType')

    def set_precision(self,num):
      return self.set_attribute('Precision',num)

    def get_precision(self):
      return int(self.get_attribute('Precision'))

    def set_format(self,format):
      return self.set_attribute('Format',format)

    def get_format(self):
      return self.get_attribute('Format')

'''
XDMF Domain Element
'''
class XdmfDomain(XdmfElementName):

    tag = 'Domain'


'''
XDMF Geometry Element
'''
class XdmfGeometry(XdmfElementName):

    tag = 'Geometry'

    def __init__(self,name,type='XYZ'):
      attr={}
      attr['Type'] = str(type)
      XdmfElementName.__init__(self,name,attr)

    def set_type(self,type):
      return self.set_attribute('Type',type)

    def get_type(self):
      return self.get_attribute('Type')

'''
XDMF Include Element
'''
class XdmfInclude(XdmfElement):

  tag = 'xi:include'

  def set_href(self,string):
    return self.set_attribute('href',string)

  def get_href(self):
    return self.get_attribute('href')

  def set_xpointer_path(self,path):
    xml_string='xpointer(/%s)'%(path)
    return self.set_attribute('xpointer',xml_string)

  def get_xpointer_path(self):
    return self.get_attribute('xpointer')

'''
XDMF Information Element
'''
class XdmfInformation(XdmfElementName):

  tag = 'Information'

  def set_value(self,value):
    return self.set_attribute('Value',value)

  def get_value(self):
    return self.get_attribute('Value')

'''
XDMF Root Element
'''
class XdmfRoot(XdmfElement):

  tag = 'Xdmf'

  def set_version(self,version):
    return self.set_attribute(self,'Version',version)

  def get_version(self):
    return self.get_attribute(self,'Version')


'''
XDMF Topology Element
'''
class XdmfTopology(XdmfElementName):

  tag = 'Topology'

  def __init__(self,name,type,attr={}):
    attr['TopologyType'] = str(type)
    attr['BaseOffset'] = str(0)
    XdmfElementName.__init__(self,name,attr)

  def set_type(self,type):
    return self.set_attribute('Topology',type)

  def get_type(self):
    return self.get_attribute('Topology',type)

  def set_offset(self,num):
    return self.set_attribute('BaseOffset',num)

  def get_offset(self):
    return int(self.get_attribute('BaseOffset'))

  def set_order(self,num):
    return self.set_attribute('Order',num)

  def get_order(self):
    return int(self.get_attribute('Order'))

  def set_npelem(self,num):
    return self.set_attribute('NodesPerElement',num)

  def get_npelem(self):
    return int(self.get_attribute('NodesPerElement'))

  def set_nelem(self,num):
    return self.set_attribute('NumberOfElements',num)

  def get_nelem(self):
    return int(self.get_attribute('NumberOfElements'))

'''
XDMF Time Element
'''
class XdmfTime(XdmfElement):

  tag = 'Time'

  def __init__(self,data=None):
    attr={}
    if isPyFloat(data):
      attr['TimeType']='Single'
      attr['Value'] = str(data)
    XdmfElement.__init__(self,attr)

  def set_value(self,t):
    self.set_attribute('TimeType','Single')
    return self.set_attribute('Value',t)

  def get_value(self):
    return float(self.get_attribute('Value'))


'''
XDMF Grid
'''

class XdmfGrid(XdmfElementName):

  tag = 'Grid'

  def set_time(self,t):
    try:
      time=self.time
    except AttributeError:
      time=XdmfTime(t)
      self.time=time
      self.append(time)
    self.time.set_value(t)

  def get_time(self):
    t=None
    try:
      t=self.time.get_value()
    except AttributeError:
      print 'Time is not set for this grid'
    return t

  def set_type(self,type):
    return self.set_attribute('GridType',type)

  def get_type(self):
    return self.get_attribute('GridType')

'''
XdmfTree -- Base class for the XML trees
'''
class XdmfTree(object):

  def __init__(self):
    self.root = XdmfRoot()

  def write(self,fh=sys.stdout):
    root_elem=self.root.get_elem_object()
    indent(root_elem)
    ElementTree.ElementTree(root_elem).write(fh)

  def append(self,xdmf_elem):
    if isinstance(xdmf_elem,XdmfElement):
      self.root.append(xdmf_elem)
    else:
      raise 'Invalid argument in append_elem'


'''
XdmfModel
'''
class XdmfModel(XdmfTree):

  def __init__(self,name=None):
    XdmfTree.__init__(self)

  def add_grid(self,name):
    d=self.add_domain('dummy')
    g=XdmfGrid(name)
    d.append(g)
    return g

  def add_domain(self,name):
    d=XdmfDomain(name)
    self.root.append(d)
    return d


'''
XdmfTemporalCollection
'''
class XdmfTemporalCollection(XdmfTree):

  def __init__(self):
    XdmfTree.__init__(self)
    self.root.set_attribute('Version','2.0')
    self.root.set_attribute('xmlns:xi','http://www.w3.org/2001/XInclude')
    self.domain=XdmfDomain(name='Test')
    self.append(self.domain)
    g_attr={"CollectionType" : "Temporal", "GridType" : "Collection"}
    self.grid=XdmfGrid('Time Steps',g_attr)
    self.domain.append(self.grid)

  def add_file_ref(self,h5file,path):
    x=XdmfInclude()
    x.set_href(h5file)
    x.set_xpointer_path(path)
    self.grid.append(x)


def DanuSimulationLocation(sim_name='MAIN'):
  return '/Simulations/'+sim_name

def DanuSeriesLocation(id,sim_name='MAIN'):
  sim_loc=DanuSimulationLocation(sim_name)
  seq_grp_name=sim_loc+'/Series Data'
  seq_name='%s %d'%(seq_grp_name,id)
  return seq_name

'''
XdmfDanuDataset
'''
class XdmfDanuDataset(XdmfAttribute):


  def __init__(self,danu_seq,data_name):
    # Verify data is a field type
    try:
      data_attr=danu_seq.data_attributes(data_name)
    except:
      msg='Failed to read dataset %s attributes'%(data_name)
      raise Exception(msg)
    try:
      field_type=data_attr['FIELDTYPE']
    except KeyError:
      msg='Dataset %s is not a field type dataset'%(data_name)
      raise KeyError(msg)
    # Convert to the Xdmf center type
    if field_type == 'CELL':
      center='Cell'
    elif field_type == 'NODE':
       center='Node'
    else:
       msg='Invalid field type %s'%(field_type)
       raise ValueError(msg)
    # Find the dataset dimensions and set the type Scalar, Vector
    dims=danu_seq.get_data_dimensions(data_name)
    if len(dims) > 1:
      type='Vector'
    else:
      type='Scalar'
    # For now assume all data types are REAL8
    data_type='Float'
    precision=8
    # For Danu, file format is always HDF5
    format='HDF5'
    # Tricky to handle scalar and vector data field names
    # need the attribute name to appear correctly in the viewer (ParaView VisIt)
    import re
    field_names_re=re.compile('FIELDNAME')
    field_names=[n for n in data_attr if field_names_re.match(n)]
    if len(field_names) == 1 :
      aname=data_attr[field_names[0]]
    elif len(field_names) == 2:
      aname='[%s,%s]'%(data_attr[field_names[0]],data_attr[field_names[1]])
    elif len(field_names) == 3:
      aname='[%s,%s,%s]'%(data_attr[field_names[0]],
                         data_attr[field_names[1]],
                         data_attr[field_names[2]])
    else:
      msg='Dataset dims (%d) greater than 3 are not supported'%(len(field_names))
      raise Exception(msg)
    # Now ready to construct the XDMF Elements
    XdmfAttribute.__init__(self,aname,type,center)
    self.data_item=XdmfDataItem(data_name,data_type,precision,format)
    self.data_item.set_dimensions(dims)
    self.append(self.data_item)

  def set_data_location(self,loc):
    return self.data_item.set_text(loc)

  def get_data_location(self):
    return self.data_item.get_text()




def buildDanuDatastLocation(dname,seqname,simname='MAIN'):
  simulation_group_name='/Simulations/'+simname
  seq_group_name=simulation_group_name+'/Series Data'
  loc='%s/%s/%s'%(seq_group_name,seqname,dname)
  return loc





if __name__ == '__main__':

  root = XdmfRoot()
  time = XdmfTime(1.0)
  root.append(time)

  ElementTree.ElementTree(root.elem).write(sys.stdout)


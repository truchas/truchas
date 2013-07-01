"""
vtkWriter.py
Written by Sriram Swaminarayan

This module reads in a TBrook file, pulls out variables from it,
modifies them, and spits them back out as a VTK format file.

Please ensure that <truchasdir>/tools is in your PYTHONPATH variable.


"""

import os, sys

try:
    # Do path magic
    truchasdir = os.path.abspath(os.path.dirname(sys.argv[0]))
    if(len(truchasdir)):
        truchasdir = truchasdir+'/../'
    else:
        truchasdir = '../'
    sys.path.append(truchasdir)
    from PythonPackages import *
    sys.path.pop()
except:
    print """
    
    Unable to import standard modules Please ensure that
    <truchasdir>/tools is in your PYTHONPATH variable.

    It could also be that you don't have Numeric installed.

    """
    sys.exit(1)

from writeStorageObject import writeStorageObject
wso = writeStorageObject()

input = usubs.input('')


def vtkOpenFile(outfile, comment, format='ascii'):
    # Write a VTK file from scratch without a builtin writer
    fp = None
    fp = open(outfile,'w')
    fp.write('# vtk DataFile Version 1.0\n')
    fp.write(comment+'\n')

    if(format == 'binary'):
        fp.write('BINARY\n\n')
    else:
        fp.write('ASCII\n\n')
    return fp

def vtkWriteMesh(fp,f,m, format='ascii'):
    fp.write('DATASET UNSTRUCTURED_GRID\n')
    
    fp.write('POINTS %d double\n'%(m.nnodes))
    if(format == 'binary'):
        fp.write(m.nodeCoords.tostring()+'\n')
    else:
        wso.writeArray(2,m.nodeCoords,'','%12.5f ',fp,nperline=3)
    # VTK uses 0-(n-1) numbering
    cc = Numeric.array(m.cellConnectivity)
    cc = cc-1
    fp.write('CELLS %d %d\n'%(m.ncells,m.ncells*9))
    if(format == 'binary'):
        eight = Numeric.array([8],'i').tostring()
        for c in range(m.ncells):
            fp.write(eight+cc[c].tostring())
        fp.write('\n')
    else:
        wso.writeArray(2,cc,'8 ','%d ',fp,nperline=8)
    
    fp.write('CELL_TYPES %d\n'%m.ncells)
    
    tmpArray = Numeric.zeros((m.ncells,1),'i')
    for i in range(m.ncells):
        tmpArray[i] = 12
    if(format == 'binary'):
        fp.write(tmpArray.tostring()+'\n')
    else:
        wso.writeArray(2,tmpArray,'','%d ',fp,nperline=1)
    tmpArrray = None
    fp.write('CELL_DATA %d\n'%m.ncells)
    return 

def vtkAddVariables(fp, label, vars, format='ascii'):
    fp.write('FIELD %s %d\n'%(label,len(vars)))
    for var in vars:
        fp.write('%s 1 %d double\n'%(var.nickname,len(var.data)))
        if(format == 'binary'):
            fp.write(var.data.tostring()+'\n')
        else:
            wso.writeArray(1,var.data,'','%12.5f ',fp,nperline=1)


# load file and primary mesh
infile = input.usetc('filename', '/Users/sriram/codes/t/bsoft/xstf/stefan.TBrook.xml')


f = suFile(infile)  # try f.help() to get info on the suFile structure
m = suMesh(f.meshes[0],f)          # try m.help() to get info on the suMesh structure

comment = input.usetc('comment for vtk file?','Created from file '+infile)

fmt = 'binary'


ifile = 0
varnames = ["VOF0001","VOF0002","VOF0003", "VOF0004","VOF0005",
            "Enthalpy","T","Density"]
for ts in f.timesteps:
    outf = 'vtkfiles/b%010d.vtk'%ifile

    fp = vtkOpenFile(outf, comment, fmt)
    vtkWriteMesh(fp, f, m, fmt)
    vars = []
    # Create the Cp variable
    for v in varnames:
        vars.append( f.findVariable(ts,v))
        
    vtkAddVariables(fp, str(ts.cycle), vars, fmt)

    fp.close()
    print '%6d'%ifile, ts.cycle, 'done.'
    ifile = ifile + 1


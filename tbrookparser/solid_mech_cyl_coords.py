

"""
solid_mech_cyl_coords.py
Written by Dave Korzekwa and Sriram Swaminarayan

This module reads in a TBrook file, pulls out the displacement and
stress and strain tensors and transforms them to cylindrical
coordinates.  It is assumed that z is the cylinder axis.  The scalar
effective stress is also calculated and put in the GMV file.

The script is executed by entering:

python /path_to_script/solid_mech_cyl_coords TBrook_file_name

The user is prompted for a cycle number.  A binary GMV file is written
with displacements Dr, Dth, Dz and stress and strain components sigrr,
sigtt, sigzz, sigrt, sigrz, sigtz, etc.  The temperature, density and
volume fraction variables are also included in the GMV file to aid in
selecting regions to plot in gmv.

Please ensure that <truchasdir>/tools is in your PYTHONPATH variable.
This script uses the scriptUtils module written by Sriram Swaminarayan
to access the TBrook file data.


"""

import os, sys, getopt
try:
    # Do path magic
    truchasdir     = os.path.abspath(os.path.dirname(sys.argv[0]))
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


argv        = sys.argv[1:]
opts, rest  = getopt.getopt(argv, 'hf:e:')
if (not len(rest)):
     print
     print ' No input file specified on command line'
     print 
     sys.exit()
    
files = rest[0]

# load file and primary mesh
f = suFile(files)  # try f.help() to get info on the suFile structure
m = suMesh(f.meshes[0],f)          # try m.help() to get info on the suMesh structure

#print files

input = usubs.input()

tn = 0

while tn != -1:

    # Get cycle number
    tn = input.useti('\n Available cycles are: ' + str(f.cycles)+'\n'+'Choose a cycle number (-1 to quit): ',0)

    if tn == -1:
        sys.exit()
        
    for cycle in f.cycles:
        if tn == cycle:
            idx = f.cycles.index(cycle)
            ts = f.timesteps[idx]

    vvar = []
    for v in ts.vlist:
        if v.name.startswith('VOF'):
            vvar = vvar + [f.findVariable(ts,v.name)]

    # Get rank one cell based variable to create new tensor components
    TVar = f.findVariable(ts,'T')
    DenVar = f.findVariable(ts,'Density')
    #dumcVar = f.findVariable(ts,'T')

    # Get existing tensors
    sigVar = f.findVariable(ts,'sigma')
    epsVar = f.findVariable(ts,'epsilon')
    eplVar = f.findVariable(ts,'e_plastic')

    # Create tensor components to be transformed 
    sigrr = []
    sigtt = [] 
    sigzz = []
    sigrt = []
    sigrz = []
    sigtz = []

    epsrr = []
    epstt = [] 
    epszz = []
    epsrt = []
    epsrz = []
    epstz = []

    eplrr = []
    epltt = [] 
    eplzz = []
    eplrt = []
    eplrz = []
    epltz = []

    # Create effective stress variable
    sigeff = []



    # loop over cells
    for i in range(m.ncells):
        # Get cell centroids
        centx = m.cellCentroids[i][0]
        centy = m.cellCentroids[i][1]
        # Calculate rotation matrix and transpose
        radius =  sqrt(centx*centx + centy*centy)
        a11 = centx/radius
        a12 = centy/radius
        a13 = 0.0
        a21 = - centy/radius
        a22 = centx/radius
        a23 = 0.0
        a31 = 0.0
        a32 = 0.0
        a33 = 1.0
        
        arotmat = array([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]])
        arottrans = transpose(arotmat)
        # Stress components
        sxx = sigVar.data[i][0]
        syy = sigVar.data[i][1]
        szz = sigVar.data[i][2]
        sxy = sigVar.data[i][3]
        sxz = sigVar.data[i][4]
        syz = sigVar.data[i][5]
        sigeff.append(sqrt(((sxx - syy)**2 + (syy -szz)**2 + (szz - sxx)**2 +
                            6.0 * (sxy * sxy + sxz * sxz + syz * syz))/2.0))

        # Transform: [srt] = [A][sxy][A]^T
        sigmat = array([[sxx,sxy,sxz],[sxy,syy,syz],[sxz,syz,szz]])
        sigtmp = Numeric.dot(arotmat,sigmat)
        sigrot = Numeric.dot(sigtmp,arottrans)

        # Copy into components to write out 
        sigrr.append(sigrot[0,0])
        sigtt.append(sigrot[1,1])
        sigzz.append(sigrot[2,2])
        sigrz.append(sigrot[0,1])
        sigrt.append(sigrot[0,2])
        sigtz.append(sigrot[1,2])

        # Total strain components    
        exx = epsVar.data[i][0]
        eyy = epsVar.data[i][1]
        ezz = epsVar.data[i][2]
        exy = epsVar.data[i][3]
        exz = epsVar.data[i][4]
        eyz = epsVar.data[i][5]
        
        epsmat = array([[exx,exy,exz],[exy,eyy,eyz],[exz,eyz,ezz]])
        epstmp = Numeric.dot(arotmat,epsmat)
        epsrot = Numeric.dot(epstmp,arottrans)
        
        epsrr.append(epsrot[0,0])
        epstt.append(epsrot[1,1])
        epszz.append(epsrot[2,2])
        epsrz.append(epsrot[0,1])
        epsrt.append(epsrot[0,2])
        epstz.append(epsrot[1,2])

        # Plastic strain components
        epxx = eplVar.data[i][0]
        epyy = eplVar.data[i][1]
        epzz = eplVar.data[i][2]
        epxy = eplVar.data[i][3]
        epxz = eplVar.data[i][4]
        epyz = eplVar.data[i][5]

        eplmat = array([[epxx,epxy,epxz],[epxy,epyy,epyz],[epxz,epyz,epzz]])
        epltmp = Numeric.dot(arotmat,eplmat)
        eplrot = Numeric.dot(epltmp,arottrans)

        eplrr.append(eplrot[0,0])
        epltt.append(eplrot[1,1])
        eplzz.append(eplrot[2,2])
        eplrz.append(eplrot[0,1])
        eplrt.append(eplrot[0,2])
        epltz.append(eplrot[1,2])

    # Displacement components
    dxyVar = f.findVariable(ts,'Displacement')
    # Get rank one node based variable to create new tensor components
    dumnVar = f.findVariable(ts,'gapdisp')
    # Rotated displacement components
    Dr = []
    Dth = []
    Dz = []

    for i in range(m.nnodes):
        x = m.nodeCoords[i][0]
        y = m.nodeCoords[i][1]
        dx = dxyVar.data[i][0]
        dy = dxyVar.data[i][1]
        dz = dxyVar.data[i][2]
        # Radius of current node
        radius = sqrt(x * x + y * y)
        if  radius > 1.0e-10:
            # Projection of displacement on the radial unit vector
            dr = (dx * x + dy * y)/radius
            # Projection onto a unit vector perpendicular to the radial vector
            dth = - dx * y/radius + dy * x/radius
        else:
            # If the node is at x=y~=0 just use the radial component
            dr = sqrt(dx * dx + dy * dy)
            dth = 0.0

        Dr.append(dr)
        Dth.append(dth)
        Dz.append(dz)



    # Create tensor components to be transformed 
    srrVar = f.createVariable('sigrr',data=sigrr,map='CELL',type='d')
    sttVar = f.createVariable('sigtt',data=sigtt,map='CELL',type='d') 
    szzVar = f.createVariable('sigzz',data=sigzz,map='CELL',type='d')
    srtVar = f.createVariable('sigrt',data=sigrt,map='CELL',type='d')
    srzVar = f.createVariable('sigrz',data=sigrz,map='CELL',type='d')
    stzVar = f.createVariable('sigtz',data=sigtz,map='CELL',type='d')

    seffVar = f.createVariable('sigeff',data=sigeff,map='CELL',type='d')

    errVar = f.createVariable('epsrr',data=epsrr,map='CELL',type='d')
    ettVar = f.createVariable('epstt',data=epstt,map='CELL',type='d') 
    ezzVar = f.createVariable('epszz',data=epszz,map='CELL',type='d')
    ertVar = f.createVariable('epsrt',data=epsrt,map='CELL',type='d')
    erzVar = f.createVariable('epsrz',data=epsrz,map='CELL',type='d')
    etzVar = f.createVariable('epstz',data=epstz,map='CELL',type='d')

    eprrVar = f.createVariable('eplrr',data=eplrr,map='CELL',type='d')
    epttVar = f.createVariable('epltt',data=epltt,map='CELL',type='d') 
    epzzVar = f.createVariable('eplzz',data=eplzz,map='CELL',type='d')
    eprtVar = f.createVariable('eplrt',data=eplrt,map='CELL',type='d')
    eprzVar = f.createVariable('eplrz',data=eplrz,map='CELL',type='d')
    eptzVar = f.createVariable('epltz',data=epltz,map='CELL',type='d')

    # Rotated displacement components
    drVar  = f.createVariable('Dr',data=Dr,map='VERTEX',type='d')
    dthVar = f.createVariable('Dth',data=Dth,map='VERTEX',type='d')
    dzVar  = f.createVariable('Dz',data=Dz,map='VERTEX',type='d')





    # Write GMV file
    fileout   ='%s.%06d'%('cyl_coords.gmv',idx)    # The output file name
    outputfmt = 'binary'                       # The output format
    vars      = [TVar,DenVar,seffVar,srrVar,sttVar,srtVar,srzVar,stzVar,szzVar,
                 errVar,ettVar,ertVar,erzVar,etzVar,ezzVar,
                 eprrVar,epttVar,eprtVar,eprzVar,eptzVar,epzzVar,
                 drVar,dthVar,dzVar]   # variables to be put into GMV file
    # Add volume fraction data
    vars      = vars + vvar
    seq_no    = 1                             # The sequence number
    flags     = {}

    # Write our variables to a binary GMV file
    GMVwriteStorageObject(fileout, outputfmt, f.meshes[0],
                          ts.time, ts.cycle, vars, seq_no, flags)            


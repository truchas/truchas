docstr = \
"""
Further testing of the mapping infrastructure

Tests the mapping MA -> MB -> MA for predefined fields [C=1,C=x,C=x+y,C=step function] 
 MA = Mesh A (Exodus II format)
 MB = Mesh B (Exodus II format)

Tests both 'constant preserving' and 'conservative' mappings.

Please consult the README_extex file for further information about
testing the mapping.

Test parameters are imported from the file given with the -f command
line option. That file is a python-syntax assignment of various
parameters.  Provide only the prefix, the suffix must be '.py'.
The parameters needed are:

  InputsDir    : the source directory for input files
  mesha_files  : a dict of 'A' input mesh files and types
  meshb_files  : a dict of 'B' input mesh files and types
  filename     : an output filename
  graphicsname : a graphics filename prefix

The parameters are accessed thusly: 'p.<param>'
"""

#------------------------------------------------------------------------
#prelims needed to perform these tests
import os, sys, string

try:
   import numpy.oldnumeric as Numeric
except ImportError:
   import Numeric
except:
   raise

if ('TRUCHAS' in os.environ.keys()):
    "pick up TRUCHAS home directory from environment by default"
    truchasdir = os.environ['TRUCHAS']
else:
    "otherwise, assume a relative location to this file. (bad form!!!)"
    __cwd      = os.getcwd()
    truchasdir = os.path.abspath(__cwd + '../../../../../../')

sys.path = sys.path + \
           [truchasdir +'/tools/PythonPackages/TBrookParser/'] + \
           [truchasdir +'/tools/PythonPackages/TBrookParser/MODutils'] 

from PYTHONutils           import uTestOpts
from MAPVariable           import mapVariable
from MODMesh               import modMesh
from GMVwriteStorageObject import GMVwriteStorageObject

# Parse command-line for input parameters file.
o = uTestOpts('fdH',
              defaults={'f' : 'param_Tet-Hex-Tet', 'H': False},
              actions ={'H' : 'store_true'},
              dir = truchasdir+\
                    'tools/PythonPackages/TBrookParser/MODutils/MAPutils')
(opt,a) = o.parse_args()

if (opt.H):
   print docstr
   sys.exit(0)
debug = opt.d
p = __import__(opt.f) # a bunch of parameters

#------------------------------------------------------------------------
# handy classes used below
class var:

    #generic variable class 

    def __init__(self):

        self.name      = 'thisvar'
        self.nickname  = self.name
        self.rank      = 1
        self.shape     = [1]
        self.type      = 'd'
        self.data      = None
        self.mesh      = 'DefaultMesh'
        self.meshspace = 'CELL'

class field:

    #generic field class 

    def __init__(self):

        self.name  = None
        self.rank  = 1
        self.shape = None
        self.type  = None
        self.data  = None

    def getData(self,name,shape,data,type='d'):

        self.name  = name
        self.shape = shape
        self.type  = type
        self.data  = data
        

#now define a writer class to present diagnostics neatly

class write:

    #class to write out diagnostics to a file

    def header(self,fp,filea,fileb,maprule,fieldtitle):

        S  = '-' * 115
        print >> fp, S
        print >> fp, 'Field: %+6s  Map: %+10s  Transfer:%+10s ' \
                     %(fieldtitle,maprule,'MA->MB->MA ' )
        print >> fp
        print >> fp, 'MA = %s         MB = %s ' %(filea,fileb)
        print >> fp, S


    def diagnostics(self,fp,vars):

        fmt  = '%15.5E'
        size = len(vars)

        for i in range(size):

            int = vars[i].newintgrl
            
            if i==0:
                #print header                
                D    = ['Ave','Stdev','Max','Min','Intgrl','L0 error', 'L2 error']
                S    = ''
                for stat in D:
                    S = S + '%+15s ' %(stat)
                print >> fp, S
                #set first field integral
                int = vars[1].oldintgrl
            
            datasize = len(vars[i].data)            
            identity = Numeric.array([1],'d')
            identity = Numeric.resize(identity,[datasize])
            ave      = float(Numeric.sum(vars[i].data))/datasize
            tmp      = (vars[i].data-ave*identity)**2
            stdev    = Numeric.sum(tmp)
            stdev    = (stdev/datasize)**(0.5)
            max      = Numeric.argmax(vars[i].data)
            max      = vars[i].data[max]
            min      = Numeric.argmin(vars[i].data)
            min      = vars[i].data[min]

            L0       = -1000.0
            L2       = -1000.0

            if i > 1:
               #calculate L0,L2 errors
               assert len(vars[0].data) == len(vars[i].data)
               tmp  = (vars[i].data - vars[0].data)**2
               tmp  = tmp**(0.5)
               arg  = Numeric.argmax(tmp)
               L0   = tmp[arg]
               tmp  = (vars[i].data - vars[0].data)**2
               L2   = Numeric.sum(tmp)/datasize
               L2   = L2**(0.5)

            D2   = [ave, stdev, max, min, int, L0, L2]

            thismesh = 'MA'
            if i==1 or i==3: thismesh = 'MB'

            S = '%+3s ' %(thismesh)
            for stat in D2:
                S = S + fmt%(stat) 

            print >> fp,S
            print >> fp
            print >> fp

#------------------------------------------------------------------------
# Begin the tests

mesh_msnames = ['CELL', 'VERTEX', 'FACE']
wdir         = os.getcwd()
meshdir      = wdir + '/../EXOutils'

#create file handler to output results to
fp       = open(p.filename, 'w')

#initialise field classes
field1 = field()
field2 = field()
field3 = field()
field4 = field()

# Loop over 'A' input files.
combos  = 0
for mFileA in p.mesha_files.keys():
    
    cmd   = 'cp %s%s .' %(p.InputsDir,mFileA)
    os.system(cmd)

    clean = 'clean'
    if combos < 1:
       clean = 'cleanest'

    try:
       mesha = modMesh(p.mesha_files[mFileA], # type (from dict)
                       mesh_msnames,          
                       mFileA,                # name
                       meshdir,
                       1.0,
                       clean,
                       debug=1)
    except:
       print
       print 'Problems occurred in generating the Exodus mesh from %s' \
             %(mFileA)
       print 'Exiting this component test'
       print
       raise

    for meshspace in mesha.mslist:
        if meshspace.name == 'CELL':
            mesha.fillCells(meshspace)
            mesha_ncells = meshspace.size
        if meshspace.name == 'VERTEX':
            mesha.fillVertices(meshspace)

    a_coords    = mesha.vertices['COORDS']
    ndim        = len(a_coords[0])
    a_connect   = mesha.cells['VERTICES']
    a_centroids = Numeric.array([0],'d')
    a_centroids = Numeric.resize(a_centroids,[mesha_ncells,ndim])

    for k in range(ndim):
        for i in range(len(a_connect)):
            sum   = 0.0
            cnt = 0
            for j in a_connect[i]:
                sum = sum + a_coords[j-1,k]
                cnt = cnt + 1
            a_centroids[i,k] = sum/cnt


    #now define some fields [C=1,C=x,C=y,C=x+y]

    #define C=1 field:
    fdata   = Numeric.array([1],'d')
    fdata   = Numeric.resize(fdata,[mesha_ncells])
    fname   = 'C=1'
    field1.getData(name=fname,shape=[mesha_ncells],data=fdata)

    #define C=x field:
    fdata   = a_centroids[:,0]+10
    fname   = 'C=x+10'
    field2.getData(name=fname,shape=[mesha_ncells],data=fdata)

    #define C=x+y field:
    fdata   = a_centroids[:,0] + a_centroids[:,1]+10
    fname   = 'C=x+y+10'
    field3.getData(name=fname,shape=[mesha_ncells],data=fdata)

    #define C=step function:
    x      = a_centroids[:,0]
    ave    = Numeric.sum(x)/len(x)
    fdata  = Numeric.greater(x,ave)
    fname  = 'C=Step'
    field4.getData(name=fname,shape=[mesha_ncells],data=fdata)

    # Loop over 'B' files, mapping to and from...
    for mFileB in p.meshb_files.keys():

        cmd = 'cp %s%s .' %(p.InputsDir,mFileB)
        os.system(cmd)
        
        try:
            meshb = modMesh(p.meshb_files[mFileB], # type (from dict)
                            mesh_msnames,
                            mFileB,                # name
                            meshdir,
                            1.0,
                            clean,
                            debug=1)
        except:
            print
            print 'Problems occurred in generating the Exodus mesh from %s' \
                  %(mFileB)
            print 'Exiting this component test'
            print
            raise

        vars_A1_CP  = []
        vars_A1_CO  = []
        vars_B_CP   = []
        vars_B_CO   = []
        vars_A3_CP  = []
        vars_A3_CO  = []

        # Rules loop
        for maprule in ['constants_preserving','conservative']:

            for field in [field1,field2,field3,field4]:

                vars             = []
                
                newvar           = var()
                newvar.name      = field.name
                newvar.nickname  = newvar.name
                newvar.rank      = field.rank
                newvar.shape     = field.shape
                newvar.type      = field.type
                newvar.meshspace = 'CELL'
                newvar.mesh      = 'ExoMesh'
                newvar.oldintgrl = 0.0
                newvar.newintgrl = 0.0

                write().header(fp,mFileA,mFileB,maprule,field.name)

                #now change the variable shape to conform to mesh 'a'
                newvar.shape[0]  = mesha_ncells
                newvar.data      = Numeric.resize(field.data,newvar.shape)

                vars.append(newvar)

                if maprule == 'conservative':
                    vars_A1_CO.append(newvar)
                elif maprule == 'constants_preserving':
                    vars_A1_CP.append(newvar)

                #now map the field variable from MA->MB->MA
        
                #first MA->MB
                newvar1   = newvar
                newvar1   = mapVariable(newvar,mesha,meshb,wdir,mustmap=1,
                                        maprule=maprule,clean=clean)
                newvar1.getData()
                vars.append(newvar1)
                
                if maprule == 'conservative':
                    vars_B_CO.append(newvar1)
                elif maprule == 'constants_preserving':
                    vars_B_CP.append(newvar1)

                #now MB->MA
                newvar2 = newvar1
                newvar2 = mapVariable(newvar1,meshb,mesha,wdir,mustmap=1,
                                      maprule=maprule,clean=clean)
                newvar2.getData()
                vars.append(newvar2)

                if maprule == 'conservative':
                    vars_A3_CO.append(newvar2)
                elif maprule == 'constants_preserving':
                    vars_A3_CP.append(newvar2)                

                #provide diagnostics
                write().diagnostics(fp,vars)

                combos = combos + 1 # N(A)*N(B)*N(rules)

        # End of Rules loop
        
        #now write out GMV files for each mesh/mesh/rule combination
        if mFileA != mFileB:
            outfmt    = 'ascii'
            extension = '.gmv'
            D         = {'A1_CP' : vars_A1_CP, 'A1_CO' : vars_A1_CO,
                         'B_CP'  : vars_B_CP , 'B_CO'  : vars_B_CO ,
                         'A3_CP' : vars_A3_CP, 'A3_CO' : vars_A3_CO}
            D2        = {'A1' : mesha,
                         'B'  : meshb,
                         'A3' : mesha}

            #strip mFileA and mFileB of their extensions
            fila = mFileA[:-4]
            filb = mFileB[:-4]
            prefix = '%s%s-%s_' %(p.graphicsname,fila,filb)
            
            for mesh in ['A1','B','A3']:
                for maptype in ['CP','CO']:
                    specifier = '%s_%s' %(mesh,maptype)
                    outfile   = prefix + specifier + extension
                    GMVwriteStorageObject(outfile, outfmt, D2[mesh],
                                          0.0, 0, D[specifier])        

    # End of 'B' loop
# End of 'A' loop

#close the file
fp.close()

#now clean up stray files
cmd = 'rm *.exo'
os.system(cmd)
cmd = 'rm diagnostics_file int_vols_file'
os.system(cmd)


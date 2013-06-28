docstr = \
"""
 Usage:
 > python test_t_to_ex_mapping.py [-Hd] [-f <params>]

    -H          : print this message and exit
    -d          : set debug flag
    -f <params> : import test parameters from <params>.py


  Testing the mapping of Truchas to Exodus II non-degenerate hex meshes.

  Tests the mapping of fields from MA -> MB
    MA = Mesh A (comes from a Truchas TBrook.xml file)
    MB = Mesh B (non degenerate HEX mesh in Exodus II format)

  Please consult the README_t2ex file for further information
  about testing this mapping.

""" 
#------------------------------------------------------------------------
#prelims needed to perform these tests
import os, sys, string

thisdir      = os.path.dirname(sys.argv[0])
truchasdir   = os.path.abspath(thisdir) +'../../../../../../'
truchasdir   = os.path.abspath(truchasdir)

sys.path     = sys.path + \
               [truchasdir +'/tools/PythonPackages/TBrookParser/'] + \
               [truchasdir +'/tools/PythonPackages/TBrookParser/MODutils']
    
from PYTHONutils import uTestOpts
"Parse command-line for input parameters file and options."
o = uTestOpts('fdH',
              defaults={'f' : 'param_t2ex', 'H': False},
              actions ={'H' : 'store_true',},
              dir = truchasdir+\
              'tools/PythonPackages/TBrookParser/MODutils/MAPutils')
(opt,a) = o.parse_args()
if (opt.H):
   print docstr
   sys.exit(0)

"Import the parameters file"
try:
   params = __import__(opt.f) # a bunch of parameters    
except:
   print "\nFailed to import paramaters from file '%s.py'\n" %(opt.f)
   raise

# End prelims
#------------------------------------------------------------------------
# Internals

class write:

    #class to write out diagnostics to a file

    def header(self,fp,filea,fileb,vara,varb,map=1):

        if not map:
            varb.maprule = 'None'
            fileb        = filea
            
        S  = '-' * 115
        print >> fp, S
        print >> fp, 'Variable: %+6s  Map: %+10s  Transfer:%+10s ' %(vara.name,varb.maprule,'MA->MB ')
        print >> fp
        print >> fp, 'MA = %s.TBrook.xml         MB = %s ' %(filea,fileb)
        print >> fp, S


    def diagnostics(self,fp,mesha,meshb,vara,varb,map=1):

        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise
       
        fmt   = '%15.5E'
        size  = 2
        count = 0

        for var in [vara,varb]:

            inta  = 0.0
            intb  = 0.0
            L0    = -1000.0
            L2    = -1000.0
            ave   = 0.0
            stdev = 0.0
            max   = 0.0
            min   = 0.0

            if map:
                inta     = varb.oldintgrl
                intb     = varb.newintgrl
                datasize = len(var.data)
                identity = Numeric.array([1],'d')
                identity = Numeric.resize(identity,[datasize])
                ave      = float(Numeric.sum(var.data))/datasize
                tmp      = (var.data-ave*identity)**2
                stdev    = Numeric.sum(tmp)
                stdev    = (stdev/datasize)**(0.5)
                max      = Numeric.argmax(var.data)
                max      = var.data[max]
                min      = Numeric.argmin(var.data)
                min      = var.data[min]


            if count == 0:
                thismesh = 'MA'
                D2       = [ave, stdev, max, min, inta, L0, L2]
            if count == 1:
                thismesh = 'MB'
                if map:
                    #calculate L0,L2 errors when meshes are identical
                    if mesha.vertices['N'] == meshb.vertices['N']:
                        if len(vara.data) == len(varb.data):
                            tmp  = (vara.data - varb.data)**2
                            tmp  = tmp**(0.5)
                            arg  = Numeric.argmax(tmp)
                            L0   = tmp[arg]
                            tmp  = (vara.data - varb.data)**2
                            L2   = Numeric.sum(tmp)/len(vara.data)
                            L2   = L2**(0.5)
                    D2 = [ave, stdev, max, min, intb, L0, L2]

            if map:
                if count==0:
                   #print header                
                   D    = ['Ave','Stdev','Max','Min','Intgrl','L0 error', 'L2 error']
                   S    = ''
                   for stat in D:
                       S = S + '%+15s ' %(stat)
                   print >> fp, S

                S = '%+3s ' %(thismesh)
                for stat in D2:
                    S = S + fmt%(stat) 

                print >> fp,S
                print >> fp
                print >> fp

            count = count + 1
#------------------------------------------------------------------------

# Begin the tests...

"First get the TBrook storage object from the Truchas simulatiom"
try:
   tsim      = params.TBprfx
   # Get storage object and mesh information
   from XMLgetStorageObject       import getStorageObject
   T_Sim = params.InputsDir + params.T_Sim
   assert os.path.exists(T_Sim), \
          'File %s does not exist. Check the path and file prefix.' %(T_Sim)

   storage   = getStorageObject(T_Sim)
   mesha    = storage.mlist[0] #assume default mesh from Truchas XML file
   for meshspace in mesha.mslist:
      meshspace.vlist = storage.getValues(meshspace.vlist)
   mesha.fillMesh()

   # Set cycle numer and retrieve time and variable list
   cycle = params.CYCLE
   if (cycle != -1):
      assert cycle in storage.timesteps.cyclelist, \
             ' \n \n Cycle %i is not part of %s \n' \
             %(cycle, str(storage.timesteps.cyclelist))
   else:
      cycle = storage.timesteps.cyclelist[-1]
   itimeout = storage.timesteps.cyclelist.index(cycle)
   time     = storage.tlist[itimeout].time
   varsa    = storage.tlist[itimeout].vlist
   storage.getValues(varsa)

   # Set coordinate scale factor
   if (params.CSF):
      csf = params.CSF
   else:
      csf = storage.specs.csf 
except:
   print 'Problem with TBrook file: %s' %(T_Sim)
   print 'Exiting this component test'
   raise

"EXODUS file/type"

exofile  = params.EXOfile.keys()[0]
cmd          = 'cp %s%s .' %(params.InputsDir,exofile)
os.system(cmd)
   
mesh_type    = params.EXOfile[exofile]
if mesh_type!='HEX':
   raise 'EXODUS meshtype error: expected HEX, got %s' %(ExoMTyp)

msnames = ['CELL', 'VERTEX', 'FACE']
wdir         = os.getcwd()
meshdir      = wdir + '/../EXOutils'

from MODMesh                   import modMesh
try:
   meshb = modMesh(mesh_type,msnames,
                   exofile,meshdir,
                   scale_factor=csf)
except:
   print
   print 'Problems occurred in generating the Exodus II mesh from %s' \
         %(exofile)
   print 'Exiting this component test'
   print
   raise

"Now map the field variables that live on Mesh A to Mesh B"

from MODFieldVariable import modFieldVariable
varsb  = []
for var in varsa:
    if var.rank > 0:
        mapvardir   = os.getcwd()
        newvar      = modFieldVariable(var, mesha, meshb,
                                       mapvardir, debug=opt.d)
        newvar.mesh = meshb.name
        varsb.append(newvar)

"""
Output:
"""
try:
    os.chdir(params.OutputDir)
    os.chdir('../')
except:
    os.mkdir(params.OutputDir)
    
diagfile = '%s/Diagnostics.txt' %(params.OutputDir)
fp       = open(diagfile, 'w')

# diagnostics of the mapping
# - only for variables that will be mapped in the restart 
mvarsa = []
mvarsb = []
for vara in varsa:
    mappedvar = None
    for varb in varsb:
        if varb.mustmap and vara.name == varb.name:
            mvarsa.append(vara)
            mvarsb.append(varb)

for vara in varsa:
    for varb in varsb:
        if vara in mvarsa and varb in mvarsb :
            if vara.name == varb.name:
                write().header(fp,tsim,exofile,vara,varb)
                write().diagnostics(fp,mesha,meshb,vara,varb)
for vara in varsa:
    for varb in varsb:
        if vara not in mvarsa and varb not in mvarsb:
            if vara.name == varb.name and vara.rank > 0:
                write().header(fp,tsim,exofile,vara,varb,map=0)
                write().diagnostics(fp,mesha,meshb,vara,vara,map=0)

#close the diagnostics file
fp.close()

"""
Write out the binary RESTART file resulting from the mapping
"""
from RESTARTwriteStorageObject import RESTARTwriteStorageObject

graphicsname = '%s/%s_to_%s' %(params.OutputDir,tsim,exofile[:-4])
extension    = '.restart.%i' %(cycle)
outfile      = graphicsname + extension
outfmt       = 'binary'
RESTARTwriteStorageObject(outfile,outfmt,storage,meshb,itimeout,time,varsb)

"""
Write out ascii GMV files for the original XML simulation
and the new mapped simulation
"""
from GMVwriteStorageObject     import GMVwriteStorageObject

#remove Solid Mechanics RHS variable from variable lists
for vara in varsa:
    if 'RHS' in vara.name:
        del varsa[varsa.index(vara)]

for varb in varsb:
    if 'RHS' in varb.name:
        del varsb[varsb.index(varb)]

graphicsname = '%s/%s' %(params.OutputDir,tsim)
outfmt       = 'ascii'
extension    = '.gmv.%i' %(cycle)
outfile      = graphicsname + extension
GMVwriteStorageObject(outfile,outfmt,mesha,time,itimeout,varsa)        
        
graphicsname = '%s/%s_to_%s' %(params.OutputDir,tsim,exofile[:-4])
outfmt       = 'ascii'
extension    = '.gmv.%i' %(cycle)
outfile      = graphicsname + extension
GMVwriteStorageObject(outfile,outfmt,meshb,time,itimeout,varsb)        

#------------------------------------------------------------------------
"All Done "
print "Done: Output written to directory %s" %(params.OutputDir)

"Clean up files not needed"
cmd = 'rm *.exo'
os.system(cmd)
cmd = 'rm diagnostics_file int_vols_file'
os.system(cmd)


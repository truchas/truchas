docstr = \
"""
 Usage:
 > python test_2t_to_ex_mapping.py [-Hd] [-f <params>]

    -H          : print this message and exit
    -d          : set debug flag
    -f <params> : import test parameters from <params>.py

 Testing the mapping of 2 Truchas meshes to an Exodus II
 non-degenerate hex meshes.

 Tests the mapping of fields from MA + MA_2 -> MB
   MA   = 1st source mesh  (from a Truchas TBrook.xml file)
   MA_2 = 2nd source mesh  (from a 2nd Truchas TBrook.xml file)
   MB   = destination mesh (from a non-degenerate HEX mesh
                            in Exodus II format)

 Please consult the README_2t2ex file for further information
 about testing this mapping.

""" 
#------------------------------------------------------------------------
# prelims
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
              defaults={'f' : 'param_2t2ex', 'H': False},
              actions ={'H' : 'store_true',},
              dir = truchasdir+\
              'tools/PythonPackages/TBrookParser/MODutils/MAPutils')
(opt,a) = o.parse_args()
if (opt.H):
   print docstr
   sys.exit()

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

    # a writer class to present diagnostics neatly

    def header(self,fp,filea,fileb,vara,varb,map=1):

        if not map:
            varb.maprule = 'None'
            fileb        = filea
            
        S  = '-' * 115
        print >> fp, S
        print >> fp, 'Variable: %+6s  Map: %+10s  Transfer:%+10s ' %(vara.name,varb.maprule,'MA->MB ')
        print >> fp
        print >> fp, 'MA = %s         MB = %s ' %(filea,fileb)
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
                   D    = ['Ave','Stdev','Max','Min','Intgrl','L0 error','L2 error']
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

def getTSim(storage,cycle):
    
    # get a source mesh and all field values from
    # the Truchas storage object

    mesha    = storage.mlist[0] #assume default mesh from Truchas XML file

    # now get the data values for mesha
    for meshspace in mesha.mslist:
        meshspace.vlist = storage.getValues(meshspace.vlist)

    # having obtained the data values
    # now fill in the following convenient dictionaries
    mesha.fillMesh()

    #now get the field data associated with this timestep

    cyclelist = storage.timesteps.cyclelist
    
    for i in cyclelist:
        if cycle == i: itimeout = cyclelist.index(i)

    time         = storage.tlist[itimeout].time
    varsa        = storage.tlist[itimeout].vlist
    storage.getValues(varsa)

    return mesha, time, varsa

# End Internals
#------------------------------------------------------------------------

# Begin the tests...

from XMLgetStorageObject       import getStorageObject

T_Sim     = params.InputsDir + params.T_Sim

"First TBrook file..."
try:
   T_Sim     = params.InputsDir + params.T_Sim
   storage   = getStorageObject(T_Sim)
   csf1      = storage.specs.csf
   cycle1 = params.CYCLE[0]
   if cycle1 != -1:
      assert cycle1 in storage.timesteps.cyclelist, \
             ' \n \n For T simulation : %s  \n cycle %i is not in %s \n'\
             %(T_Sim,cycle1,str(storage.timesteps.cyclelist))
   else:
      cycle1 = storage.timesteps.cyclelist[-1]

   mesha,timea,varsa  = getTSim(storage,cycle1)
   tsim1 = '%s.TBrook.xml'%(params.TBprfx[0])
except:
   print 'Problem with first TBrook file: %s' %(params.T_Sim)
   print 'Exiting this component test'
   raise

"Second TBrook file..."
try:
   T_Sim2    = params.InputsDir + params.T_Sim2
   storage2  = getStorageObject(T_Sim2)
   csf2      = storage2.specs.csf
   cycle2 = params.CYCLE[1]
   if cycle2 != -1:
      assert cycle2 in storage2.timesteps.cyclelist, \
             ' \n \n For T simulation : %s \n Cycle %i is not in %s \n' \
             %(T_Sim2,cycle2,str(storage2.timesteps.cyclelist))
   else:
      cycle2 = storage2.timesteps.cyclelist[-1]
   mesha2,timea2,varsa2 = getTSim(storage2,cycle2) 
   tsim2 = '%s.TBrook.xml'%(params.TBprfx[1])
except:
   print 'Problem with second TBrook file: %s' %(params.T_Sim2)
   print 'Exiting this component test'
   raise

tsim  = '%s_and_%s' %(params.TBprfx[0],params.TBprfx[1]) 

"Scale factors"
if csf1 != csf2 :
   print
   print 'Note scaling factors are different (%s,%s)' %(csf1,csf2)
   print 'Choose an appropriate csf variable in the params file.'
   print
if params.CSF: csf = params.CSF
else:          csf = csf1

#now get Mesh B and its data values

"EXODUS file/type" 
exofile  = params.EXOfile.keys()[0]
cmd  = 'cp %s%s .' %(params.InputsDir,exofile)
os.system(cmd)

mesh_type = params.EXOfile[exofile]
if mesh_type != 'HEX':
   raise 'EXODUS meshtype error: expected HEX, got %s' %(ExoMTyp)
   
mesh_msnames = ['CELL', 'VERTEX', 'FACE']
wdir         = os.getcwd()
meshdir      = wdir + '/../EXOutils'

from MODMesh                   import modMesh
try:
    meshb = modMesh(mesh_type,mesh_msnames,
                    exofile,meshdir,
                    scale_factor=csf)
except:
    print
    print 'Problems occurred in generating the Exodus II mesh from %s' \
          %(exofile)
    print 'Exiting this component test'
    print
    raise

"""
Now map the field variables that live on Mesh A + Mesh A2 -> Mesh B
"""

"""
First perform the mapping Mesh A -> Mesh B
 (external regions are filled with default values)
"""
from MODFieldVariable import modFieldVariable
varsb1 = []
for var in varsa:
    if var.rank > 0:
        mapvardir   = os.getcwd()
        newvar      = modFieldVariable(var, mesha, meshb,
                                       mapvardir, debug=opt.d)
        newvar.mesh = meshb.name
        varsb1.append(newvar)
"""
Now perform the mapping Mesh A2 -> Mesh B
(external regions are filled with default values)
"""
varsb2 = []
for var in varsa2:
    if var.rank > 0:
        mapvardir    = os.getcwd()
        newvar2      = modFieldVariable(var, mesha2, meshb,
                                        mapvardir, debug=opt.d)
        newvar2.mesh = meshb.name
        varsb2.append(newvar2)
"""
Finally, add the two variables' data together to create a new mapped variable
"""
from CREATEVariable import createVariable
varsb = []
for var in varsb1:
    if var.rank > 0:
        shape    = []
        shape[:] = var.shape[:]
        for var2 in varsb2:
            if var.name == var2.name:
                data     = var.data + var2.data
                newvar   = createVariable(var.name,var.nickname,var.rank,
                                          shape,var.type,data,var.mesh,
                                          var.meshspace,var.mustmap,
                                          var.maprule)
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
    for varb in varsb1:
        if varb.mustmap and vara.name == varb.name:
            mvarsa.append(vara)
            mvarsb.append(varb)
for vara in varsa:
    for varb in varsb1:
        if vara in mvarsa and varb in mvarsb :
            if vara.name == varb.name:
                write().header(fp,tsim1,exofile,vara,varb)
                write().diagnostics(fp,mesha,meshb,vara,varb)
for vara in varsa:
    for varb in varsb1:
        if vara not in mvarsa and varb not in mvarsb:
            if vara.name == varb.name and vara.rank > 0:
                write().header(fp,tsim1,exofile,vara,varb,map=0)
                write().diagnostics(fp,mesha,meshb,vara,vara,map=0)
mvarsa2 = []
mvarsb2 = []
for vara in varsa2:
    for varb in varsb2:
        if varb.mustmap and vara.name == varb.name:
            mvarsa2.append(vara)
            mvarsb2.append(varb)

for vara in varsa2:
    for varb in varsb2:
        if vara in mvarsa2 and varb in mvarsb2 :
            if vara.name == varb.name:
                write().header(fp,tsim2,exofile,vara,varb)
                write().diagnostics(fp,mesha2,meshb,vara,varb)
for vara in varsa2:
    for varb in varsb2:
        if vara not in mvarsa2 and varb not in mvarsb2:
            if vara.name == varb.name and vara.rank > 0:
                write().header(fp,tsim2,exofile,vara,varb,map=0)
                write().diagnostics(fp,mesha2,meshb,vara,vara,map=0)

#close the diagnostics file
fp.close()


"""
Write out the binary RESTART file resulting from the mapping
"""
from RESTARTwriteStorageObject import RESTARTwriteStorageObject
itimeout     = 0
time         = 0.0
cycle        = 0
tmp          = tsim + '_to_' + exofile[0:-4]
graphicsname = '%s/%s' %(params.OutputDir,tmp)
extension    = '.restart.%i' %(cycle)
outfile      = graphicsname + extension
outfmt       = 'binary'
RESTARTwriteStorageObject(outfile,outfmt,storage,meshb,itimeout,time,varsb)

"""
Write out ascii GMV files for the original XML simulations
and the new mapped simulation
"""
from GMVwriteStorageObject     import GMVwriteStorageObject

graphicsname = '%s/%s' %(params.OutputDir,params.TBprfx[0])
outfmt       = 'ascii'
extension    = '.gmv.%i' %(cycle1)
outfile      = graphicsname + extension
GMVwriteStorageObject(outfile,outfmt,mesha,time,itimeout,varsa)

graphicsname = '%s/%s' %(params.OutputDir,params.TBprfx[1])
outfmt       = 'ascii'
extension    = '.gmv.%i' %(cycle2)
outfile      = graphicsname + extension
GMVwriteStorageObject(outfile,outfmt,mesha2,time,itimeout,varsa2)  
        
graphicsname = '%s/%s' %(params.OutputDir,tsim)
outfmt       = 'ascii'
extension    = '.gmv.%i' %(cycle)
outfile      = graphicsname + extension
GMVwriteStorageObject(outfile,outfmt,meshb,time,itimeout,varsb)        

#------------------------------------------------------------------------
"All Done "
print "Done: Output written to directory %s" %(params.OutputDir)

# "Clean up files not needed"
# cmd = 'rm *.exo'
# os.system(cmd)
# cmd = 'rm diagnostics_file int_vols_file'
# os.system(cmd)
#------------------------------------------------------------------------

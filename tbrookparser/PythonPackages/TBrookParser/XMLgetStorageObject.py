"""

 XMLgetStorageObject

 -----------------------------------------------------------------------------
  Purpose:

     Instantiates a storage object and initializes some of its
     internal Metadata.

     A special method (getValues) is provided to be used
     when the actual data is needed (more than the Metadata).

  Public Interface(s):

   SO = getStorageObject(filename)
   SO.getValues(vars)
   SO.str()
   
  Contains:
    class getStorageObject
        __init__(self,filename,format='brookXML',parse=MiniDom,debug=0)
        __nameMap(self,oldpattern,newpattern)
        getValues(self,vars)
        str(self)

    Unit Test Block

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""
import sys
if __name__=='__main__':
    " Set sys.path for component testing mode "
    import os
    thisdir   = os.path.dirname(__file__)
    parsedir  = os.path.abspath(thisdir+'../') 

    sys.path.append(parsedir)

from PYTHONutils import getpath
from XMLutils    import MiniDom,getSimSpecs,getMesh,getTimeSteps,getProbes,getAborts

class getStorageObject:

    def __init__(self,filename,fp=sys.stdout,format='brookXML',parse=MiniDom,debug=0):

        assert format == 'brookXML'
        
        "initialise the storage object"

        self.specs          = None #the simulation specifications
        self.mlist          = []   #a set of meshes
        self.tlist          = []   #a set of timestep structures
        self.aborts         = []   #a set of abort structures
        self.timesteps      = None #a timesteps object (id, cycle, times lists)
        self.plist          = []   #a set of probes
        self.varlist        = []   #a set of user specified variables
        self.namelu         = {}   #maps nicknames to names
        self.spacelu        = {}   #maps nickname spacenames to actual spacenames
        self.fp             = fp   #file pointer to direct screen output to
        
        "establish space nickname lookup"

        self.spacelu = { None      :'CELL',
                         'cell'    :'CELL',
                         'cells'   :'CELL',
                         'NONE'    :'NONE',
                         'CELL'    :'CELL',
                         'CELLS'   :'CELL',
                         'FACES'   :'FACE',
                         'FACE'    :'FACE',
                         'VERTICES':'VERTEX',
                         'VERTEX'  :'VERTEX',
                         'NODES'   :'VERTEX',
                         'NODE'    :'VERTEX',
                         'node'    :'VERTEX',
                         'nodes'   :'VERTEX',
                         'EDGES'   :'EDGE',
                         'EDGE'    :'EDGE',
                         'MESH'    :'MESH'}

        """
        establish variable nickname lookup
        """

        self.namelu = {'Z_RHO'              :'Density',
                       'Z_RHOOLD'           :'Density_Old',
                       'Z_TEMP'             :'T',
                       'Z_TEMP_OLD'         :'T_Old',
                       'Z_ENTHALPY'         :'Enthalpy',
                       'Z_ENTHALPYOLD'      :'Enthalpy_Old',
                       'Z_P'                :'P',
                       'Z_VC'               :'Velocity',
                       'Z_VCOLD'            :'Velocity_Old',
                       'Z_VF'               :'V_Face',
                       'VC'                 :'Velocity',
                       'VCOLD'              :'Velocity_Old',
                       'VF'                 :'V_Face',
                       'RHS'                :'SM_RHS',
                       'Plastic_Strain_Rate':'epsdot',
                       'Thermal_Strain'     :'epstherm',
                       'PC_Strain'          :'epspc',
                       'Displacement'       :'disp',
                       'DISPLACEMENT'       :'DISP',
                       'Gradient_T'         :'GradT',
                       'Grad(T)'            :'Grad_T_',
                       'del-rho'            :'delrho'}

        self.debug          = debug

        "get the simulation specifications"

	if self.debug: print '  Debug > getStorageObject: getSimSpecs'
        self.specs          = getSimSpecs(filename,parse,self.fp)
        
        "get all meshes"

        elems               = parse().getElements(self.specs.thedoc,'MESH')

        for element in elems:
            
            filenode        = parse().getElement(element,'FILE')
            filefmt         = parse().getElemAttribute(filenode,'FORMAT')
            assert filefmt  == 'xml'
            meshfile        = str(parse().getValue(filenode))
            #given path of filename, ensure actual path of meshfile is correct
            corrmeshfile    = getpath(filename,meshfile)
            mesh            = getMesh(corrmeshfile,self.spacelu,parse)
            
            self.mlist.append(mesh)

        "get all probes"
	if self.debug: print '  Debug > getStorageObject: getProbes'
        probes             = getProbes(filename,self.spacelu,self.specs.thedoc,self.fp,parse)
        self.plist         = probes.plist
        
        "get all timesteps"
	if self.debug: print '  Debug > getStorageObject: getTimeSteps'
        self.timesteps     = getTimeSteps(filename,self.spacelu,self.specs.thedoc,self.fp,parse)
        self.tlist         = self.timesteps.tlist

        "now set nicknames for variables"

        for mesh in self.mlist:
            for meshspace in mesh.mslist:
                for var in meshspace.vlist:
                    if self.namelu.has_key(var.name): 
                        var.nickname = self.namelu[var.name]
                    else:
                        var.nickname = var.name

        "now update the self.namelu table to to cater for sensitivity variables"

        if len(self.tlist):
            for var in self.tlist[0].vlist:
                self.__nameMap(var.name,'temp_sens_','tsens') 
                self.__nameMap(var.name,'enthalpy_sens_','esens') 
                self.__nameMap(var.name,'volume_fraction_sens_','vfsens') 
            
        for thistime in self.tlist:
            for var in thistime.vlist:
                if self.namelu.has_key(var.name): 
                    var.nickname = self.namelu[var.name]
                else:
                    var.nickname = var.name

        for thisprobe in self.plist:
            for var in thisprobe.vlist:
                if self.namelu.has_key(var.name): 
                    var.nickname = self.namelu[var.name]
                else:
                    var.nickname = var.name

        "get all aborts"
	if self.debug: print '  Debug > getStorageObject: getAborts'
        self.aborts  = getAborts(filename,self.spacelu,self.specs.thedoc,self.fp,parse)

	if self.debug: print '  Debug > getStorageObject: done'


    def __nameMap(self,oldname,oldpattern,newpattern):
        "provide 'newpattern' nicknames for names with 'oldpattern'" 

        if oldpattern in oldname:
            u                         = oldname.replace(oldpattern,newpattern)
            self.namelu[str(oldname)] = str(u.replace('0',''))

    def getValues(self,vars,options=None):

        """
        provides user data values for a list of variables [vars]
        
        we seperate this from obtaining other variable attributes
        so that the data values are obtained only when the user requests them
        """

        for var in vars:
            
	    if self.debug: print '  Debug > gSO.getValues: var.name = ',var.name
            if var.data == None:
                
                "first get variable data in truchas coordinates"
                
	        if self.debug: print '    Debug > gSO.getValues: getData'
                var.getData()
                
                """
                we need to unpermute some data to user coordinates
                """
                
                onmesh = (var.mesh == 'DefaultMesh')
                
                if onmesh and var.rank > 0: # don't bother for rank = 0
                    """
                    unpermute based on the mesh and the meshspace
                    that the variable lives on
                    """
                    
                    for m in range(len(self.mlist)):
                        
                        mesh = self.mlist[m]

                        if mesh.name == var.mesh:
                            """
                            established the mesh this variable lives on 
                            we now can obtain unpermute data for each meshspace
                            and hence alter variable data to be in global
                            context rather than truchas context
                            """
                            if var.name not in ('UNPERMUTE','PERMUTE'):

                                tmp = var.data.copy()
                                
                                if     var.name == 'VERTICES' \
                                   and self.spacelu[var.meshspace] != 'VERTEX':
                                    """
                                    for the VERTICES (connectivity),
                                    variable must obtain map data for both the
                                    vertex space and this particular space
                                    """
                                    for k in range(var.shape[1]):
                                        for j in range(var.shape[0]):
                                            tmp[mesh.updata[var.meshspace][j]-1,k] = \
                                               mesh.updata['VERTEX'][var.data[j,k]-1]

                                elif var.name == 'RHS':
                                    """
                                    Special permutation treatment
                                    for solid mechanics RHS variable
                                    It is of size ndims*nverts
                                    """
                                    ndims  = 3
                                    nsize  = var.shape[0]/ndims
                                    for nd in range(ndims):
                                        for j in range(nsize):
                                            jj  = j*ndims + nd
                                            pj  = mesh.updata[var.meshspace][j]-1
                                            tmp[nd + ndims*pj] = var.data[jj]
                                            
                                else:

				    if var.name not in ('SIDESETS', 'NUMSIDESETS'):
                                        for j in range(var.shape[0]):
                                            tmp[mesh.updata[var.meshspace][j]-1] = \
                                                var.data[j]

                    "now reassign variable data to tmp"

                    var.data = tmp

        return vars

                    
    def str(self):

        print 'name lookup:'
        print self.namelu

        print 'simulation specs:'
        self.specs.str()

        print 'meshes:'
        for mesh in self.mlist:            
            for meshspace in mesh.mslist:
                self.getValues(meshspace.vlist)
            mesh.fillMesh()
            mesh.str()

        print 'timesteps:'
        for timestep in self.tlist:
            self.getValues(timestep.vlist)
            timestep.str()

        print 'probes:'
        for probe in self.plist:
            self.getValues(probe.vlist)
            probe.str()

        print 'aborts:'
        print 'non-linear residual'
        for nlr in self.aborts.nlrlist:
            self.getValues(nlr.vlist)
        print 'linear residual'
        for lr in self.aborts.lrlist:
            self.getValues(lr.vlist)

        self.aborts.str()

if __name__=='__main__':

    'test creating XML storage object'

    from PYTHONutils import uTestOpts
    dfltdir = '../../scripts/test_TBrookParse/samples/'

    opts = uTestOpts('fd', dir=dfltdir )
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)
    
    try:
        xmlstorage = getStorageObject(opt.f,debug=1)
        xmlstorage.str()
    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise

                
                
                                    
            
                
        

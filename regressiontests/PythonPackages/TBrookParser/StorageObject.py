"""

 StorageObject

 -----------------------------------------------------------------------------
  Purpose:

     Provides a generic base class for the datastore - defined as a StorageObject.
     This object class is not intended to be directly used.
     It is used as the base class for (XML)getStorageObject and MODStorageObject classes.

  Public Interface(s):

    None

  Contents

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""

class baseStorageObject:

    def __init__(self):
        "define the storage object"

        self.specs          = None #the simulation specifications
        self.mlist          = []   #a set of meshes
        self.tlist          = []   #a set of timestep structures
        self.aborts         = []   #a set of abort structures
        self.timesteps      = None #a timesteps object (see baseAborts)
        self.plist          = []   #a set of probes
        self.varlist        = []   #a set of user specified variables
        self.namelu         = {}   #maps nicknames to names
        self.spacelu        = {}   #maps nickname spacenames to actual spacenames

        "establish meshspace nickname lookup"

        self.spacelu = { None     :'CELL',
                        'cell'    :'CELL',
                        'cells'   :'CELL',
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
                        'MESH'    :'MESH',
                        'NONE'    :'NONE'
                        }

        "establish variable nickname lookup"


        self.namelu = {'Z_RHO'             :'Density',
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
                      'RHS'                :'Solid_Mechanics_RHS',
                      'Plastic_Strain_Rate':'epsdot',
                      'Thermal_Strain'     :'epstherm',
                      'PC_Strain'          :'epspc',
                      'Displacement'       :'disp',
                      'DISPLACEMENT'       :'DISP',
                      'Gradient_T'         :'GradT',
                      'temp_sens_0001'     :'tsens1',
                      'temp_sens_0002'     :'tsens2',
                      'temp_sens_0003'     :'tsens3',
                      'enthalpy_sens_0001' :'esens1',
                      'enthalpy_sens_0002' :'esens2',
                      'enthalpy_sens_0003' :'esens3',
                      'volume_fraction_sens_00010001':'vfsens11',
                      'volume_fraction_sens_00010002':'vfsens12',
                      'volume_fraction_sens_00020001':'vfsens21',
                      'volume_fraction_sens_00020002':'vfsens22'                       
                      }


    def fillValues(self,vars):
        """
        fills in data values given a list of variables [vars]
        """

    def metastr(self):
        "prints out meta values of attributes in the storage object"
        
        print '/n simulation specs: /n'
        self.specs.str()

        print '/n meshes: /n'
        for mesh in self.mlist:
            mesh.str()

        print '/n timesteps: /n'
        for timestep in self.tlist:
            timestep.str()

        print '/n probes: /n'
        for probe in self.plist:
            probe.str()

        print '/n aborts /n'
        self.aborts.str()


    def str(self):
        "prints out meta + data values of all attributes in the storage object"

        for mesh in self.mlist:
            for meshspace in mesh.mslist:
                self.fillValues(meshspace.vlist)
            mesh.fillMesh()
            
        for timestep in self.tlist:
            self.fillValues(timestep.vlist)

        for probe in self.plist:
            self.fillValues(probe.vlist)

        for nlr in self.aborts.nlrlist:
            self.fillValues(nlr.vlist)
        for lr in self.aborts.lrlist:
            self.fillValues(lr.vlist)
            
        self.metastr()

    def debug(self,verbosity=0):
        """
        for debugging this object;
        two verbosity levels exist:
        1. provide meta values (verbosity=0)
        2. provide meta + data values (verbosity=1)
        """

        if verbosity > 0:
            self.str()
        else:
            self.metastr()
        

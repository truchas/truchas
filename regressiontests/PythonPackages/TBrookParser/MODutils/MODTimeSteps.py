"""

 MODTimeSteps

 -----------------------------------------------------------------------------
  Purpose:

     <purpose>

  Public Interface(s):

   mTSs = modTimesteps(tstorages,currtimesteps)
   mTSs.str()

  Contains:
    class modTimeSteps
        __init__(self,tstorages,currtimesteps,stepids=[0])
        str(self)

    Unit Test Block

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
             Erin Iesulauro Barker (eibarker@lanl.gov)
 -----------------------------------------------------------------------------
"""

if __name__=='__main__':
    " Set sys.path for component testing mode "
    import os, sys
    thisdir   = os.path.dirname(__file__)
    # modify dot-dots appropriately
    parsedir  = os.path.abspath(thisdir+'../') 

    sys.path.append(parsedir)

from BASEutils   import baseTimeSteps
from MODTimeStep import modTimeStep

class modTimeSteps(baseTimeSteps):

    def __init__(self,tstorages,currtimesteps,stepids=[0]):

        "initialise TimeSteps object"
	baseTimeSteps.__init__(self)

	if currtimesteps != None:            
            #append these new timesteps to the existing set of timesteps 'currtimesteps'

            self.idlist     = currtimesteps.idlist
            self.timeslist  = currtimesteps.timeslist    
            self.tlist      = currtimesteps.tlist    
            self.cyclelist  = currtimesteps.cyclelist
            self.vlist      = currtimesteps.vlist    
            self.clist      = currtimesteps.clist    

        #first accumulate all timestep information

        if len(tstorages) > 1:
            
            #mapping multiple T meshes to a single Exodus mesh
            #the new ids and cycles are then simply incremented
            #the new times are set to zero

            #new timestep object created with variables

            thisid       = len(self.idlist)
            thistime     = float(0.0)
            thiscycle    = len(self.cyclelist)
            thisvlist    = tstorages[0].timesteps.tlist[0].vlist
            thisdt       = tstorages[0].timesteps.tlist[0].dt

            self.idlist.append(thisid)
            self.timeslist.append(thistime)
            self.cyclelist.append(thiscycle)

            thistimestep = modTimeStep(id=thisid,cycle=thiscycle,
                                       time=0.0,dt=thisdt,
                                       dtcourant=0.0,vlist=thisvlist)

            self.tlist.append(thistimestep)
                                
        else:

            #mapping from a single Truchas mesh to a single Exodus mesh
            #the new ids, cycles and times are taken from the 'stepids'
            #chosen from the original Truchas simulation
            
            for j in stepids:

                thisid        = len(self.idlist)
                self.idlist.append(thisid)
                self.timeslist.append(tstorages[0].timesteps.timeslist[j])
                self.cyclelist.append(tstorages[0].timesteps.cyclelist[j])
                ttimestep    = tstorages[0].timesteps.tlist[j]
                thistimestep = modTimeStep(id = thisid, cycle = ttimestep.cycle,
                                           time = ttimestep.time,
                                           dt = ttimestep.dt,
                                           dtcourant = ttimestep.dtcourant,
                                           vlist = ttimestep.vlist)
                self.tlist.append(thistimestep)


        #now obtain updated vlist and clist dictionaries 

        for j in self.idlist:

            id             = j
            timestep       = self.tlist[j]
                
            self.vlist[id] = []
            velems         = timestep.vlist
            for velem in velems:
                self.vlist[id].append(velem)
                "initialise the clist dictionary"
                self.clist[velem.name] = []

        cnt = 0
        for id in self.idlist:
            thiscycle = self.cyclelist[cnt]
            for velem in self.vlist[id]:
                self.clist[velem.name].append(thiscycle)
            cnt   = cnt + 1

        
if __name__=='__main__':

    # test modifying timesteps

    dfltdir   = '../../../scripts/test_TBrookParse/samples/'
    T_Sim     = 'map_output/map.TBrook.xml'

    from PYTHONutils import uTestOpts
    opts = uTestOpts('fd',
		     dir      = dfltdir,
                     defaults = {'f': dfltdir + T_Sim})
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    try:

        from XMLgetStorageObject import getStorageObject
        
        modtimesteps  = None
        tstorages     = []

        tstorage      = getStorageObject(opt.f,debug=1)
        tstorages.append(tstorage)
     
        modtimesteps  = modTimeSteps(tstorages,modtimesteps,stepids=[0])
        modtimesteps.metastr()

    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise



        


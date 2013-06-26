"""

 MODTimeStep

 -----------------------------------------------------------------------------
  Purpose:

     <purpose>

  Public Interface(s):

    mTS = modTimeStep(id,cycle,time,dt,dtcourant,vlist)

  Contains:
    class modTimeStep
         __init__(self,id,cycle,time,dt,dtcourant,vlist)
         ungetData(self)
         str(self)

    ---> NO Unit Test Block
         Testing through MODTimeSteps.py unit tests
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

from BASEutils      import baseTimeStep
from CREATEVariable import createVariable
     
class modTimeStep(baseTimeStep):

     def __init__(self,id,cycle,time,dt,dtcourant,vlist):

        "initialise the timestep object"
	baseTimeStep.__init__(self)

        self.id         = id
        self.cycle      = cycle
        self.time       = time
        self.dt         = dt
        self.dtcourant  = dtcourant
        self.vlist      = [] #a list of variable objects   

        #create a set of variables based on the variables in vlist
        for velem in vlist:
             newvar     = createVariable(velem.name,velem.nickname,
                                         velem.rank,velem.shape,velem.type,
                                         None,velem.mesh,velem.meshspace)
             self.vlist.append(newvar)


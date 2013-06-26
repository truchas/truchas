#!/usr/bin/env python
"""
DataCreator

-----------------------------------------------------------------------------
   Purpose:
  
      Creates a DataStore (such as TruchasDataStore). The DataStore contains a
      -simulation storage object (which may contain simulation specs, mesh, timesteps)  
      -list of regions (typically created by the developer)
      -getTimeSteps,getRegion,getField,getProbe methods 
  
   Public Interface:

      T = getDataCreator(dir,map)
      T.getStorage()
      T.getAborts(type)
      T.getRegion(selector,desc)
      T.getField(field,time,cycle,region)
      T.getTimeSteps(timerange,cyclerange,dtstep,cyclestep)
      T.getProbe(name,field,timerange)

   Contains:
      class TruchasDataCreator(dir,map)
            getStorage()
            getAborts(type)
            getRegion(selector,desc)
            getField(field,time,cycle,region)
            getTimeSteps(timerange,cyclerange,dtstep,cyclestep)
            getProbe(name,field,timerange)
            
   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
-----------------------------------------------------------------------------

"""

import os, sys, fnmatch
try:
   import numpy.oldnumeric as Numeric
except ImportError:
   import Numeric
except:
   raise


if __name__ == '__main__':
    print "\n for component test in %s \n" %(__file__)
    thisdir     = os.path.abspath(os.path.dirname(__file__))
    parserdir   = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)

from POSTPROCESSORutils import getRegionObject
from MODutils           import createVariable
from PYTHONutils        import unique, getAbortString

class TruchasDataCreator:
    "defines a creator of a Truchas data store"

    def __init__(self,dir,map=0):

	#initialise the resulting data store that comes out of thisTruchasDataCreator
	self.regions     = []
	self.storage     = None

        self.abortstring = getAbortString()

        #check directory input
        errstring         = '\n \n You did not specify a directory in %s \n \n' %(self.__class__.__name__)
        errstring        += self.abortstring
        assert dir != None, errstring

	#set variables particular to this class
	self.dir           = os.path.abspath(dir)
        self.map           = map
	self.pattern       = '*.TBrook.xml'
	self.filename      = None

    def getStorage(self):
	"factory method to choose a particular way of obtaining a Truchas StorageObject"
	"""
	currently the two choices we have are to 
	    -read in an XML file (XMLgetStorageObject)
	    -map XML fields to an existing exodus mesh (MODStorageObject)
	"""
	for filename in os.listdir(self.dir):
            if fnmatch.fnmatch(filename, self.pattern):
		self.filename = filename

	errstring    = '\n \n No *.TBrook.xml file in %s directory \n \n' %(self.dir) 
        errstring   += self.abortstring
	assert self.filename != None, errstring

        self.xmlfile = os.path.join(self.dir,self.filename)

        if self.storage == None:

            fp = open("/dev/null", 'w')
            #check if we get the storage object from mapping or from reading in XML file
            if self.map :
                from MODStorageObject import modStorageObject  
                self.storage = modStorageObject(self.xmlfile) 
            else:
                from XMLgetStorageObject import getStorageObject
                self.storage = getStorageObject(self.xmlfile,fp)
            fp.close()

        #now ensure mesh field values (such as centroids, posns) are filled in
        for mesh in self.storage.mlist:
            for meshspace in mesh.mslist:
                self.storage.getValues(meshspace.vlist)
            mesh.fillMesh()


        self.meshes = self.storage.mlist

        #check smallest dt value...get time input tolerance based on this value
        dtmin   = 10000000.0
        if len(self.storage.tlist) > 1:
           lastime = self.storage.tlist[0].time
           for thistep in self.storage.tlist[1:]:
              dtmin   = min(dtmin,thistep.time-lastime)
              lastime = thistep.time 

        self.dttol = 0.10*dtmin
                    
        return 

    def getAborts(self,type='NONLINEAR'):

        if type.upper() == 'NONLINEAR':
           aborts = self.storage.aborts.nlrlist
           for nlr in aborts:
              self.storage.getValues(nlr.vlist)

        if type.upper() == 'LINEAR':
           aborts = self.storage.aborts.lrlist
           for lr in aborts:
              self.storage.getValues(lr.vlist)

        return aborts

    def getRegion(self,selector=['ID','CELL','1-10'],desc=None):
	"user employs this function to create a region based on criteria - returns a region object" 

        newselector = selector
        
        #wrapper to redefine the selector for 'VAR' cases -
        #developer provides variable name rather than the actual variable 
        if selector[0] == 'VAR':
            
            thevar     = None

            nostepserrstring   = "\n \n           In region creation with the 'VAR' selector.                 \n"
            nostepserrstring  += "            No timestep output exists in %s                              \n" %(self.filename)
            nostepserrstring  += "            Please check the BasicRun.log file for successful completion \n"
            nostepserrstring  += "            Also ensure XML_Output_Dt_Multiplier > 0 in your input fil e \n"
            nostepserrstring  += self.abortstring

            assert len(self.storage.tlist) > 0, nostepserrstring

            for v in self.storage.tlist[0].vlist:
                if v.name == selector[1] or v.nickname == selector[1]:
                    thevar = v

            varerrstring  = "\n \n             In region creation with the 'VAR' selector. No variable with   \n"
            varerrstring += "             the name or nickname = %s in %s                                 \n" %(selector[1],self.filename)
            varerrstring += "             Please check the BasicRun.log file for a list of valid variables \n"
            varerrstring += self.abortstring

            assert thevar != None, varerrstring

            newselector[1] = thevar 

        region = getRegionObject(userinput = None,
                                 name      = desc,
                                 file      = self.xmlfile,
                                 storage   = self.storage,
                                 mesh      = self.storage.mlist[0],
                                 selector  = newselector)

        #we will automatically get the region indices for the first cycle (i.e cycleid = 0) 
        #for cycleids > 0 we invoke region.getIndices in the getField method
        
	self.regions.append(region)

	return region 


    def getField(self,field,time=None,cycle=None,region=None):
        "returns a field (variable) at either time or cycle with data defined on a particular region"

        nostepserrstring  = "\n \n        No timestep output exists in %s  \n" %(self.filename)
        nostepserrstring  += "            Please check the BasicRun.log file for successful completion \n"
        nostepserrstring  += "            Also ensure XML_Output_Dt_Multiplier > 0 in your input file \n"
        nostepserrstring  += self.abortstring

        assert len(self.storage.tlist) > 0, nostepserrstring

        thestep = None
        
        if time == None and cycle==None:
            #get field from first output timestep
            thestep = self.storage.tlist[0]
        
        if time == None and cycle != None:
            #get timestep based on cycle number
            for thistep in self.storage.tlist:
                if thistep.cycle == cycle:
                    thestep = thistep

            cycerrstring  = "\n \n        Cycle = %d does not exist in %s \n" %(cycle,self.filename)
            cycerrstring += "            Please check the BasicRun.log file for a list of valid cycles \n"
            cycerrstring += self.abortstring

            assert thestep != None, cycerrstring

        if time != None and cycle == None:
            #get timestep based on time value

            for thistep in self.storage.tlist:
                if abs(thistep.time-time) < self.dttol:
                    thestep = thistep
                    
            timerrstring  = "\n \n        Time = %12.5e with tolerance = %12.5e does not exist in %s \n" %(time,self.dttol,self.filename)
            timerrstring += "            Please check the BasicRun.log file for a list of valid times \n"
            timerrstring += self.abortstring

            assert thestep != None, timerrstring

        meshfieldslu = {'CENTROIDS'  :'CENTROIDS',
                        'VOLUMES'    :'VOLUMES',
                        'CELLVOLUMES':'VOLUMES',
                        'COORDS'     :'COORDS',
                        'NODES'      :'COORDS',
                        'VERTEX'     :'COORDS',
                        'VERTICES'   :'COORDS',
                        'VERTEXES'   :'COORDS'}

        if meshfieldslu.has_key(field.upper()):
            #the field is a mesh variable rather than a physical variable

            tname   = meshfieldslu[field.upper()]
            mesh    = self.storage.mlist[0]
            for meshspace in mesh.mslist:
               if meshspace.name == 'VERTEX': 
                  for var in meshspace.vlist:
                     if var.name == tname:
                        thevar  = var
               elif meshspace.name == 'CELL': 
                  for var in meshspace.vlist:
                     if var.name == tname:
                        thevar  = var
        else:
            #the field is a physical variable
        
            thevar = None
            for v in thestep.vlist:
                if v.name == field or v.nickname == field:
                    thevar = v

            varerrstring  = "\n \n        No variable with the name or nickname = %s in %s\n" %(field,self.filename)
            varerrstring += "            Please check the BasicRun.log file for a list of valid variables \n"
            varerrstring += self.abortstring

            assert thevar != None, varerrstring

            self.storage.getValues([thevar])

        #if a region is specified get the field values on this region, otherwise
        #get the field values on the entire mesh

        newvar       = createVariable(thevar.name,
                                      thevar.nickname,
                                      thevar.rank,
                                      [],
                                      thevar.type,
                                      thevar.data,
                                      thevar.mesh,
                                      thevar.meshspace)

        #the field variable returned from this method has 2 additional attributes - stepid, cycle and time
        newvar.stepid = thestep.id
        newvar.cycle  = thestep.cycle
        newvar.time   = thestep.time

        if region != None:
            #get field with data defined on the region only
            region.getIndices(region.selector,cycleid=thestep.id)
            theindices = region.indices[thevar.meshspace][thestep.id]

            if len(theindices):
                shape           = []
                for i in range(len(thevar.shape)):
                    shape.append(thevar.shape[i])
                shape[0]        = len(theindices)
                newvar.shape[:] = shape[:]
                newvar.data     = Numeric.array(0.0,thevar.type)
                newvar.data     = Numeric.resize(newvar.data,newvar.shape)
                for i in theindices:
                    newvar.data[theindices.index(i)] = thevar.data[i-1]
                
            else:
                newvar.shape[:] = thevar.shape[:]
                newvar.shape[0] = 0
                newvar.data     = None
        else:
            #get field with data defined on the entire mesh
            newvar.shape[:] = thevar.shape[:]
            newvar.data     = thevar.data

        return newvar

    def getProbe(self,name,field,timerange=None):
        "returns values of a field at a probe location within a specified timerange"

        probes       = self.storage.plist
        #thevar       = None
        theprobevar  = None
        theprobename = None

        for probe in probes:
            if probe.name == name:
                theprobename = name 
                #found the probe now get probe field values...

                for var in probe.vlist:
                    if var.name == field or var.nickname == field:
                        self.storage.getValues([var]) #([thevar])                       
                        theprobevar = createVariable(var.name,
                                                     var.nickname,
                                                     var.rank,
                                                     var.shape,
                                                     var.type,
                                                     var.data,
                                                     var.mesh,
                                                     var.meshspace)
                        #thevar          = var
                        #theprobevar.str()
                        #theprobevar.data[:,:] = var.data[:,:]

        proberrstring  = "\n \n        No probe with the name = %s in %s \n" %(name,self.filename)
        proberrstring += "            Please check the BasicRun.log file for a list of valid probes \n"
        proberrstring += self.abortstring

        varerrstring   = "\n \n        Probe %s has no field variable with the name or nickname = %s in %s\n" %(name,field,self.filename)
        varerrstring  += "            Please check the BasicRun.log file for a list of valid variables \n"
        varerrstring  += self.abortstring

        assert theprobename != None, proberrstring
        assert theprobevar  != None, varerrstring

        timerrstring   = '\n \n You did not specify the timerange in the form [a,b] in getProbe method \n \n'
        timerrstring  += self.abortstring

        tmp     = Numeric.array(0.0,theprobevar.type)
        tmp     = Numeric.resize(tmp,[len(theprobevar.data)])
        tmp[:]  = theprobevar.data[:,1] 
        mintime =  1.0e+14
        maxtime = -1.0e+14
        mintime = min(tmp[Numeric.argmin(tmp)],mintime) 
        maxtime = max(tmp[Numeric.argmax(tmp)],maxtime)

        tolerance = abs(maxtime-mintime)*0.000001
        mintime  -= 1.5*tolerance
        maxtime  += 1.5*tolerance
        tmp[:]   += tolerance #to ensure non-zero entries in tmp array

        if timerange != None:
            if '[' == timerange[0] and ']' == timerange[len(timerange)-1]:
                timerange = timerange[1:len(timerange)-2]
                timerange = timerange.split(',')

            tolerance = abs(maxtime-mintime)*0.000001
            mintime   = float(timerange[0])-tolerance
            maxtime   = float(timerange[1])+tolerance

            assert len(timerange)==2, timerrstring

        x       = Numeric.choose(Numeric.less_equal(tmp,mintime),(tmp,0))
        y       = Numeric.choose(Numeric.greater_equal(x,maxtime),(x,0))
        nz      = Numeric.nonzero(y)

        theprobevar.data  = theprobevar.data[nz[0]:nz[len(nz)-1]+1,0:] #2:]

        return theprobevar 
 
            
    def getTimeSteps(self,timerange=None,cyclerange=None,dtstep=None,cyclestep=None):
	"returns a list of Truchas timesteps with timestep.time in timerange or timestep.cycle in cyclerange"

	timesteps     = self.storage.tlist

        timerrstring   = '\n \n You did not specify the timerange in the form [a,b] in getTimeSteps method \n \n'
        timerrstring  += self.abortstring 

        cyclerrstring  = '\n \n You did not specify the cyclerange in the form [a,b] in getTimeSteps method \n \n'
        cyclerrstring += self.abortstring

        #if no timerange or cyclerange specified - return the entire timestep range
        if timerange == None and cyclerange == None:
            reqdsteps = timesteps

        #if timerange specified....
        if timerange != None:
            assert len(timerange) == 2, timerrstring
            reqdsteps = []
            for timestep in timesteps:
                if timestep.time <= (timerange[1]+self.dttol) and (timestep.time >= timerange[0]-self.dttol):
                     if (dtstep != None) and (timestep.time % dtstep < self.dttol): 
                        reqdsteps.append(timestep)

        #if cyclerange specified....
        if cyclerange != None:
            assert len(cyclerange) == 2, cyclerrstring
            reqdsteps = []
            for timestep in timesteps:
                if timestep.cycle <= (cyclerange[1]) and timestep.cycle >= (cyclerange[0]):
                    if (cyclestep != None) and (timestep.cycle % cyclestep == 0): 
                        reqdsteps.append(timestep)
        
	return reqdsteps

def getDataCreator(dir=None):
    "factory method to choose a particular DataCreator based on dir value - right now we have only a TruchasDataCreator"

    return TruchasDataCreator(dir)

if __name__=='__main__':
    testdata  = getDataCreator(dir='../TestCases/static_drop/static_drop_golden')
    testdata.getStorage()
    #testdata.storage.str()
    timesteps = testdata.getTimeSteps(cyclerange=[0,1],cyclestep=1)
    region    = testdata.getRegion(selector=['ID','VERTEX','50-80'],desc='50<CELL IDs<80')
    newregion = testdata.getRegion(selector=['VAR','Density',[0.95,1.0]],desc='0.95<Density<1.0')
    newregion = testdata.getRegion(selector=['S','[0.0,0.0,0.0->0.75,0.75,0.75]'],
                                   desc='Box region:[0.0,0.0,0.0->0.75,0.75,0.75]')
    newregion.getIndices(newregion.selector,cycleid=1)
    newregion.str()
    thefield  = testdata.getField(field='Velocity',time=0.00105,region=newregion)
    thefield.str()

    mapbdir   = thisdir + '/../../../scripts/test_TBrookParse/samples/map-b_output'
    testdata  = getDataCreator(mapbdir)
    testdata.getStorage()
    newregion = testdata.getRegion(selector=['MB',4,5],desc='Mesh Block Number = 4')
    for region in testdata.regions:
        region.str()


    



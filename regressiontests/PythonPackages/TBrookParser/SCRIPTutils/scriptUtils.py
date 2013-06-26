

"""
scriptUtils.py
written by Sriram Swaminarayan

This module has utilities that will allow one to
parse the TBrook XML files and manipulate the data.

provides classes suFile and suMesh

try:
  from scriptutils import *
  f = suFile('path/to/file.TBrook.xml')
  f.help()

  m = suMesh(f.meshes[0],f)
  m.help()
  
"""


import os, sys
try:
    from pathmagic import *
except:
    print """
    Unable to import path magic.
    Did you move this script?
    """
    sys.exit(2)
    
try:
    import numpy.oldnumeric as Numeric
except ImportError:
    import Numeric
except:
    print "Unable to import numpy or Numeric"
    sys.exit(2)
    
class suVariable:
    """
    The class for variables created by the suFile function createVariable.
    We subclass (or is it superclass) the brookVariable class overriding
    only the create and getData() functions.


    Class attributes are:
      name        name of variable
      nickname    same as name
      type        i/r/d/l/c
      rank        rank of variable
      data        pointer to data
      mesh        name of mesh
      meshspace   space on which variable lives (cell/vertex/...)
      element     None
      offset      None
      arrorder    None
      parse       None
      dformat     None
      tfile       None

    """


    def __init__(self,name, data=[], type=None, 
                 map='CELL',meshName='DefaultMesh'):
        """
        Initializes a variable.
        The only required arument is name.
        """

        self.typeLookup={'integer':'i', 'int':'i','i':'i',                         
                         'real':'f', 'float':'f','r':'f','f':'f',
                         'double precision':'d','real*8':'d','double':'d','d':'d',
                         'logical':'l', 'bool':'l','l':'l',
                         'string':'c', 'character':'c','char':'c','c':'c',
                         None:'f'
                         }
        self.name         = name
        self.nickname     = name
            
        self.rank         = None
        self.data         = None
        self.mesh         = meshName
        self.meshspace    = map
        self.element      = None
        self.offset       = None
        self.arrorder     = None
        self.parse        = None
        self.dformat      = None
        self.tfile        = None

        self.setData(data,type)
    
    def getMETAVALUES(self):

        "return all values describing the variable except for the actual data values" 
        return self.name,self.type,self.rank,self.shape

    def setData(self, data, type):
        """
        The method for setting the data for the variable
        """
        
        if ( data == None ): return


        self.rank   = 0
        self.shape  = []

        t = data;
        try:
            while (1):
                self.shape = self.shape + [ len(t) ]
                t = t[0]
                self.rank = self.rank + 1
        except:
            # do nothing
            aHack = None
                
        if(len(self.shape) == 0 ):
            self.shape = [1]
            self.rank  = 1

        
        self.type = self.typeLookup[type]
        self.data   = Numeric.array(data,self.type)

        return

    def getValues(self):
        """
        Returns a pointer to self.data
        """
        return  self.data


    def getData(self):
        """
        Returns a pointer to self.data
        """
        return  self.data
        
    def help(self):
        """
        prints this message
        """
        print help(self)

class suMesh:
    """
    The mesh structure for class suFile.
    This class is a baby mesh that has information I think is
    of interest to the general person. You initialize it with a
    XMLMesh structure (from suFile.meshes[0] or some such).
    You can then directly access the different aspects of the
    mesh to manipulate whatever quantity you wish to modify

    Creataion: It is generally created in the following manner:
        from scriptUtils import *
        f = suFile('path/to/xmlFile.TBrook.xml')
        myMesh = suMesh(f.meshes[0], f)

    The structure has the following components:
        name: name of the mesh
        type: HEX or TET or whatever
        ncells: number of cells
        nnodes: number of nodes
        nfaces: number of faces (not initialized)
        nedges: number of edges (not initialized)

        cellCentroids:      double list of centroids of cells
        cellConnectivity:   what vertices make up this cell in truchas order
        cellVolumes:        the volume of the cell
        cellNumNeighbors:   number of neighbors of this cell
        cellBoundaryFlag:   boundary flag put out by truchas

        nodeCoords:   coordinates of the different nodes
        nodeRSumRVol: the RSumRVol value of the node as computed by truchas

    Unfortunately, little if any error checking is done
    """
    name=''
    type=''
    ncells=0
    nnodes=0
    nedges=0
    nfaces=0
    mesh=None
    cellCentroids=[[]]
    cellConnectivity=[[]]
    cellVolumes=[[]]
    cellNumNeighbors=[[]]
    cellBoundaryFlag=[[]]
    nodeCoords=[[]]
    nodeRSumRVol=[[]]

    def help(self):
        print help(suMesh)
        return

            
    def __init__(self,mesh,suFile):
        """
        Initializes a suMesh structure by reading in values from the
        mesh file and filling them in.

        suFile is used as the source and mesh is the mesh of interest.
        """
        if ( suFile.storage == None ): return # no file means no mesh
        if ( mesh == None ) : return # no mesh means no suMesh
        self.name = mesh.name
        self.mesh = mesh
        #first get all data values associated with this mesh
        for x in mesh.mslist:
            name = suFile.storage.spacelu[x.name]
            suFile.storage.getValues(x.vlist)
            if(name == 'CELL'):
                mesh.fillCells(x)
                self.ncells           = x.size
                self.cellVolumes      = mesh.cells['VOLUMES']
                self.cellCentroids    = mesh.cells['CENTROIDS']
                self.cellConnectivity = Numeric.array(mesh.cells['VERTICES'])
                self.cellNumNeighbors = mesh.cells['NUMNEIGHBORS']
                self.cellBoundaryFlag = mesh.cells['BOUNDARY']
            elif(name == 'FACE'):
                mesh.fillFaces(x)
                self.nfaces = x.size
                if(self.nfaces == None): self.nfaces=0
            elif(name == 'VERTEX'):
                mesh.fillVertices(x)
                self.nnodes       = x.size
                self.nodeCoords   = Numeric.array(mesh.vertices['COORDS'])
                self.nodeRSumRVol = mesh.vertices['RSUMRVOL']
            elif(name == 'EDGE'):
                mesh.fillEdges(x)
                self.nedges = x.size
                if(self.nedges == None): self.nedges=0


class suFile:
    """
    The basic class that has the functions for
    navigating the XML file.
    
    Help on class scriptUtils.suFile
    
    This is the heart of the scriptUtils module and provides an easy interface
    to the TBrook xml files.  I expect people will use this as an entry point
    into the many different modules written for the TBrookParser.
    
    Typically one would initialize a structure are:
      from scriptUtils import *
      f = suFile('path/to/file.TBrook.xml')

      f then has the following properties:
           f.init:          Is the suFile class initialized?
           f.file:          name of xml file
           f.meshes:        list of meshes in f ( see also the suMesh class)
           f.storage:       storage object with the actual pointers to data
           f.times:         list of times of different timesteps
           f.cycles:        list of cycle IDs corresponding to those times
           f.timesteps:     list of timestep structures for above times and cycles
           f.varCycles:     dictionary of variables and cycles they appear in

        and Methods:
           f.findVariable(ts, name): returns variable of name (or nickname) 'name' in timestep ts
           f.copyVariable(var):      returns a copy of variable var that you can modify
           f.newVariable(rank, name, nickname,data):returns a new variable
           f.findClosestTimestep(t): returns the timestep in f closest to time t
           f.getValues(varlist):     gets values for all variables in list varlist
        --
    """
    
    file=None
    storage=None
    times = []
    timesteps=[]
    cycles = []
    regions={}
    init = 0
    meshes=[]
    varCycles = []
    ppobject   = getPostProcessorObject('','',None,debug=0)

    def help(self):
        print help(suFile)
    def __init__(self, file):
        """
        Read the file, and load data into an internal data structure.
        """
        Estring="file assign error"
        try:
            self.file = file
            Estring="ppobject load error"
            X  = self.ppobject.load(self.file)
            Estring="get storage object error"
            self.storage = self.ppobject.files_storages[X.file][0]
        except:
            self.init = 0
            print
            print "ERROR loading",file
            print "    "+Estring
            print 
            return
        self.meshes = self.storage.mlist
        self.timesteps = self.storage.timesteps.tlist
        self.cycles = self.storage.timesteps.idlist
        self.varCycles = self.storage.timesteps.clist
        self.times = self.storage.timesteps.timeslist
        self.tsteps = None
        self.init = 1

    def copyVariable(self,v,name=None):
        """
        Returns a copy of variable v that can be modified wihout
        affecting the original variable.

        Both name and nickname are set to the input variable name
        If name is not supplied, we append '_copy' to the name
        and nickname of the original variable
        """
        
        import copy
        import XMLutils
                                
        if(name == None):
            tmpName =  v.name+'_copy'
        else:
            tmpName = name

        self.storage.getValues([v])
        tmpData = v.data.tolist()
        tmpVar = self.createVariable(name,tmpData,v.meshspace,v.mesh, v.type)
        return(tmpVar)

    def createVariable(self, name, data, map='CELL', mesh='DefaultMesh', type='d'):

        """
        Creates a variable with name name, data, mesh, and type

        type is one of: 
           'integer','int','i'
           'real','float','f'
           'double precision','real*8','double','d'
           'character','string','char','c'
           'logical','bool','l'

        type defaults to None if none is specified

        Note that type is not used for any of the graphics write.
        Only the restart writer uses the type, but you shouldn't
        be messing with the variables for the restart writer anyway!

        """
        tmpVar = suVariable(name,data=data, map=map, meshName=mesh, type=type)
        return tmpVar

    def findClosestTimestep(self,t):
        """
        Find timestep in suFile that is closest to t
        """
        if(not self.init): return(None)

        tsNow = None
        tsPrev = self.timesteps[0]
        tPrev = tsPrev.time
        for tsNow in self.timesteps:
            tNow = tsNow.time
            if(tNow == t): return(tsNow)
            if(tNow >= t):
                if(tNow-t > t-tPrev):
                    return(tsPrev)
                else:
                    return(tsNow)
        # If we are here, then it means that
        # the highest t in the suFile was less
        # than t.  So, return tsNow
        #
        # This also has the desirable effect of returning None if there
        # are no timesteps in suFile
        return(tsNow)
    def getValues(self,varList):
        if(not self.init): return(None)
        self.storage.getValues(varList)
        
    def findVariable(self,ts,name):
        if(not self.init): return(None)
        for v in ts.vlist:
            tmpVar=None
            if(v.nickname == name or v.name == name):
                tmpVar = v
                break
        if(tmpVar == None):
            print "Unable to find variable"
            return(None)

        try:
            self.storage.getValues([tmpVar])

        except:
            print "Unable to get data for variable.  Perhaps a different map?"
        return(tmpVar)




def getLegacy(filename,outFile=None,mode='o'):
    """

    In GetLegacy
    Obtains the legacy strings.  Default is to get the .out file.
    mode='o' gets .out, 'e' gets .err, and 't' gets 'tty' output.
    Returns legacy string on success, None on failure.

    Legacy string is written to outFile, if outFile is provided.
    """
    #
    # Initialize some variables first
    import sys
    ll   = [];
    tags = []
    idx  = 0
    l    = ''
    lstr = ''

    if ( mode == 'o'):
        legacyTag='OUT'
    elif (mode == 'e'):
        legacyTag='ERR'
    elif (mode == 't'):
        legacyTag='TTY'
    else:
        print """
        Invalid mode in getLegacy.
        Valid values are: 'o', 'e', and 't'.
        """
        return None

    try:

        Estring = 'file read error, file='+filename
        print 'reading',; sys.stdout.flush()

        fp    = open(filename,'r')
        lines = fp.readlines()
        fp.close()

        flag   = 0
        nlines = len(lines)

        Estring = 'file process error, file='+filename
        for l in lines:
            # flag determines wheter a line will  be included
            # or not.
            # flag = 0: include line
            # flag = 1: legacy line
            # flag = 2: inside a <CYCLE> tag
            #

            idx += 1

            if(idx%(nlines/10) == 0):
                print '.',;
                sys.stdout.flush()
            if(flag == 0 and l.startswith('[')):
                if (l.find(legacyTag)):
                    flag = 0
                else:
                    flag = 1
                continue

            if (flag < 2 ) and (l.startswith('<')):
                flag = 0


            if( (len(l) > 1) and l.startswith('<')):
                if(l.startswith('<?xml') or l.startswith('<!')):
                    # do not add xml and cdata tag lines
                    donotaddtotagshack=' '
                    continue
                elif (l.startswith('</')):
                    # pop the tags array if it exists
                    if (len(tags)):
                        tn = tags.pop()
                elif(l.startswith('<') and (not l.endswith('>\n'))):
                    # an incomplete tag
                    # do nothing only if last line of file
                    if(idx == nlines):
                        # last line of file
                        continue
                elif (l.startswith('<') and l.endswith('>\n') and (not l.find('</')>0)):
                    # Only add tags that do not end on the same line
                    tagname = (l[1:-1][:l.find(' ')])
                    if(tagname.endswith(' ')): tagname = tagname[:-1]
                    tags.append(tagname)
                continue

            if((flag == 0) and (len(tags) == 1) ):
                if ((l.find('      END') < 0 ) and ( not l.endswith('                      '))):
                    ll.append(l)
                else:
                    ll.append('\n')

        Estring = 'string join error'
        lstr  = ''.join(ll)

        lines = None
        l     = None
        ll    = None

        if ( outFile != None ):
            Estring = 'File write error'
            fpo = open (outFile,'w')
            fpo.write(lstr)
            fpo.close()

    except:
        print
        print
        print '___PROCESSING ERROR:'
        print '    ',Estring
        print '     file =',filename
        print '     line =',idx,l
        print
        print

    return lstr








if __name__ =='__main__':

    # Load the file
    # Note that we depend on pathdir from pathmagic
    f = suFile(pathdir+'/../../scripts/test_TBrookParse/samples/phasechange_mixed_output/phasechange_mixed.TBrook.xml')
   
    # Get the primary mesh structure
    mesh = f.meshes[0]  # The primary mesh file, note this is 
                        # what the different output writers want
   
    m    = suMesh(mesh, f) # a convenience structure
   
    # Info on what timesteps are in file f
    print "\n  Summary of timesteps in",f.file
    print "   cycle: t  H_total"
    for ts in f.timesteps:
        h_total= 0.0
        hVar = f.findVariable(ts,'Enthalpy')
        rhoVar = f.findVariable(ts,'Density')

        for i in range(m.ncells):
            h_total = h_total + hVar.data[i]*rhoVar.data[i]*m.cellVolumes[i]
        print "  %6d: %.2f %.2f"%(ts.cycle, ts.time,h_total)

    # Info on variables in the last timestep
    # Note that we use ts from above loop
    #
    # Pick up the enthalpy and temperature variables
    # Note we use nicknames.  You can check nicknames
    # by looping through vlist and printing out the
    # names and nicknames
    hVar = f.findVariable(ts,'Enthalpy')
    rhoVar = f.findVariable(ts,'Density')
    TVar = f.findVariable(ts,'T')
    
    # Create a copy of the enthalpy variable
    cpVar = f.copyVariable(hVar,'Cp') # make a deepcopy of the hVar
    cpVar.data = hVar.data/TVar.data
   
    #
    # Write a GMV file
    fileout   ='test.gmv'    # The output file name
    outputfmt = 'ascii'      # The output format
    seq_no    = 1            # The sequence number
    vars = [hVar, cpVar, TVar]
    flags     = { 'testRegion':[1,2,3,4,99], #comma separated dictionary of cell groups
              'test2':[4,20,35,56]        #that will show up in groups in the GMV file
              }
    # Write our variables to an ascii GMV file
    from GMVwriteStorageObject import GMVwriteStorageObject
    GMVwriteStorageObject('test.gmv', outputfmt, f.meshes[0],
                          ts.time, ts.cycle, vars, seq_no, flags)            
   


    # Write a binary restart file:
    f.storage.getValues(ts.vlist) # Important to get values beforehand
    from RESTARTwriteStorageObject import RESTARTwriteStorageObject
    RESTARTwriteStorageObject('test.rst', 'binary', f.storage, mesh,
                              0, ts.time, ts.vlist, {})


    # Some helpful commands
    # uncomment for help.
    #
    # f.help()
    # m.help()
    # print help(f)
    # print help(m)


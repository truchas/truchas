

"""

 getRegionObject

 -----------------------------------------------------------------------------
  Purpose:
  
     To instantiate and provide methods for Region objects.
  
  Public Interface(s):

    ro = getRegionObject(userinput)
    i  = ro.getIndices(selection)                               
    ro.str()
    
  Contains:
    class getRegionObject
        __init__(self,userinput,file=None,storage=None,
                 mesh=None,name= '',selector=['Default'],debug=0)
        getIndices(self,selection,cycleid=0)
          __getIndicesFromPromptInterface(self,selection,cycleid)
          __getIndicesFromDeveloper(self,criteria,cycleid)
          __getIDIndices(self,input,ninds)
          __getVARIndices(self,thisvar,cycleid,lower,upper,dim=None,face=None)
          __getMBIndices(self,mblock)
          __getSIndices(self,mspace,posns)
        str(self)                                                          

    Unit Test Block
  
  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""
import os, sys, string
if __name__ != '__main__':
    from PYTHONutils import unique, getAbortString

    thisdir   = os.path.dirname(__file__)
    parsedir  = os.path.abspath(thisdir+'../') 
    sys.path.append(parsedir)

try:
   import numpy.oldnumeric as Numeric
except ImportError:
   import Numeric
except:
   raise

class getRegionObject:

    def __init__(self,userinput,file=None,storage=None,
                 mesh=None,name= '',selector=['Default'],fpwatch=sys.stdout,debug=0):

        "initialise the region object"

        self.name           = name      # name of the region
        self.file           = file      # filename to use for creating region
        self.storage        = storage   # storage object to use for creating region
        self.mesh           = mesh      # mesh that this region will live on
        self.meshspace      = []        # meshspace(s) the region will live on 
        self.selector       = selector  # list of possible methods to create a region
        self.input          = userinput # user input mechanism
        self.describe       = ''        # string descriptor of the region
        self.indices        = {}        # provides the indices that live on the region for all valid meshspaces
        self._IDMBSindices  = {}        #local variable that provides indices from the 'ID,MB and S' selector(s)
        self._var0indices   = []        #local variable that provides indices from the 'VAR' selector at cycle 0        

        "now start defining the above attributes in the region object"

        self.debug            = debug
        self.fpwatch          = fpwatch

        self.selectorCriteria = ''
        self.varformat        = 3

        nonevalid             = 1

        if self.input == None:
            "region defined by developer inputted list of the form [ID|VAR|S|MB,*]"

            selectionisvalid = self.getIndices(selector[0])
            nonevalid        = nonevalid*(not selectionisvalid)

        else:
            "region defined by prompt interface"

            for meshspace in self.mesh.mslist:
                if not self.indices.has_key(meshspace.name):
                    self.indices[meshspace.name]          = {}

                if not self.indices[meshspace.name].has_key(0):
                    self.indices[meshspace.name][0]       = []

                if not self._IDMBSindices.has_key(meshspace.name):
                    self._IDMBSindices[meshspace.name]    = {}

                if not self._IDMBSindices[meshspace.name].has_key(0):
                    self._IDMBSindices[meshspace.name][0] = []   

            for i in self.selector:

                if 'Default' in i:
            
                    #creating a default region - required when a file is first loaded
                    selection   = 'Default'
                    cnt         = 1
                    for m in self.storage.mlist:
                        if self.mesh.name == m.name:
                            count         = cnt
                        cnt              += 1
                        self.name         = 'allMESH%i' %(count)

                        selectionisvalid  = 1
                        selectionisvalid  = self.getIndices(selection)

                        tmp               = self.file.split('/')
                        fname             = str(tmp[len(tmp)-1])

                        mspaces           = []
                        for ms in m.mslist:
                            mspaces.append(ms.name)
                        self.meshspace.append(mspaces)
                        self.selectorCriteria = 'all indices defined'


                elif 'ID' in i:
                    
                    #creating a region based on user index ID input
                    #- need the relevant meshspace
                    selection        = 'ID'
                    selectionisvalid = 1
                    selectionisvalid = self.getIndices(selection)

                elif 'VAR' in i:

                    #creating a region based on variable values
                    selection        = 'VAR'
                    selectionisvalid = 1
                    selectionisvalid = self.getIndices(selection)

                elif 'S' in i:

                    #creating a region based on spatial extents
                    selection        = 'S'
                    selectionisvalid = 1
                    selectionisvalid = self.getIndices(selection)


                elif 'MB' in i:

                    #creating a region based on mesh block choice
                    selection        = 'MB'
                    selectionisvalid = 1
                    selectionisvalid = self.getIndices(selection)

                else :
                
                    print >> self.fpwatch
                    print >> self.fpwatch, '%s is not an appropriate selector ' %(i)
                    print >> self.fpwatch, 'Please choose a combination of [ID, VAR, S, MB]'
                    print >> self.fpwatch
                    return

                nonevalid            = nonevalid*(not selectionisvalid)

                if selectionisvalid:

                    if len(self.meshspace):
                        themeshspace = self.meshspace[-1] # the last one
                    tmp              = self.file.split('/')
                    fname            = str(tmp[len(tmp)-1])
                    self.describe   += ' %-20s'%self.name + \
                                       ': defined for ' + fname + \
                                       ' on mesh ' + self.mesh.name + \
                                       ' as '+str(themeshspace)+ \
                                       '(s) with '+self.selectorCriteria + '\n'


                    #care must be taken when combinations of ID,MB,S and VAR are chosen in the
                    #prompt driven interface - for a given cycle we need to keep the ID,MB,S
                    #indices and update the VAR indices

                    if 'Default' not in self.selector:
                        for meshspace in self.mesh.mslist:
                            length = len(self._IDMBSindices[meshspace.name][0])
                            L      = self._IDMBSindices[meshspace.name][0][0:length]
                            L      = unique(L)
                            L.sort()
                            self.indices[meshspace.name][0] = L

            
                    if 'VAR' in selection:
                        for meshspace in self.mesh.mslist:
                            if meshspace.name == self.var.meshspace:
                                length = len(self._IDMBSindices[meshspace.name][0])
                                L      = self._IDMBSindices[meshspace.name][0][0:length]
                                L.extend(self._var0indices)
                                L = unique(L)
                                L.sort()
                                self.indices[meshspace.name][0] = L
                
        if nonevalid:
            """
            entire region selection is not valid so discard all
            properties of this region
            """
            self.name        = ''
            self.file        = ''
            self.storage     = None
            self.mesh        = None
            self.meshspace   = None
            self.selector    = None
            self.describe    = ''
            self.indices     = None

    def getIndices(self,selection,cycleid=0):
        """
        factory method that gets indices of the region
        - depends on the form of the user inputted selection criteria
        """
        if self.input == None:
            "selection from developer"
            validselection = self.__getIndicesFromDeveloper(self.selector,cycleid)

        else:
            "selection from prompt driven interface"
            validselection = self.__getIndicesFromPromptInterface(selection,cycleid)

        return validselection

    def __getIndicesFromPromptInterface(self,selection,cycleid):
        "get indices when selector is the prompt driven interface"

        if cycleid == 0:
            """
            when first defining a region get all indices from all
            chosen selection criteria
            """
            
            if selection == 'Default':
                #provide entire indices for all meshspaces on the chosen mesh

                for meshspace in self.mesh.mslist:

                    if meshspace.size != None:
                        x = meshspace.size+1
                        L = range(x)
                        T = L[1:len(L)]
                        self.indices[meshspace.name][cycleid] = T                    

            if selection == 'ID':
                #provide indices based on user input

                prompt = 'Region based on  "CELL" / "VERTEX" / "FACE" / "EDGE" '
                mspace = self.input.usetc(prompt, 'CELL')
                mspace = mspace.upper()
                self.meshspace.append(mspace)

                for meshspace in self.mesh.mslist:
                    if mspace.upper() == meshspace.name:
                        tmp           = int(meshspace.size)
                        ninds         = str(tmp)
                        nindsh        = str(int(0.5*tmp))
                try:
                    x = ninds
                except:
                    print >> self.fpwatch
                    print >> self.fpwatch, '%s meshspace not in the mesh' %(mspace)
                    print >> self.fpwatch, 'Will not define the region based on this selection'
                    print >> self.fpwatch
                    return 0
            
                prompt = '    IDs: (e.g.: "1-' + ninds +',' + nindsh + ',' + ninds + ')'
                self.selectorCriteria = 'IDs: ' + self.input.usetc(prompt, nindsh)

                prompt = 'Region is '+ mspace+'(s) with  '+ self.selectorCriteria +'.'
                accept = self.input.usetc(prompt+' Accept? (y/n)','y')
                accept.lower()

                if accept == 'n' or accept=='no':
                    #region is not wanted so discard it
                    print >> self.fpwatch
                    print >> self.fpwatch, 'Region is not acceptable so discard it'
                    print >> self.fpwatch
                    return 0

                string                   = self.selectorCriteria
                indices                  = self.__getIDIndices(string,ninds)

                self._IDMBSindices[mspace][cycleid].extend(indices)

                self._IDMBSindices[mspace][cycleid]     = unique(self._IDMBSindices[mspace][cycleid])

            if selection == 'MB':

                #provide indices based on mesh block
                if self.mesh.cells.has_key('BLOCKID'):
                
                    blockdata             = self.mesh.cells['BLOCKID']
                    def_mblock            = blockdata[0]
                    val_mblocks           = str(unique(blockdata))
                    prompt                = 'Existing mesh blocks are %s. Provide mesh block number :' %(val_mblocks)
                    mblock                = self.input.useti(prompt,def_mblock)
                    mspace                = 'CELL'
                    self.meshspace.append(mspace)
                    self.selectorCriteria = 'mesh block number:' + str(mblock)

                    indices = self.__getMBIndices(mblock)

                    self._IDMBSindices[mspace][cycleid].extend(indices)
                    
                    self._IDMBSindices[mspace][cycleid]     = unique(self._IDMBSindices[mspace][cycleid])
                    
                else:
                    print >> self.fpwatch
                    print >> self.fpwatch, 'No mesh blocks specified in output file'
                    print >> self.fpwatch, 'This region will not contain indices based on mesh blocks'
                    print >> self.fpwatch
                    return 0
            
            if selection == 'S':
                    
                #provide indices based on spatial domain input

                format                = '%10.2E'
                thecelldata           = self.mesh.cells['CENTROIDS']
                thevertexdata         = self.mesh.vertices['COORDS']
                min_extents           = []
                min_extents.append(thecelldata[Numeric.argmin(thecelldata[:,0]),0])
                min_extents.append(thecelldata[Numeric.argmin(thecelldata[:,1]),1])
                min_extents.append(thecelldata[Numeric.argmin(thecelldata[:,2]),2])
                max_extents           = []
                max_extents.append(thecelldata[Numeric.argmax(thecelldata[:,0]),0])
                max_extents.append(thecelldata[Numeric.argmax(thecelldata[:,1]),1])
                max_extents.append(thecelldata[Numeric.argmax(thecelldata[:,2]),2])
                def_extents           = str(min_extents)+'->'+str(max_extents)
                posns                 = self.input.usetc('Provide spatial extents: (i.e [0,0,0]->[1,1,1])',def_extents)
                self.selectorCriteria = 'spatial domain: ' + posns

                posns = posns.replace('->',',')
                posns = posns.replace('[','')
                posns = posns.replace(']','')
                posns = posns.replace(' ','')
                posns = posns.split(',')

                posnspecerr  = "\n \n      For S criteria, spatial box specification needs to be of the form '[x1,y1,z1->x2,y2,z2]' \n \n"                 
                assert len(posns) == 6, posnspecerr

                mspaces          = []
                for ms in self.mesh.mslist:
                    mspaces.append(ms.name)
                self.meshspace.append(mspaces)

                prompt                = 'Region is '+ str(mspaces)+'(s) with  '+ self.selectorCriteria +'.'
                accept                = self.input.usetc(prompt+' Accept? (y/n)','y')
                accept.lower()
                if accept == 'n' or accept=='no':
                    #region is not wanted so discard it
                    print >> self.fpwatch
                    print >> self.fpwatch, 'Region is not acceptable so discard it'
                    print >> self.fpwatch
                    return 0

                for meshspace in mspaces:

                    indices = self.__getSIndices(meshspace,posns)
                    self._IDMBSindices[meshspace][cycleid].extend(indices)

                    self._IDMBSindices[meshspace][cycleid] = unique(self._IDMBSindices[meshspace][cycleid])


            if selection == 'VAR':
                
                #provide indices based on variable value

                count       = 0
                for i in self.storage.tlist:
                    if cycleid == i.id:
                        theid = count
                    count = count + 1

                count       = 0
                #assemble variable list for the user to choose from 
                varlist     = '\n'
                for var in self.storage.tlist[theid].vlist:
                    s       = '%+25s' %(str(var.nickname))
                    s       = s + ','
                    varlist = varlist + s
                    count   = count + 1
                    if (count%self.varformat == 0):
                        varlist = varlist + '\n'
                varlist = varlist[0:len(varlist)-1]
                if count < 1 :
                    print >> self.fpwatch
                    print >> self.fpwatch, 'No variables live in this file'
                    print >> self.fpwatch, 'Will not define the region based on this selection'
                    print >> self.fpwatch
                    return 0
                else:
                    prompt  = 'Which variable? (%s)' %(varlist)
                    varname = self.input.usetc(prompt, var.nickname)
                
                try:
                    for var in self.storage.tlist[theid].vlist:
                        if varname == var.nickname:
                            thisvar = var
                    x = thisvar.rank
                except:
                    print >> self.fpwatch
                    print >> self.fpwatch, 'Invalid entry.  Must pick one of: ',varlist
                    print >> self.fpwatch, 'Will not define the region based on this selection'
                    print >> self.fpwatch
                    return 0

                self.meshspace.append(thisvar.meshspace)

                #now get data for this variable

                vlist            = [thisvar]
                self.storage.getValues(vlist)
                vdata            = thisvar.data

                #check dimensions of the variable to see if more user input is required
                dims             = 1
                desc1            = ''
                desc2            = ''

                dim              = 1
                face             = 1
                if thisvar.rank > 1:
                    dims         = thisvar.shape[len(thisvar.shape)-1]
                if dims < 2 :                     
                    vmin         = vdata[Numeric.argmin(vdata)]
                    vmax         = vdata[Numeric.argmax(vdata)]
                else:
                    #dealing with vectors, arrays so need additional user input
                    dimlist      = range(dims+1)
                    dimlist      = str(dimlist[1:len(dimlist)])
                    print >> self.fpwatch, 'This variable is a vector of ' + str(dims) + ' dimensions.'
                    dim          = self.input.usetc('Specify which component you would like to choose your region from ' + dimlist, '1')
                    dim          = int(dim)
                    desc1        = '_' + str(dim)
                    if thisvar.rank == 2:
                        vmin     = vdata[Numeric.argmin(vdata[:,dim-1]),dim-1]
                        vmax     = vdata[Numeric.argmax(vdata[:,dim-1]),dim-1]
                    if thisvar.rank == 3:
                        print >> self.fpwatch, 'This variable is a face-centred vector.'
                        facelist = range(thisvar.shape[len(thisvar.shape)-2]+1)
                        facelist = str(facelist[1:len(facelist)])
                        face     = self.input.usetc('Specify which face you would like to choose your region from ' + facelist, '1')
                        face     = int(face)
                        desc2    = '_' + str(face)
                        vmin     = vdata[Numeric.argmin(vdata[:,face-1,dim-1]),face-1,dim-1]
                        vmax     = vdata[Numeric.argmax(vdata[:,face-1,dim-1]),face-1,dim-1]

                prompt = 'Lower limit (inclusive) for %s '%(thisvar.name)
                vmin2  = '%10.3e' %(vmin)
                lower  = self.input.usetc(prompt, vmin2)
        
                prompt = 'Upper limit (inclusive) for %s '%(thisvar.name)
                vmax2  = '%10.3e' %(vmax) 
                upper  = self.input.usetc(prompt, vmax2)  

                self.selectorCriteria = thisvar.nickname + desc1 + desc2 + ' between ' + lower + ' and ' + upper
                prompt                = 'Region is '+ thisvar.meshspace+'(s) with  '+ self.selectorCriteria +'.'
                accept                = self.input.usetc(prompt+' Accept? (y/n)','y')
                accept.lower()
                if accept == 'n' or accept=='no':
                    #region is not wanted so discard it
                    print >> self.fpwatch
                    print >> self.fpwatch, 'Region is not acceptable so discard it'
                    print >> self.fpwatch
                    return 0
            
                self._var0indices = self.__getVARIndices(thisvar,cycleid,lower,upper,dim,face)

                #for the VAR selector store choice of variable,dim,face,lower,upper for
                #future indices storage when cycleid > 0

                self.var          = thisvar
                self.dim          = dim
                self.face         = face
                self.lower        = lower
                self.upper        = upper


        elif (cycleid > 0 and 'VAR' in self.selector):
            """
            For cases where the original region indices alter with time.
            These indices must be stored for the relevant cycleid.
            Perform this selection process only when indices for this
             cycleid do not exist.
            """
            if not self.indices[self.var.meshspace].has_key(cycleid):

                count     = 0
                for i in self.storage.tlist:
                    if cycleid == i.id:
                        theid = count
                    count = count + 1

                for var in self.storage.tlist[theid].vlist:
                    if self.var.nickname == var.nickname:
                        thisvar = var

                varindices = self.__getVARIndices(thisvar,cycleid,self.lower,self.upper,self.dim,self.face)

                for meshspace in self.mesh.mslist:
                    if meshspace.name == thisvar.meshspace:
                        self.indices[meshspace.name][cycleid] = varindices


                """
                Now concatentate indices resulting from the ID,S,MB
                 selectors at cycleid 0 with the indices resulting from
                 this VAR selector.
                For this cycleid, return a set of indices resulting
                 from all the region's selectors.
                """

                for ms in self._IDMBSindices:
                    if ms == thisvar.meshspace:
                        length = len(self._IDMBSindices[thisvar.meshspace][0])
                        L      = self._IDMBSindices[thisvar.meshspace][0][0:length]
                        L.extend(self.indices[thisvar.meshspace][cycleid])
                        L = unique(L)
                        L.sort()
                        self.indices[ms][cycleid] = L

        """
        If all choices in the selection process are valid return a success
        """
        return 1


    def __getIndicesFromDeveloper(self,criteria,cycleid):
        "get indices when criteria is directly from developer"

        self.abortstring = getAbortString()
        
        criterrstring    = "\n \n criteria is %s and must be of the " \
                           %(str(criteria))
        criterrstring   += "form [selector,*] where selector is either "
        criterrstring   += "'ID','VAR','S','MB', 'Default' \n \n" 
        criterrstring   += self.abortstring

        assert len(criteria) > 0, criterrstring
        assert criteria[0] in ['ID','VAR','S','MB','Default'], criterrstring

        meshspacelist = ['CELL','VERTEX']

        for meshspace in meshspacelist:
            if not self.indices.has_key(meshspace) :
                self.indices[meshspace] = {}
            if not self.indices[meshspace].has_key(cycleid) :
                self.indices[meshspace][cycleid] = []

        iderrstring    = "\n \n For ID criteria must choose either "
        iderrstring   += "'CELL' or 'VERTEX' \n \n"
        iderrstring   += self.abortstring

        varerrstring   = "\n \n For VAR criteria must choose a valid XML variable \n"
        varerrstring  += "      Please consult your BasicRun.log file for a "
        varerrstring  +=       "list of valid XML variables. \n"
        varerrstring  += self.abortstring

        varangerr      = "\n \n For VAR criteria must choose a field value range"
        varangerr     +=       "in the form [x,y] \n \n"
        varangerr     += self.abortstring

        if criteria[0] == 'ID':
            #for ID choice criteria is of the form [ID,CELL|VERTEX,x-y|x,y,z]
            
            assert criteria[1] in meshspacelist, iderrstring

            self.meshspace.append(criteria[1])
            
            for meshspace in self.mesh.mslist:
                if meshspace.name == 'CELL'  : size=int(meshspace.size)
                if meshspace.name == 'VERTEX': size=int(meshspace.size)

            self.indices[criteria[1]][cycleid] = \
                 self.__getIDIndices('IDs: ' + criteria[2],size) 
            
        elif criteria[0] == 'VAR':
            #for VAR choice criteria is of the form [VAR,variable,[a,b],dim,face]

            assert criteria[1].rank > 0, varerrstring
            assert len(criteria[2]) == 2, varangerr
            
            varspecerr     = "\n \n      For VAR criteria, variable %s rank > 1" \
                             %(criteria[1].name)
            varspecerr    +=      " but you did not specify   \n" 
            varspecerr    += "      which dimension and face the field value range"
            varspecerr    +=       "is applicable for.\n"
            varspecerr    += "      Your criteria should be of the form "
            varspecerr    +=       "[VAR, variable, [x,y], dim, face].\n"
            varspecerr    += self.abortstring

            vardimerr      = "\n \n      For VAR criteria, variable %s rank = 2" \
                             %(criteria[1].name)
            vardimerr     +=       "but you did not specify  \n" 
            vardimerr     += "      which dimension the field value range "
            vardimerr     +=       "is applicable for.\n"
            vardimerr     += "      Your criteria should be of the form "
            vardimerr     +=       "[VAR, variable, [x,y], dim].     \n"
            vardimerr     += self.abortstring

            varfacerr      = "\n \n      For VAR criteria, variable %s rank = 3" \
                              %(criteria[1].name)
            varfacerr     += "but you did not specify    \n"
            varfacerr     += "      which dimension and which face the field "
            varfacerr     +=        "value range is applicable for.\n"
            varfacerr     += "      Your criteria should be of the form "
            varfacerr     +=        "[VAR, variable, [x,y], dim, face]  \n"
            varfacerr     += self.abortstring

            if len(criteria) == 3:
                assert criteria[1].rank < 2, varspecerr
                dim  = None
                face = None
            if len(criteria) == 4:
                assert criteria[1].rank == 2, vardimerr
                dim  = int(criteria[3])
                face = None
            if len(criteria) == 5:
                assert criteria[1].rank == 3, varfacerr
                dim  = int(criteria[3])
                face = int(criteria[4])

            lower = float(criteria[2][0])
            upper = float(criteria[2][1])

            self.meshspace.append(criteria[1].meshspace)

            count = 0
            for i in self.storage.tlist:
                if cycleid == i.id:
                    theid = count
                count += 1
                
            for var in self.storage.tlist[theid].vlist:
                if criteria[1].nickname == var.nickname:
                    thisvar = var
            
            self.var = thisvar
            self.indices[criteria[1].meshspace][cycleid] = \
                 self.__getVARIndices(thisvar,cycleid,lower,upper,dim,face)

        elif criteria[0] == 'S':
            #for S choice, criteria is of the form [S,'[x1,y1,z1->x2,y2,z2]']

            posns        = criteria[1]
            posnspecerr  = "\n \n      For S criteria, spatial box specification"
            posnspecerr += "needs to be of the form '[x1,y1,z1->x2,y2,z2]' \n \n" 
            posnspecerr += self.abortstring

            assert '['  in posns, posnspecerr
            assert '->' in posns, posnspecerr

            posns = posns.replace('->',',')
            posns = posns.replace('[','')
            posns = posns.replace(']','')
            posns = posns.replace(' ','')
            posns = posns.split(',')

            assert len(posns) == 6, posnspecerr

            for thisMS in meshspacelist:
                self.meshspace.append(thisMS)
                self.indices[thisMS][cycleid] = self.__getSIndices(thisMS,posns)

        elif criteria[0] == 'MB':
            #for MB choice, criteria is of the form [MB, number]'


            meshblockerr  =  "\n \n      Chosen the MB criteria,"
            meshblockerr += " but mesh %s does not have any mesh blocks outputted.\n \n" %(self.mesh.name)
            meshblockerr += self.abortstring

            assert self.mesh.cells.has_key('BLOCKID'), meshblockerr

            blockdata     = self.mesh.cells['BLOCKID']
            def_mblock    = blockdata[0]
            val_mblocks   = str(unique(blockdata))

            meshblocknumerr  = "\n \n      For MB criteria, chosen MB = "
            meshblocknumerr += "%s but available mesh block numbers are %s.\n \n" %(str(criteria[1]),val_mblocks)
            meshblocknumerr += self.abortstring

            for mblock in criteria[1:len(criteria)-1]:
                assert str(mblock) in val_mblocks, meshblocknumerr

            self.meshspace.append('CELL')

            L = self.__getMBIndices(int(criteria[1]))
            for mblock in criteria[2:len(criteria)]:
                L.extend(self.__getMBIndices(int(mblock)))
            L = unique(L)
            L.sort()
            self.indices['CELL'][cycleid] = L
            
        elif criteria[0] == 'Default':

            for meshspace in self.mesh.mslist:
                self.meshspace.append(meshspace.name)

                if not self.indices.has_key(meshspace.name):
                    self.indices[meshspace.name]          = {}

                if not self.indices[meshspace.name].has_key(cycleid):
                    self.indices[meshspace.name][cycleid] = []
                                            
                if meshspace.size != None:
                    x = meshspace.size+1
                    L = range(x)
                    T = L[1:len(L)]
                    self.indices[meshspace.name][cycleid] = T

       #if all choices in the selection process are valid return a success
        return 1
        

    def __getIDIndices(self,input,ninds):

        indices               = []

        string                = input.replace(' ',',')
        while (string.find(',,') > -1):
            string            = string.replace(',,',',')
        if (string[0] == ','):
            string = string[1:]
        sl                    = string.split(',')

        for c in sl[1:len(sl)]:
            cl = c.split('-')
            if(len(cl)):
                i1 = int(cl[0])
                if(len(cl)>1):
                    i2 = int(cl[1])+1
                else:
                    i2 = i1+1
                indices = indices + range(i1,i2)
        indices.sort()

        newindices = []
        for i in range(len(indices)):
            if(0 < indices[i] <= int(ninds) and indices[i]>0):
                newindices.append(indices[i])
        
        return newindices


    def __getVARIndices(self,thisvar,cycleid,lower,upper,dim=None,face=None):

        indices          = []

        vlist            = [thisvar]
        self.storage.getValues(vlist)
        v                = thisvar.data

        for i in range(len(v)):
            
            if thisvar.rank == 1:
                if(v[i]>=float(lower) and v[i]<=float(upper)):
                    indices.append(i+1)
            elif thisvar.rank == 2:
                if(v[i][dim-1]>=float(lower) and v[i][dim-1]<=float(upper)):
                    indices.append(i+1)
            else:
                if(v[i][face-1][dim-1]>=float(lower) and v[i][face-1][dim-1]<=float(upper)):
                    indices.append(i+1)

        
        return indices

    def __getMBIndices(self,mblock):

        indices   = []

        blockdata = self.mesh.cells['BLOCKID']
        
        for i in range(len(blockdata)):
            if blockdata[i] == mblock:
                indices.append(i+1)

        return indices


    def __getSIndices(self,mspace,posns):

        indices = []

        if mspace == 'CELL':
            thedata = self.mesh.cells['CENTROIDS']
        elif mspace == 'VERTEX':
            thedata = self.mesh.vertices['COORDS']
        else:
            return indices
            
        xlength = max(thedata[:,0]) - min(thedata[:,0])
        ylength = max(thedata[:,1]) - min(thedata[:,1])
        zlength = max(thedata[:,2]) - min(thedata[:,2])
        length  = max(xlength,ylength,zlength)
        
        tol     = 0.0001*length/(len(thedata)+1)
        
        for i in range(len(thedata)):
            
            xposn  = float(posns[0]) - tol
            xposn2 = float(posns[3]) + tol
            if (xposn <= thedata[i,0] <= xposn2):

                yposn  = float(posns[1]) - tol
                yposn2 = float(posns[4]) + tol
                if (yposn <= thedata[i,1] <= yposn2):

                    zposn  = float(posns[2]) - tol
                    zposn2 = float(posns[5]) + tol
                    if (zposn <= thedata[i,2] <= zposn2):

                        indices.append(i+1)

        return indices
    
    def str(self):

        print >> self.fpwatch, self.name
        print >> self.fpwatch, self.file
        print >> self.fpwatch, self.mesh.name
        print >> self.fpwatch, self.meshspace
        print >> self.fpwatch, self.selector
        print >> self.fpwatch, self.describe
        print >> self.fpwatch, self.indices

if __name__=='__main__':

    import os, sys
    thisdir   = os.path.dirname(__file__)
    parsedir  = os.path.abspath(thisdir+'../') 
    sys.path.append(parsedir)
    from PYTHONutils import uTestOpts

    dir  = '../../../../tools/scripts/test_TBrookParse/samples/'
    file = 'map-b_output/map-b.TBrook.xml'

    opts = uTestOpts('fd', defaults={'f': dir+file})
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fpwatch = sys.stdout
    try:

        import usubs
        inBuffer     = ''
        #input        = usubs.input(inBuffer)
        input        = None
                
        from XMLgetStorageObject import getStorageObject
        tstorage     = getStorageObject(opt.f)
        tstorage.str()
        
        from PYTHONutils import unique, getAbortString

        "Three tries..."
        
        newregion = getRegionObject(input, file=opt.f,
                                    storage = tstorage, mesh=tstorage.mlist[0],
                                    selector= ['ID','VERTEX','50-80'],
                                    fpwatch = fpwatch,
                                    debug   = opt.d)
        newregion.getIndices(newregion.selector,cycleid=0)
        newregion.str()

        newregion = getRegionObject(input, file=opt.f,
                                    storage = tstorage, mesh=tstorage.mlist[0],
                                    selector= ['S','[0.045,-0.05,0.0->0.05,-0.01,0.03]'],
                                    fpwatch = fpwatch, 
                                    debug   = opt.d)
        newregion.getIndices(newregion.selector,cycleid=0)
        newregion.str()

        newregion = getRegionObject(input, file=opt.f,
                                    storage = tstorage, mesh=tstorage.mlist[0],
                                    selector= ['MB',4,5],
                                    fpwatch = fpwatch,
                                    debug   = opt.d)
        newregion.getIndices(newregion.selector,cycleid=0)
        newregion.str()
        
    except:
        print >> fpwatch, "---> Test failed in some aspect <---"
        print >> fpwatch,  "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise


    







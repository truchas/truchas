"""

 XMLgetVariable

 -----------------------------------------------------------------------------
  Purpose:
  
     To instantiate and populate a getVariable object.
  
  Public Interface(s):
    v = getVariable(srcPointer,element)
    v.getData()
    v.str()

  Contains:
    class getVariable
        __init__(self,srcPointer,element,parse=MiniDom,debug=0)
        getData(self)
        str(self)

    Unit Test Block
  
  Author(s): Sharen Cummins (scummins@lanl.gov)
             Erin Iesulauro Barker (eibarker@lanl.gov)
 -----------------------------------------------------------------------------
"""

if __name__=='__main__':
    # Set sys.path for component testing mode
    import os, sys
    thisdir    = os.path.dirname(__file__)
    parserdir  = os.path.abspath(thisdir+'../')
    sys.path.append(parserdir)

from BASEutils      import baseVariable
from PYTHONutils    import SystemExit,getpath
from FORTRANutils   import fortransupport
from XMLutils       import MiniDom,BrookUtilities,BrookVariable
import string

class getVariable(baseVariable):

    def __init__(self,srcPointer,element,parse=MiniDom,debug=0):

        #srcpointer is xml pointer

        "initialise the variable object"
	baseVariable.__init__(self)

	self.__parse      = parse
        self.__pointer    = {}   # all the attributes particular to the XML parser

        "now create the variable"
        self.mesh      = parse().getElemAttribute(element,'Mesh')
        self.meshspace = parse().getElemAttribute(element,'Map')
        if len(self.meshspace) == 0:
            self.meshspace = 'CELL' # default

        spacelu        = srcPointer['spacelu']
        if spacelu.has_key(self.meshspace):
            self.meshspace = spacelu[self.meshspace]

        datafmt        = srcPointer['datafmt']
        tfile          = srcPointer['tfile']        
        
        name,datat,rank,shape = BrookVariable(datafmt,tfile,element,parse).\
                                getMETAVALUES() #don't get data
        self.name      = name
        self.nickname  = ''  #lookup a dictionary for this given name
        self.rank      = rank
        self.shape     = shape
        self.type      = datat
        
        if 'PHASECONC' in name:
            if '_' in name:
                #probe phase concentration
                tmp = string.split(name,'_')
                for j in tmp[1]:
                    w = string.find(j,'0')
                    if w < 1:
                        tmp[1] = tmp[1][w+1:len(tmp[1])]
                s = tmp[1]
                for j in tmp[2]:
                    w = string.find(j,'0')
                    if w < 1:
                        tmp[2] = tmp[2][w+1:len(tmp[2])]
                n = tmp[2]
            else:
                #field variable concentration
                s = parse().getElemAttribute(element,'s')
                n = parse().getElemAttribute(element,'n')
                s = s[len(s)-1:len(s)]
                n = n[len(n)-1:len(n)]
            if 'OLD' in name:
                name = 'PCO'
            else:
                name = 'PC'
            self.name = name + '%s_%s' %(s,n)
        

        """
        additional items stored for convenience
        - these are particular to the XML parser
        """
        self.__pointer['node']    = element
        self.__pointer['datafmt'] = srcPointer['datafmt']
        self.__pointer['tfile']   = srcPointer['tfile']
        self.__pointer['offset']  = int(parse().getElemAttribute(element,'Offset'))

    def getData(self):

        dformat = self.__pointer['datafmt']
        tfile   = self.__pointer['tfile']
        offset  = self.__pointer['offset']
        element = self.__pointer['node']

        probe   = (self.mesh == 'NoMesh')
        
        if self.data == None:
            if probe:
                # probe variable -
                # data and shape are obtained by closer inspection of the binary file
                self.shape,self.data = BrookVariable(dformat,tfile,element,self.__parse).\
                                       getPROBEVALUES(self.name,self.type,self.shape,offset)
            else:
                self.data = BrookVariable(dformat,tfile,element,self.__parse).\
                            getVALUES(self.name,self.type,self.shape,offset)
                

    def metastr(self):

        print self.name, self.nickname, self.rank, self.shape,\
              self.type, self.mesh, self.meshspace, self.__pointer['offset']
        #print self.data


if __name__=='__main__':

    'for testing this component'

    from PYTHONutils import uTestOpts
    opts = uTestOpts('fd')
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)
    
    try:
        parse      = MiniDom
        node       = parse().getDoc(opt.f)
        nodes      = parse().getElements(node,'MESH')
        for node in nodes:
            filenode     = parse().getElement(node,'FILE')
            meshfile     = str(parse().getValue(filenode))
            corrmeshfile = getpath(opt.f,meshfile)
        node     = parse().getDoc(corrmeshfile)
        meshnode = parse().getElement(node,'MESH')
        filenode = parse().getElement(meshnode,'FILE')
        file     = str(parse().getValue(filenode))
        corrfile = getpath(opt.f,file)
        file     = BrookUtilities().fixFileName(corrfile)
        tfile    = fortransupport.T_Binary_File(file)
        
        spacelu = {None:'CELL','NONE':'NONE', 'CELL':'CELL',
                   'CELLS':'CELL', 'FACES':'FACE', 'FACE':'FACE',
                   'VERTICES':'VERTEX', 'VERTEX':'VERTEX',
                   'NODES':'VERTEX','NODE':'VERTEX',
                   'EDGES':'EDGE','EDGE':'EDGE', 'MESH':'MESH'}

        pointer    = {'node':meshnode, 'datafmt':'binary',
                      'tfile':tfile, 'spacelu':spacelu}

        elems      = parse().getElements(meshnode,'FILEVAR')
        for element in elems:
            var  = getVariable(pointer,element,parse)
            var.metastr()

    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise

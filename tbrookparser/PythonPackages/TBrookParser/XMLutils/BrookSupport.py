

"""

 BrookSupport

 -----------------------------------------------------------------------------
  Purpose:
  
     The Great One only knows...

  Public Interface(s):
  
    Usage for BrookVariable
     mv = BrookVariable(datafmt,tfile,element).getMETAVALUES()
     v  = BrookVariable(datafmt,tfile,element).getVALUES(name,type,shape,offset)     
     pv = BrookVariable(datafmt,tfile,element).getPROBEVALUES(name,type,shape,offset)

   Usage for BrookUtilities
     file  = BrookUtilities().fixFileName(corrfile)
     
  Contains:
    class BrookVariable
        __init__(self,dformat,tfile,element,parse=MiniDom)
        getMETAVALUES(self)
        getVALUES(self,name,type,shape,offset)
        getPROBEVALUES(self,name,type,shape,offset)

    class BrookUtilities
        __init__(self)
        fixFileName(self, fileName, debug=0)

    ---> NO Unit Test Block
  
  Author(s): Sriram Swaminarayan (sriram@lanl.gov)
 -----------------------------------------------------------------------------
"""

from PYTHONutils import SystemExit

try:
    try:
       import numpy.oldnumeric as Numeric
    except ImportError:
       import Numeric
    except:
       raise
except:
    print "Numeric need to parse XML data, please install Python Numeric before continuing"
    SystemExit()
    
from Parsing import MiniDom

class BrookVariable:

    def __init__(self,dformat,tfile,element,parse=MiniDom):

        self.name     = None
        self.rank     = None
        self.shape    = None
        self.type     = None
        self.offset   = None
        self.mesh     = None
        self.map      = None
        self.arrorder = None
        self.data     = None
        self.element  = element
        self.parse    = parse
        self.dformat  = dformat
        self.tfile    = tfile

    def getMETAVALUES(self):

        "return all values describing the variable except for the actual data values" 

        if self.name == None:

            element     = self.element
            Parse       = self.parse
            self.name   = Parse().getElemAttribute(element,'Name')
            self.rank   = int(Parse().getElemAttribute(element,'Rank'))
            shape       = Parse().getElemAttribute(element,'Shape')
            if self.rank == 0:
                shape   = "1"
            self.shape  = [int(d) for d in shape.split()]
            self.type   = str(Parse().getElemAttribute(element,'DataType'))
            ptype       = Parse().typeConverter(self.type)
            self.offset = int(Parse().getElemAttribute(element,'Offset'))
            self.mesh   = Parse().getElemAttribute(element,'Mesh')
            self.map    = Parse().getElemAttribute(element,'Map')
            self.arrorder = Parse().getElemAttribute(element,'ArrayOrder')
            if self.arrorder != 'C':
                self.shape.reverse()
            assert self.dformat == 'binary'

        return self.name,self.type,self.rank,self.shape


    def getVALUES(self,name,type,shape,offset):

        file   = self.tfile
        data   = file.get(name,type,shape,offset)
        if self.arrorder != 'C':
            #data   = Numeric.reshape(data,shape)
	    try:
              data   = Numeric.reshape(data,shape)
	    except ValueError:
	      try:
		data = Numeric.reshape(data[0],shape)
	      except:
		raise ValueError('Attempted to correct Numeric.reshape error and failed')
        return data

    def getPROBEVALUES(self,name,type,shape,offset):

        file   = self.tfile
        info   = file.getInfo(name,type,shape,offset)

        #this loop is performed to check for inconsistent binary files due to Truchas
        #segmentation faults
        count  = 1 
        for i in range(len(info)-1):
            if info[i+1] != info[i]:
                count -= 1
                break
            count +=1
                
        L         = [count] 
        L.append(shape[0])
        newshape  = L
        data   = file.get(name,type,newshape,offset)
        if self.arrorder != 'C' :
            data   = Numeric.reshape(data,newshape)

        return newshape,data

class BrookUtilities:

    def __init__(self):
        
        self.debug = 0

    def fixFileName(self, fileName, debug=0):
        
        if(self.debug): print 'file in fixFileName is: ...', fileName,'...'
        p  = fileName.split('/')   # split the filename
        #p2 = fileName.split('/')
        #t  = p2.pop()
        t = p.pop()                   # Pop off the filename
        if(self.debug + debug): print '   popped file is: ', t, fileName
        file = ''
        for d in p:
            file = file + d + '/'
        file = file + t
        if(self.debug): print '  Fixed file is: ', file
        return(file)


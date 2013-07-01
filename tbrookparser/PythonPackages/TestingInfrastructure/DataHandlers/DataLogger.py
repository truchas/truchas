#!/usr/bin/env python
"""
 DataLogger

-----------------------------------------------------------------------------
   Purpose:
  
      Provides convenient methods for the developer to record data
      to a log file resulting from their TestCase. 
  
   Public Interface:
  
      T = DataLogger(filename,mode,buffering)
      T.writeField(field,region,desc)
      T.writeScalar(scalar,desc)
      T.write(string)
      T.close()
        
   Contains:
      class DataLogger
        __init__(filename,mode,buffering)
        writeField(field,region,desc)
        writeScalar(scalar,desc)
        write()
        close()
        __logThisField(field,face,hdr)
        __logThisFieldAndRegion(field,indices,posns,face,hdr)
        __rjustln(column)
        
   Unit Test Block
  
   Author: Sharen Cummins (scummins@lanl.gov)
           Erin Iesulauro Barker (eibarker@lanl.gov)
-----------------------------------------------------------------------------
"""

import os, sys
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

from writeStorageObject import writeStorageObject

class DataLogger(writeStorageObject):
    "class that allows developers to log fields in their TestCase"

    def __init__(self,filename='Capabability.log',mode='w',buffering=0,debug=0):

	self.filename    = filename
	self.logger      = open(self.filename,mode,buffering) 
        self.iformat     = '%8d'
        self.dformat     = '%15.3e'
        self.column      = 5
        self.column1     = 20
        self.fpwatch     = self.logger
        self.module_name = 'asciiwriter'
        self.debug       = debug

    def writeField(self,field,region=None,desc=None):

        self.logger.write('\n')
        st   = 'LOGGING FIELD  : %s' %(field.nickname)
        self.__rjustln(self.column1,st)
        if desc != None:
            st  = 'MY DESCRIPTION : ' + str(desc)
            self.__rjustln(self.column1,st)
            self.logger.write('\n')

        s   = '#' + '-' * 100
        self.__rjustln(self.column, s)

        s   = '#' + '%+15s %+25s %+25s' %('VAR', 'TIME', 'CYCLE')
        self.__rjustln(self.column, s)
        s   = '#' + '%+15s %25.5e %23d '%(field.nickname, field.time, field.cycle)
        self.__rjustln(self.column, s)
        s   = '#' + '-' * 100
        self.__rjustln(self.column, s)

        if field.rank == 3:
            #log values for all faces:
            facelist = range(field.shape[len(field.shape)-2]+1)
            facelist = facelist[1:]
        else:
            facelist = [0]

        if region == None:
            #no region specified...just write the field array out
            for face in facelist:
                if face > 0:
                    s   = '#' + '-' * 100
                    self.__rjustln(self.column, s)
                    s   = '%s %+35s %d'%('#','Value(s) For Face: ',face)
                    self.__rjustln(self.column, s)
                    s   = '#' + '-' * 100
                    self.__rjustln(self.column, s)
                else:
                    s = ''
                    self.__rjustln(self.column, s)
                self.__logThisField(field,face)
                
        else:
            #region specified ... write indices, positions, field array
            indices = region.indices[field.meshspace][field.stepid]

            #find co-ordinates for the variables
            if field.meshspace == 'CELL':
                #if the field lives on the cell space use CELL CENTROIDS for positions
                posns   = region.mesh.cells['CENTROIDS']
            if field.meshspace == 'VERTEX':
                #if the field lives on the vertex space use VERTEX COORDS for positions
                posns   = region.mesh.vertices['COORDS']


            for face in facelist:
                if face == 0:
                    s  = '#' + '-' * 100
                    self.__rjustln(self.column, s)
                    s  = '%s%s %+25s %+32s'%('#', 'Index', 'Position',  'Value(s)')
                    self.__rjustln(self.column, s)
                    s  = '#' + '-' * 100
                    self.__rjustln(self.column, s)
                else:
                    s  = '#' + '-' * 100
                    self.__rjustln(self.column, s)
                    s  = '%s%s %+25s %+50s %d'%('#', 'Index', 'Position',  'Value(s) at Face: ',face)
                    self.__rjustln(self.column, s)
                    s  = '#' + '-' * 100
                    self.__rjustln(self.column, s)
                self.logger.write('\n')
                self.__logThisFieldAndRegion(field,indices,posns,face)
            

    def writeScalar(self,scalar,desc=None):
	
        if desc != None:
            self.logger.write(desc + '\n')
	self.logger.write('Scalar to go here \n')

    def write(self,string):

	self.logger.write(string)

    def close(self):

	self.logger.close()


    def __logThisField(self,field,face,hdr=''):

        myshape    = Numeric.array(field.shape)
        nperline   = 1
        if len(myshape) > 1:
            nperline = myshape[1]
        if face == 0:
            myfield    = Numeric.array(field.data)
            myrank     = field.rank
        else:
            errstring  = "\n In %s face = %d but field rank = %d \n" %(self.__class__.__name__,face,field.rank)
            assert field.rank > 2, errstring
            myfield    = Numeric.array(field.data[:,face-1,:])
            myrank     = 2

        self.writeArray(myrank,myfield,hdr,self.dformat,self.logger,nperline=nperline)

    def __logThisFieldAndRegion(self,field,indices,posns,face,hdr=''):

        myindices   = Numeric.array(indices,'i')
        
        self.writeReducedFieldAndRegion(face,field,myindices,posns,hdr,self.logger)


    def __rjustln(self,column,s):

        tmp = s.rjust(column + len(s))
        self.logger.write(tmp)
        self.logger.write('\n')

if __name__== '__main__':
    " DataLogger component test"
    
    from DataCreator import getDataCreator
    from PYTHONutils import uTestOpts
    import Clean
    
    # file, debug, output prefix, binary, ascii, clean
    dfltdir  = '../TestCases/static_drop/static_drop_golden/'
    dfltfile = 'static_drop'
    opts = uTestOpts('fdc', 
                     defaults = {'d' : False,
                                 'c' : False,
				 'f' : dfltdir+dfltfile+'TBrook.xml'},
                     actions  = {'d' : 'store_true',
                                 'c' : 'store_true'},
                     dir      = dfltdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    try:

        if opt.c: #clean

	    print "Cleaning up from a past test run...\n"
	    Clean.Clean()

	else:
		
	    if opt.d: debug = 1
	    else:     debug = 0

            testdata  = getDataCreator(dir='../TestCases/static_drop/static_drop_golden')
            testdata.getStorage()
            newregion = testdata.getRegion(selector=['VAR','Density',[0.95,1.0]],desc='0.95<Density<1.0')
            newregion.getIndices(newregion.selector,cycleid=1)
            newregion.str()
            thefield = testdata.getField(field='Velocity',time=0.00105,region=newregion)
            logger   = DataLogger('Test.log',debug=debug)
            logger.writeField(field=thefield,region=newregion,desc='\n Component Test \n')
            logger.close()

    except:

	print "---> Test failed in some aspect <---"
	print "\nNature-Of-Error:",sys.exc_info()[0],"\n"
	if opt.d: raise
   

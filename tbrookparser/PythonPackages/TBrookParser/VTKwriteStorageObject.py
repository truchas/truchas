

"""

 VTKwriteStorageObject

 -----------------------------------------------------------------------------
  Purpose:

     Defines a subclass of writeStorageObject for the writing of VTK
     graphics files from Truchas simulation data
     
  Public Interface(s):

    VTKwSO = VTKwriteStorageObject(fileName, outputfmt, mesh, t, cycle, bVariables)

    VTKwSO.openFile(comment)
    VTKwSO.writeMesh()        
    VTKwSO.writeVariables()   

  Contains:
    class VTKwriteStorageObject(writeStorageObject)
        __init__(self, fileName, outputfmt, mesh, t, cycle, bVariables,
                 seq_no   = 0,
                 bFlags   = None,
                 bRegions = None,
                 dformat  = '%20.5e',
                 fpwatch  = sys.stdout,
                 debug    = 0)
        openFile(self, comment)
        writeMesh(self)
        writeVariables(self)   

    Unit Test Block

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
             Erin Iesulauro Barker (eibarker@lanl.gov)
 -----------------------------------------------------------------------------
"""

from writeStorageObject import writeStorageObject
import os, sys, re, struct
try:
   import numpy.oldnumeric as Numeric
except ImportError:
   import Numeric
except:
   raise

class VTKwriteStorageObject(writeStorageObject):

    def __init__(self, fileName, outputfmt, mesh, t, cycle, bVariables,
                 seq_no   = 0,
                 bFlags   = None,
                 bRegions = None,
                 dformat  = '%20.5e',
                 fpwatch  = sys.stdout,
                 debug    = 0):

        writeStorageObject.__init__(self,fpwatch,debug)
        self.debug    = debug
        self.fileName = fileName
        self.vars     = bVariables
        self.t        = t
        self.cycle    = cycle
        self.mesh     = mesh
        self.bFlags   = bFlags
        self.bRegions = bRegions
        self.seq_no   = seq_no
        self.mapLookup  = {None:0, 'CELL':0, 'NODE':1, 'VERTEX':1, 'FACE':2}
        self.outputfmt  = outputfmt

        self.fpwatch    = fpwatch
        self.dformat    = dformat

        self.typeLookup = {'i':'d', 'd':'f', 'l':'d', 'c':'s'}
        
        if(outputfmt == 'ascii'):
            self.binary = 0
	    self.asciiWrite()
        else:
            self.binary = 1
	    self.binaryWrite()

	#print >> self.fpwatch, 'Writing ', outputfmt,self.fileName,
        #self.FP = self.openFile(comment='cycle %d, t=%f'%(cycle,t))
        #self.writeMesh()
        #self.writeVariables()
        #self.FP.close()
	#print >> self.fpwatch, ' done'

        return

    def asciiWrite(self):
        "Overwrite the ascii write function in writeStorageObject"
	sys.stdout.flush()
	if self.debug: print >> self.fpwatch, 'Using asciiWrite'

	print >> self.fpwatch, 'Writing ascii ' + self.fileName

	# Open file
	self.FP = open(self.fileName,'w')

	# Write the header
	self.writeHeader(comment='cycle %d, t=%f'%(self.cycle,self.t))

	# Write the mesh
	self.writeMesh(endline='\n')

	# Write the variables
	self.writeVariables()

	# Close the file
	self.finalizeAndClose()

	print >> self.fpwatch, ' done'
	return

    def binaryWrite(self):
        "Overwrite the binary write function in writeStorageObject"

	sys.stdout.flush()
	if (self.debug): print >> self.fpwatch, 'Using binaryWrite'

	print >> self.fpwatch, 'Writing binary' + self.fileName

	# Open file
	self.FP = open(self.fileName,'w')

	# Write the header
	self.writeHeader(comment='cycle %d, t=%f'%(self.cycle,self.t))

	# Write the mesh
	self.writeMesh(endline='')

	# Write the variables
	self.writeVariables()

	# Close the file
	self.finalizeAndClose()

	print >> self.fpwatch, ' done'
	return
	    
	    
    def writeHeader(self, comment):
	"Write the VTK header information"

	if (self.debug): print >> self.fpwatch, '  Writing header'
	fp = self.FP
        fp.write('# vtk DataFile Version 1.0\n')
        fp.write(comment+'\n')
        
        if(self.binary):
            fp.write('BINARY\n')
        else:
            fp.write('ASCII\n\n')

	fp = None
        return 

    def writeMesh(self,endline):

        if(self.debug): print >> self.fpwatch, '  Initializing fp and mesh variables'
	
        fp = self.FP
        m = self.mesh
        if(self.debug): print >> self.fpwatch, '    writing header'
        fp.write('DATASET UNSTRUCTURED_GRID\n')
        if(self.debug): print >> self.fpwatch, '    writing vertices'
        v = self.mesh.vertices
        fp.write('POINTS %d double\n'%(v['N']))
        if(self.binary):
            #fp.write(v['COORDS'].tostring()+'\n')
	    for num in range(len(v['COORDS'])):
		for i in range(0,3):
	            data = struct.pack('!d',v['COORDS'][num][i])
	            fp.write(data)
        else:
            self.writeArray(2,v['COORDS'],'',self.dformat,fp,nperline=3)

        if(self.debug): print >> self.fpwatch, '    writing cells'
        # VTK uses 0-(n-1) numbering
        c = self.mesh.cells
        ncells = c['N']
        cc = Numeric.array(c['VERTICES'],'i')
        cc = cc-1
        fp.write('\nCELLS %d %d\n'%(ncells,ncells*9))
        if(self.binary):
	    for c in range(ncells):
		data = struct.pack('!i',8)
		for item in cc[c]: 
		    data += struct.pack('!i',item)
		fp.write(data)
        else:
	    # nperline really 1, but subtract 1 since the header is not empty
            self.writeArray(2,cc,'8 ','%d ',fp,nperline=8)

        fp.write('CELL_TYPES %d\n'%ncells)
        tmpArray = Numeric.zeros((ncells,1),'i')
        if ( self.mesh.type == 'HEX'):
            tmpArray = tmpArray + 12
        elif (self.mesh.type == 'TET'):
            tmpArray = tmpArray + 10

        if(self.binary):
	    for item in tmpArray:
		data = struct.pack('!i',item[0])
		fp.write(data)
        else:
            self.writeArray(2,tmpArray,'','%d'+endline,fp,nperline=10)

        tmpArrray = None
        fp = None
        return 

    def writeVariables(self):

        if(self.debug): print >> self.fpwatch, '  Writing out variables'
        fp = self.FP
        ncf = 0
        nvf = 0
	Cnamelist = []
	Vnamelist = []
	
        for v in self.vars:
            nr = v.rank
            if(nr == 2 ): nr = v.shape[nr-1]
	    writevar = self.writeThisVar(v)
                
            if (writevar) :
                if (v.meshspace == 'CELL' and v.nickname not in Cnamelist) :
                    ncf = ncf + 1 
		    Cnamelist.append(v.nickname)
                elif (v.meshspace == 'VERTEX' and v.nickname not in Vnamelist) :
                    nvf = nvf + 1 
		    Vnamelist.append(v.nickname)

	namelist = []
        if(ncf > 0 ):
            fp.write('\n')
            fp.write('CELL_DATA %d\n'%self.mesh.cells['N'])
            fp.write('FIELD %d %d\n'%(self.cycle,ncf))
            for v in self.vars:
	        writevar = self.writeThisVar(v)
		if (v.nickname not in namelist and writevar):
		    namelist.append(v.nickname)
                    nr = v.rank
                    if(nr > 0 and v.meshspace == 'CELL'):
		        if (nr == 1):
                            fp.write('%s %d %d double\n'%(v.nickname,nr,len(v.data)))
		        else:
                            fp.write('%s %d %d double\n'%(v.nickname,v.shape[nr-1],len(v.data)))
                        if(self.binary):
	                    for item in v.data:
				data = ''
				if nr > 1:
				    for n in item:
	                	        data += struct.pack('!d',n)
				else:
	                	    data = struct.pack('!d',item)
	                	fp.write(data)

                        else:
                            self.writeArray(nr,v.data,'',self.dformat,fp,nperline=1)

	namelist = []
        if(nvf > 0 ):
            fp.write('\n')
            fp.write('POINT_DATA %d\n'%self.mesh.vertices['N'])
            fp.write('FIELD %d %d\n'%(self.cycle,nvf))
            for v in self.vars:
	        writevar = self.writeThisVar(v)
		if (v.nickname not in namelist and writeVar):
                    if(v.rank > 0 and v.meshspace == 'VERTEX'):
                        fp.write('%s %d %d double\n'%(v.nickname,v.rank,len(v.data)))
                        if(self.binary):
	                    for item in v.data:
				if v.rank > 1:
				    data = ''
				    for n in item:
	                	        data += struct.pack('!d',n)
				else:
	                	    data = struct.pack('!d',item)
	                	fp.write(data)
                        else:
                            self.writeArray(1,v.data,'',self.dformat,fp,nperline=1)
                        fp.write('\n')
            
if __name__ =='__main__':
   
    'for component testing VTK write component'

    dfltdir = '../../scripts/test_TBrookParse/samples/'
    
    from PYTHONutils import uTestOpts

    # file, debug, output prefix, binary, ascii, clean
    opts = uTestOpts('fdoc', 
                     defaults = {'c' : False,
                                 'o' : 'VTKwriteUnitTest'},
                     actions  = {'c' : 'store_true'},
                     dir      = dfltdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fpwatch    = sys.stdout
    try:

        if opt.c:
            print >> fpwatch, "Cleaning up from a past test run...",
            cmd = 'rm -f %s' %(dfltdir+opt.o+'.binary.vtk')
            os.system(cmd)
            cmd = 'rm -f %s' %(dfltdir+opt.o+'.ascii.vtk')
            os.system(cmd)
            print >> fpwatch, "Cleanup successful!!"

        else:
            from XMLgetStorageObject import getStorageObject
            storage   = getStorageObject(opt.f)	
            mesh      = storage.mlist[0] #assume default mesh
            
            #now get the data values for the mesh fields
            for meshspace in mesh.mslist:
                meshspace.vlist = storage.getValues(meshspace.vlist)

            # fill in actual data values
            mesh.fillMesh()
            
            itimeout  = 1
            t         = storage.tlist[itimeout].time
            vars      = storage.tlist[itimeout].vlist
            vars      = storage.getValues(vars)

	    # ASCII test
            print >> fpwatch, '*'*60
            print >> fpwatch, 'Testing ascii write...'
            fileout = dfltdir + opt.o + '.ascii.vtk'
            writer  = VTKwriteStorageObject(fileout, 'ascii', mesh, t,
                                            itimeout, vars, fpwatch=fpwatch, 
					    debug=opt.d)
            print >> fpwatch, 'VTK file %s\n  written to %s' % \
			      (opt.o+'.ascii.vtk',dfltdir)

	    # Binary test
            print >> fpwatch, '*'*60
            print >> fpwatch, 'Testing binary write...'
            fileout = dfltdir + opt.o + '.binary.vtk'
            writer  = VTKwriteStorageObject(fileout, 'binary', mesh, t,
                                            itimeout, vars, fpwatch=fpwatch, 
					    debug=opt.d)
            print >> fpwatch, 'VTK file %s\n  written to %s' % \
			      (opt.o+'.binary.vtk',dfltdir)
                
            print >> fpwatch, '*'*60
            print >> fpwatch, \
	         '--> Note: to clean up, re-execute same command line plus "-c"'

    except:
        print >> fpwatch, "---> Test failed in some aspect <---"
        print >> fpwatch, "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise


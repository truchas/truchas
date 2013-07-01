

"""

 ENSIGHTwriteStorageObject

 -----------------------------------------------------------------------------
  Purpose:

     Defines a subclass of writeStorageObject for the writing of ENSIGHT
     graphics files from Truchas simulation output.

  Public Interface(s):

    EwSO = ENSIGHTwriteStorageObject(fileName, outputfmt, mesh, t, cycle, bVariables)
    EwSO.writebCase()
    EwSO.writebMesh()
    EwSO.writebVariables()
    EwSO.writeMesh()
    EwSO.writeCase()
    EwSO.writeVariables()
    EwSO.finalizeAndClose()
    EwSO.writeByComponent(data, rank, fp)

  Contains:
    class ENSIGHTwriteStorageObject(writeStorageObject)
        __init__(self, fileName, outputfmt, mesh, t, cycle, bVariables,
                 seq_no   = 0,
                 bFlags   = None,
                 bRegions = None,
                 dformat  = '%20.5e',
                 fpwatch  = sys.stdout,
                 debug    = 0)
        writebCase(self)
        writebMesh(self)
        writebVariables(self)
        writeMesh(self)
        writeCase(self)
        writeVariables(self)
        finalizeAndClose(self) 
        writeByComponent(self, data, rank, fp,
                         transpose=0,ascii=0,format='%12.5f',nperline=10)

    Unit Test Block

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""

from writeStorageObject import writeStorageObject
from ENSIGHTutils       import getBlockIndicies
from FORTRANutils       import fortransupport
from PYTHONutils        import unique,getdir
import os, sys, string
try:
   import numpy.oldnumeric as Numeric
except ImportError:
   import Numeric
except:
   raise

class ENSIGHTwriteStorageObject(writeStorageObject):

    def __init__(self, fileName, outputfmt, mesh, t, cycle, bVariables,
                 seq_no   = 0,
                 bFlags   = None,
                 bRegions = None,
                 dformat  = '%20.5e',
                 fpwatch  = sys.stdout,
                 debug    = 0):

        writeStorageObject.__init__(self,fpwatch,debug)
        self.fileName   = fileName
        self.vars       = bVariables
        self.t          = t
        self.mesh       = mesh
        self.bFlags     = bFlags
        self.bRegions   = bRegions
        self.debug      = debug
        self.strsize    = 80
        self.mesh_type  = 'hexa8'
        self.mapLookup  = {None:0, 'CELL':0, 'NODE':1, 'VERTEX':1, 'FACE':2}
        self.typeLookup = {'i':'d', 'd':'f', 'l':'d', 'c':'s'}

        ibeg            = string.rfind(fileName,"/") + 1
        self.simname    = fileName[ibeg:]

        self.dirname    = '%s.ensight'%(fileName)

        self.fpwatch    = fpwatch
        self.dformat    = dformat

        if (seq_no == 0):
            try:
               tmpdir = os.getcwd()
               os.chdir(self.dirname)
               os.chdir(tmpdir)
            except:
               if (self.debug): print >> self.fpwatch, "directory name issue"
               os.mkdir(self.dirname)
        self.FP       = None
        self.cycle    = cycle
        self.seq_no   = seq_no

        if (outputfmt == 'ascii'):
           self.asciiWrite()
        if (outputfmt == 'binary'):
           self.binaryWrite()
        
        return	

    #define all the (b)inary routines specific to ENSIGHT

    def writebCase(self):

        #  Write EnSight Gold Case file
        global time_val
        global file_val

        filename = '%s.ensight.CASE'%(self.fileName)
        fp       = open(filename, 'w')
        if self.seq_no == 0:
            time_val = []
            file_val = []
        time_val.append(self.t)
        file_val.append(self.cycle)

        #  Write header section
        print >> fp, 'FORMAT'
        print >> fp, 'type:                   ensight gold'
        print >> fp

        #  Write geometry section
        print >> fp, 'GEOMETRY'
        print >> fp, 'model:                                  ',\
                     '%s.ensight/geometry' %self.simname
        print >> fp

        #  Write variable section
        print >> fp, 'VARIABLE'
        for x in self.vars:
            if(x.rank > 0):
	        #print out variables that live on a mesh only	
		if (x.mesh != ''):
		    #print out variables belong to a predefined meshspace	
                    if (self.mapLookup.has_key(x.meshspace)) :
                        writevar     = self.writeThisVar(x)
                        if writevar:
	                    newname  = string.replace(x.nickname,'(','_')
        	      	    newname  = string.replace(newname,')',"_")
                	    filename = '%s.ensight/%s.*****'%(self.simname,newname)
	                    if (self.mapLookup[x.meshspace] == 0):
            		        varfor = 'element:'
                	    else:
                    		varfor = 'node:   '
                	    if (x.rank == 1):
                    		vartype = 'scalar'
                	    elif (x.rank == 2 and x.shape[1] == 3):
                    		vartype = 'vector'
                	    else:
                    		vartype = 'tensor'
                	    nb = 16 - len(newname)
                	    if (nb < 0):
                    		nb = 0
                	    varname = x.nickname + ' '*nb
                	    print >> fp, '%s per %s    %s %s' \
                              %(vartype,varfor,varname,filename)
        print >> fp

        if self.seq_no == 0:
            increment = self.cycle
            start     = increment
        else:
            increment = self.cycle/self.seq_no
            start     = 0

        #  Write time section
        #  Write time section
        print >> fp, 'TIME'
        print >> fp, 'time set:               1'
        print >> fp, 'number of steps:      %3d' %(self.seq_no+1)

        i = 0
        file_val = unique(file_val)
        file_val.sort()
        for x in file_val:
            if (i == 0):
                print >> fp, 'filename numbers:       %05d' %(x),
            else:
                print >> fp, ' %05d' %(x),
            i += 1

        print >> fp

        i = 0
        time_val = unique(time_val)
        time_val.sort()
        for x in time_val:
            if (i == 0):
                print >> fp, 'time values:            %12.5E' %(x),
            else:
                print >> fp, ' %12.5E' %(x),
            i += 1

        fp.close()


    def writebMesh(self):

        #Write EnSight geometry (mesh) file.
        global mesh_type
        global ncblocks
        global blocksmin
        global cellindices
        global vertexindices
        global gapelemblock

        filename = '%s/geometry'%(self.dirname)

        self.FP  = fortransupport.F95Binary()
        fp       = self.FP
        fp.open(filename,'w',byteSwap=0)

        # Note hardcoded mesh name and type
        mesh = self.mesh
        if (mesh.type == 'HEX'): mesh_type = 'hexa8'
        if (mesh.type == 'TET'): mesh_type = 'tetra4'

        if mesh.cells.has_key('BLOCKID'):
            cblocks   = mesh.cells['BLOCKID']
            ncblocks  = cblocks[Numeric.argmax(cblocks)]
            blocksmin = cblocks[Numeric.argmin(cblocks)]
        else:
            ncblocks  = 1
            cblocks   = Numeric.array([1],typecode='i')
            cblocks   = Numeric.resize(cblocks,[mesh.cells['N']])
            blocksmin = 1

        # Write geometry file info..

        desc = 'Fortran Binary'
        S    = self.blankfill(desc,self.strsize)
        if self.debug: print >> self.fpwatch, type(S),S
        self.setNandWrite(S,1)

        desc = 'geometry (Truchas: "%s.ensight/geometry")' %(self.simname)
        S    = self.blankfill(desc,self.strsize)
        if self.debug: print >> self.fpwatch, type(S),S
        self.setNandWrite(S,1)

        desc = '%s.inp' %(self.simname)
        S    = self.blankfill(desc,self.strsize)
        if self.debug: print >> self.fpwatch, type(S),S
        self.setNandWrite(S,1)

        desc = 'node id assign'
        S    = self.blankfill(desc,self.strsize)
        if self.debug: print >> self.fpwatch, type(S),S
        self.setNandWrite(S,1)

        desc = 'element id assign'
        S    = self.blankfill(desc,self.strsize)
        if self.debug: print >> self.fpwatch, type(S),S
        self.setNandWrite(S,1)

        #dictionary to store cell and vertex indices for each mesh block (part) number
        cellindices   = {}
        vertexindices = {}
        gapelemblock  = []

        c      = mesh.cells
        cVerts = c['VERTICES']
        v      = mesh.vertices
        vt     = Numeric.transpose(v['COORDS'])

        if ncblocks > 1:
            blocks = []
            for i in range(blocksmin,ncblocks+1):
                cellindices[i]   = []
                vertexindices[i] = []
                thisblock        = getBlockIndicies(cblocks,i,v['N'],cVerts,fp=self.fpwatch,debug=self.debug)
                thisblock.getData()
                blocks.append(thisblock)
                cellindices[i]   = thisblock.cellids
                vertexindices[i] = thisblock.vtxids
                checkgap         = len(cellindices[i])<1 or len(vertexindices[i])<1
                gapelemblock.append(checkgap)
        else:
            cellindices[1]   = range(c['N']+1)[1:len(range(c['N']+1))]
            vertexindices[1] = range(v['N']+1)[1:len(range(v['N']+1))]
            gapelemblock.append(0)

        cnt = 0
        for i in range(blocksmin,ncblocks+1):

            # do not write out gap element blocks
            if not gapelemblock[cnt]:
                
                desc = 'part'
                S    = self.blankfill(desc,self.strsize)
                if self.debug: print >> self.fpwatch, type(S),S
                self.setNandWrite(S,1)

                S    = Numeric.array([i],typecode='i')
                if self.debug: print >> self.fpwatch, type(S),S
                self.setNandWrite(S,1)

                desc = '%s %s %i' %(self.simname,'block',i)
                S    = self.blankfill(desc,self.strsize)
                if self.debug: print >> self.fpwatch, type(S),S
                self.setNandWrite(S,1)

                desc = 'coordinates'
                S    = self.blankfill(desc,self.strsize)
                if self.debug: print >> self.fpwatch, type(S),S
                self.setNandWrite(S,1)

                # write coordinates and cell-vertices in this mesh block
                # need total number of cells, and also total number of vertices for all cells.

                if ncblocks > 1:

                    S = Numeric.array([len(vertexindices[i])],typecode='i')
                    if self.debug: print >> self.fpwatch, 'vertexindices',i,type(S),S
                    self.setNandWrite(S,1)
                    for imv in range(len(vt)):
                        tmp = []
                        for j in vertexindices[i]:
                            tmp.append(vt[imv,j-1])
                        tmparray = Numeric.array(tmp,typecode='f')
                        fp.write(tmparray,endRecord=1,debug=1)

                    desc        = mesh_type
                    S           = self.blankfill(desc,self.strsize)
                    if self.debug: print >> self.fpwatch, type(S), S                
                    self.setNandWrite(S,1)

                    S    = Numeric.array([len(cellindices[i])],typecode='i')
                    if self.debug: print >> self.fpwatch, 'cellindices',i,len(cellindices[i]),type(S), S        
                    self.setNandWrite(S,1)

                    thisblock        = blocks[cnt]
                    fp.write(thisblock.blckconnectivity,endRecord=1,debug=1)
                
                else:

                    S    = Numeric.array([v['N']],typecode='i')
                    if self.debug: print >> self.fpwatch, type(S), S
                    self.setNandWrite(S,1)
                    vt   = Numeric.transpose(v['COORDS'])
                    for imv in range(len(vt)):
                        tmparray = vt[imv].astype('f')
                        fp.write(tmparray,endRecord=1,debug=1)
                    desc        = mesh_type
                    S           = self.blankfill(desc,self.strsize)
                    if self.debug: print >> self.fpwatch, type(S), S                
                    self.setNandWrite(S,1)
                    S    = Numeric.array([c['N']],typecode='i')
                    if self.debug: print >> self.fpwatch, type(S), S        
                    self.setNandWrite(S,1)
                    tmparray = cVerts.astype('i')
                    fp.write(tmparray,endRecord=0,debug=1)


            cnt += 1

        fp.close()


    def writebVariables(self):

        s = ''
        for x in self.vars:
            if (x.rank > 0 ):
	        #print out variables that live on a mesh only	
		if (x.mesh != ''):
		    #print out variables belong to a predefined meshspace	
                    if (self.mapLookup.has_key(x.meshspace)) :
                        writevar = self.writeThisVar(x)
                        if writevar:
                	   newname  = string.replace(x.nickname,'(','_')
                	   newname  = string.replace(newname,')',"_")
                	   filename = '%s/%s.%05d'%(self.dirname,newname,self.cycle)
                           self.FP  = fortransupport.F95Binary()
                           fp       = self.FP
                           fp.open(filename,'w',byteSwap=0)
                           desc = '%s (Truchas: '"'%s_output/%s.inp'"')' \
                                  %(x.nickname,self.simname,self.simname)
                           desc = str(desc)
                           S    = self.blankfill(desc,self.strsize)
                           self.setNandWrite(S,1)
                           cnt  = 0
                           for i in range(blocksmin,ncblocks+1):
                               
                               # do not write out variables on gap element blocks
                               if not gapelemblock[cnt]:

                                   desc   = 'part'
                                   S      = self.blankfill(desc,self.strsize)
                                   self.setNandWrite(S,1)
                                   partno = i
                                   S      = Numeric.array([partno],typecode='i')
                                   self.setNandWrite(S,1)
                                   if (self.mapLookup[x.meshspace] == 0):
                                       desc   = mesh_type
                                       S      = self.blankfill(desc,self.strsize)
                                       self.setNandWrite(S,1)
                                       if self.debug: print >> self.fpwatch, type(S),S
                                       tmp = []
                                       for j in cellindices[i]:
                                           tmp.append(x.data[j-1])
                                       tmparray = Numeric.array(tmp,typecode=x.type)
                                   elif (self.mapLookup[x.meshspace] == 1):
                                       desc   = 'coordinates'
                                       S      = self.blankfill(desc,self.strsize)
                                       self.setNandWrite(S,1)
                                       if self.debug: print >> self.fpwatch, type(S),S
                                       tmp = []
                                       for j in vertexindices[i]:
                                           tmp.append(x.data[j-1])
                                       tmparray = Numeric.array(tmp,typecode=x.type)
                                   tmparray = Numeric.transpose(tmparray)
                                   self.writeByComponent(tmparray,x.rank,fp)

                               cnt += 1

                                   
                           fp.close()

        return


    #now define any ENSIGHT ascii routines


    def writeMesh(self):

        #Write EnSight geometry (mesh) file.

        global mesh_type
        global ncblocks
        global blocksmin
        global cellindices
        global vertexindices
        global gapelemblock

        filename = '%s/geometry'%(self.dirname)
        fp       = open(filename, 'w')

        # Note hardcoded mesh name and type
        mesh = self.mesh
        if (mesh.type == 'HEX'): mesh_type = 'hexa8'
        if (mesh.type == 'TET'): mesh_type = 'tetra4'

        if mesh.cells.has_key('BLOCKID'):
            cblocks   = mesh.cells['BLOCKID']
            ncblocks  = cblocks[Numeric.argmax(cblocks)]
            blocksmin = cblocks[Numeric.argmin(cblocks)]
        else:
            ncblocks  = 1
            cblocks   = Numeric.array([1],'i')
            cblocks   = Numeric.resize(cblocks,[mesh.cells['N']])
            blocksmin = 1

        # Write geometry file info.
        print >> fp, 'geometry (Truchas: "%s.ensight/geometry")' %(self.simname)
        print >> fp, '%s.inp' %(self.simname)
        print >> fp, 'node id assign'
        print >> fp, 'element id assign'

        #dictionary to store cell and vertex indices for each mesh block (part) number
        cellindices   = {}
        vertexindices = {}
        gapelemblock  = []

        v      = mesh.vertices
        c      = mesh.cells
        cVerts = c['VERTICES']
        vt     = Numeric.transpose(v['COORDS'])

        if ncblocks > 1:
            blocks = []
            for i in range(blocksmin,ncblocks+1):
                cellindices[i]   = []
                vertexindices[i] = []
                thisblock        = getBlockIndicies(cblocks,i,v['N'],cVerts,fp=self.fpwatch,debug=self.debug)
                thisblock.getData()
                blocks.append(thisblock)
                cellindices[i]   = thisblock.cellids
                vertexindices[i] = thisblock.vtxids
                checkgap         = len(cellindices[i])<1 or len(vertexindices[i])<1
                gapelemblock.append(checkgap)
        else:
            cellindices[1]       = range(c['N']+1)[1:len(range(c['N']+1))]
            vertexindices[1]     = range(v['N']+1)[1:len(range(v['N']+1))]
            gapelemblock.append(0)

        cnt = 0
        for i in range(blocksmin,ncblocks+1):

            # do not write out gap element blocks
            if not gapelemblock[cnt]:

                print >> fp, 'part'
                print >> fp, '%10d' %(i)

                desc = '%s %s %i' %(self.simname,'block',i)
                print >> fp, desc
            
                print >> fp, 'coordinates'

                # write coordinates and cell-vertices in this mesh block
                # need total number of cells, and also total number of vertices for all cells.

                if ncblocks > 1:
        
                    print >> fp, '%5d' %(len(vertexindices[i]))

                    for imv in range(len(vt)):
                        tmp = []
                        for j in vertexindices[i]:
                            tmp.append(vt[imv,j-1])
                        tmparray = Numeric.array(tmp,typecode='f')
                        self.writeArray(1,tmparray,'',self.dformat,fp,nperline=1)
                
                    desc        = mesh_type
                    tmp = '%10d' %(len(cellindices[i]))
                    print >> fp, desc + '\n' + tmp

                    thisblock = blocks[cnt]
                    self.writeArray(2, thisblock.blckconnectivity, '','%10d', fp, nperline=8)

                else:

                    print >> fp, '%10d' %(v['N'])

                    for imv in range(len(vt)):
                        tmp = []
                        for j in vertexindices[i]:
                            tmp.append(vt[imv,j-1])
                        tmparray = Numeric.array(tmp,typecode='f')
                        self.writeArray(1,tmparray,'',self.dformat,fp,nperline=1)
                    
                    desc     = mesh_type
                    print >> fp, desc

                    print >> fp, '%10d' %(c['N'])

                    self.writeArray(2,cVerts, '','%10d', fp, nperline=8)

            cnt += 1
 
        fp.close()        


    def writeCase(self):

        #  Write EnSight Case file
        global time_val
        global file_val

        filename = '%s.ensight.CASE'%(self.fileName)
        fp       = open(filename, 'w')
        if self.seq_no == 0:
            time_val = []
            file_val = []
        time_val.append(self.t)
        file_val.append(self.cycle)

        #  Write header section
        print >> fp, 'FORMAT'
        print >> fp, 'type:                   ensight gold'
        print >> fp

        #  Write geometry section
        print >> fp, 'GEOMETRY'
        print >> fp, 'model:                                  %s.ensight/geometry' \
                     %self.simname
        print >> fp

        #  Write variable section
        print >> fp, 'VARIABLE'
        for x in self.vars:
            if(x.rank > 0):
	        #print out variables that live on a mesh only	
		if (x.mesh != ''):
		    #print out variables belong to a predefined meshspace	
                    if (self.mapLookup.has_key(x.meshspace)) :
                        writevar     = self.writeThisVar(x)
                        if writevar:
	                    newname  = string.replace(x.nickname,'(','_')
        	      	    newname  = string.replace(newname,')',"_")
                	    filename = '%s.ensight/%s.*****'%(self.simname,newname)
	                    if (self.mapLookup[x.meshspace] == 0):
            		        varfor = 'element:'
                	    else:
                    		varfor = 'node:   '
                	    if (x.rank == 1):
                    		vartype = 'scalar'
                	    elif (x.rank == 2 and x.shape[1] == 3):
                    		vartype = 'vector'
                	    else:
				if (x.shape[1] == 6): vartype = 'tensor symm'
				else: vartype = 'tensor'
                	    nb = 16 - len(newname)
                	    if (nb < 0):
                    		nb = 0
                	    varname = x.nickname + ' '*nb
                	    print >> fp, '%s per %s    %s %s' \
                              %(vartype,varfor,varname,filename)
        print >> fp

        if self.seq_no == 0:
            increment = self.cycle
            start     = increment
        else:
            increment = self.cycle/self.seq_no
            start     = 0

        #  Write time section
        #  Write time section
        print >> fp, 'TIME'
        print >> fp, 'time set:               1'
        print >> fp, 'number of steps:      %3d' %(self.seq_no+1)

        i = 0
        file_val = unique(file_val)
        file_val.sort()
        for x in file_val:
            if (i == 0):
                print >> fp, 'filename numbers:       %05d' %(x),
            else:
                print >> fp, ' %05d' %(x),
            i += 1

        print >> fp

        i = 0
        time_val = unique(time_val)
        time_val.sort()
        for x in time_val:
            tmp = self.dformat%(x)
            if (i == 0):               
                print >> fp, 'time values:            ' + tmp,
            else:
                print >> fp, tmp,
            i += 1

        fp.close()


    def writeVariables(self):

        s = ''
        for x in self.vars:
            if(x.rank > 0):
	        #print out variables that live on a mesh only	
		if (x.mesh != ''):
		    #print out variables belong to a predefined meshspace	
                    if (self.mapLookup.has_key(x.meshspace)) :
                        writevar     = self.writeThisVar(x)
                        if writevar:
                	   newname  = string.replace(x.nickname,'(','_')
                	   newname  = string.replace(newname,')',"_")
                	   filename = '%s/%s.%05d'%(self.dirname,newname,self.cycle)
                           fp       = open(filename, 'w')
                           print >> fp, '%s (Truchas: '"'%s_output/%s.inp'"')' \
                                 %(x.nickname,self.simname, self.simname)
                           
                           cnt      = 0
                           for i in range(blocksmin,ncblocks+1):

                               # do not write out variables on gap element blocks
                               if not gapelemblock[cnt]:
                               
                                   print >> fp, 'part'
                                   print >> fp, '%10d' %(i)
                                   if (self.mapLookup[x.meshspace] == 0):
                                       desc   = mesh_type
                                       print >> fp, desc
                                       tmp = []
                                       for j in cellindices[i]:
                                           tmp.append(x.data[j-1])
                                       tmparray = Numeric.array(tmp,typecode=x.type)
                                   elif (self.mapLookup[x.meshspace] == 1):
                                       desc   = 'coordinates'
                                       print >> fp, desc
                                       tmp = []
                                       for j in vertexindices[i]:
                                           tmp.append(x.data[j-1])
                                       tmparray = Numeric.array(tmp,typecode=x.type)
                                   tmparray = Numeric.transpose(tmparray)
                                   self.writeByComponent(tmparray,x.rank,fp,transpose=0,ascii=1,
                                                         format = self.dformat,nperline=1)

                               cnt += 1
                               
                           fp.close()

        return

    def writeByComponent(self, data, rank, fp,
                         transpose=0,ascii=0,format='%12.5f',nperline=10):

        try:
           import numpy.oldnumeric as Numeric
        except ImportError:
           import Numeric
        except:
           raise
        if(transpose):
            tmp = Numeric.transpose(data)
        else:
            tmp = data

        if(rank <= 1):
            if ascii:
                self.writeArray(rank,data,'',format,fp,nperline)
                #print >> fp, ''
            else:
                tmparray = data.astype('f')
                self.FP.write(tmparray,endRecord=1,debug=1)
        else:
            for x in tmp:
                self.writeByComponent(x, rank-1, fp,
                                      transpose=0, ascii=ascii,
                                      format=format, nperline=nperline)
        return

    def finalizeAndClose(self):

        self.seq_no = self.seq_no + 1

        if (self.seq_no == 1):
            print >> self.fpwatch
            print >> self.fpwatch,'EnSight ascii file %s.ensight.CASE written.'\
	                           % (self.simname)
            print >> self.fpwatch, 'Directory %s.ensight created.' % \
			            (self.simname)
            print >> self.fpwatch
        print >> self.fpwatch, 'Contains files for time %f.' %(self.t)

if __name__ =='__main__':

    'component testing ENSIGHT write component'

    dfltdir = '../../scripts/test_TBrookParse/samples/'
    
    from PYTHONutils import uTestOpts

    # file, debug, output prefix, binary, ascii, clean
    opts = uTestOpts('fdoc', 
                     defaults = {'c' : False,
                                 'o' : 'ENSwriteUnitTest'},
                     actions  = {'c' : 'store_true'},
                     dir      = dfltdir)
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fpwatch    = sys.stdout
    try:
        if opt.c:
            import os
            print >> fpwatch, "Cleaning up from a past test run...\n",\
                  "  (removing %s(*) directories)..."  %opt.o
            cmd = 'rm -rf %s' %(dfltdir + opt.o + '.ascii.*')
            os.system(cmd)
            cmd = 'rm -rf %s' %(dfltdir + opt.o + '.binary.*')
            os.system(cmd)
            print >> fpwatch, "Cleanup successful!!"

        else:
            from XMLgetStorageObject import getStorageObject
        
            storage   = getStorageObject(opt.f,debug=1)	
            mesh      = storage.mlist[0] #assume default mesh

            " first get all data values associated with this mesh"
            for meshspace in mesh.mslist:
                storage.getValues(meshspace.vlist)
                
                
            "fill the mesh with is actual coorindate values"
            mesh.fillMesh()
        

            itimeout = 0
            t        = storage.tlist[itimeout].time
            vars     = storage.tlist[itimeout].vlist
            storage.getValues(vars)

            #ASCII write
            print >> fpwatch, '*'*60
            print >> fpwatch, 'Testing ascii write...' 
            fileout = dfltdir + opt.o + '.ascii'
            writer  = ENSIGHTwriteStorageObject(fileout, 'ascii', mesh, t,
                                                itimeout, vars,
                                                dformat = '%12.4e',
                                                fpwatch = fpwatch,
						debug   = opt.d)
            print 'ASCII ENSIGHT files %s\n  written to %s' \
                      %(opt.o+'.ascii',dfltdir)
                
            #BINARY write
            print >> fpwatch, '*'*60
            print >> fpwatch, 'Testing binary write...'
            fileout = dfltdir + opt.o + '.binary'
            writer  = ENSIGHTwriteStorageObject(fileout, 'binary', mesh, t,
                                                itimeout, vars, fpwatch=fpwatch,
						debug=opt.d)
            print >> fpwatch, 'Binary ENSIGHT files %s\n  written to %s' \
                              % (opt.o+'.binary',dfltdir)
                
            print >> fpwatch, '*'*60
            print >> fpwatch, \
	        '--> Note: to clean up, re-execute same command line plus "-c"'

    except:
        import sys
        print >> fpwatch, "---> Test failed in some aspect <---"
        print >> fpwatch, "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise
    



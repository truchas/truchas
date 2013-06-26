"""

 writeStorageObject

 -----------------------------------------------------------------------------
  Purpose:
 
  Write out a storage object in a standard ascii and binary form

    The arguments needed to write out this storage object are:

      fileName:   Name of output gmv file
      storage:    The storage structure (typically contains mesh, timesteps information)
      mesh:       The mesh to be outputted
      t:          Problem time
      bVariables: A list of variables that we want printed out
      bFlags:     A list of variables that are to be used as flags
                   in an as yet undetermined format
      bRegions:   A list of regions in an as yet
                 undetermined format
    Currently:
      * bFlags and bRegions are ignored.
      * The mesh is assumed to be hex
      * only python is used to write ascii gmv files

  Public Interface(s):

    None.  This is a baseClass used by the ENSIGHT, GMV, TECPLOT, RESTART and VTK
    storage object writers.
    
  Contains:
    class writeStorageObject
        bImport(self)
        binaryWrite(self)
        writebHeader(self)
        writebCase(self)
        writebMesh(self)
        writebVariables(self)
        writebFlags(self,flags)
        writebGroups(self,flags)
        bfinalizeAndClose(self)

        blankfill(self,str,stringsize)
        setNandWrite(self,SS,iend)
        writeByComponent(self, data, rank, transpose=0)
        findAndWrite(self, varname, transpose=0, byComponent=0,debug=0)

        asciiWrite(self)
        writeHeader(self)
        writeCase(self)
        writeMesh(self)
        writeVariables(self)
        writeFlags(self, flags)
        writeGroups(self, flags)
        finalizeAndClose(self)

        writeArray(self,rank,data,hdr,format,fp,nperline=7,sin='',
                   tprint=0,aprint=1,ibuflen=4096)
        writeCompsArray(self,data,comps,format,fp,nperline=7)
        writeReducedFieldAndRegion(self,face,field,indices,posns,hdr,fp,debug=0)
        writeFieldAndRegion(self,face,field,indices,posns,hdr,fp,debug=0)
        writeStatistics(self,storage,column=15)
        consistencyCheck(self)
        writeThisVar(self,x)
        rjustln(self,column,s)

    ---> NO Unit Test Block

  Version:
    $ID$

  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""

import os, sys
try:
   import numpy.oldnumeric as Numeric
except ImportError:
   import Numeric
except:
   raise

class writeStorageObject:

    triedascii = 0

    def __init__(self,fpwatch=sys.stdout,debug=0): 

       self.fpwatch     = fpwatch
       self.debug       = debug
       self.module_name = 'asciiwriter'
       iform            = '7d'
       self.iformat     = '%%%s' %(iform)

    def binaryWrite(self):

        print >> self.fpwatch, 'Writing binary '+ self.fileName

        sys.stdout.flush()

        # First check all variables live on the chosen mesh
        if self.consistencyCheck(): 
           return	

        "Need to import binaryWriter from the particular utils directory"
        self.bImport()
        
        "Write binary the header"
        self.writebHeader()
        
        "Write the binary mesh"
        if self.seq_no == 0:
           self.writebMesh()

        "Write a case file (always in ascii format)"
        self.writebCase()

        "Write the binary variables"
        self.writebVariables()

        "Write the binary flags and groups"
        self.writebFlags(self.bFlags)
        self.writebGroups(self.bFlags)
        
        "Close the binary file"
        self.bfinalizeAndClose()

        # return success
        return

    def bImport(self):

        return

    def writebHeader(self):

        return

    def writebCase(self):

        return

    def writebMesh(self):

        return

    def writebVariables(self):

        return

    def writebFlags(self,flags):

        return

    def writebGroups(self,flags):

        return
    

    def bfinalizeAndClose(self):

	self.FP.close()

	return

    def blankfill(self,str,stringsize):

        #Fill string variables to length with blanks

        format = '%%-%ds'%stringsize

        if len(str)>stringsize:
           str = str[0:stringsize-1]

        return format%str

    def setNandWrite(self,SS,iend):

        fp = self.FP
        SN = Numeric.array(SS)
        fp.write(data=SN,endRecord=iend,debug=self.debug)
        return

    def writeByComponent(self, data, rank, transpose=0):
        
        if(transpose):
            tmp = Numeric.transpose(data)
        else:
            tmp = data
            
        if(rank <= 1):
            self.FP.write(data,endRecord=1)
        else:
            for x in tmp:
                self.writeByComponent(x, rank-1, transpose=0)
        return
            

    def findAndWrite(self, varname, transpose=0, byComponent=0, debug=0):
        
        for x in self.vars:
            if (x.name == varname):
                if (transpose):
                    tmp = Numeric.transpose(x.data)
                else:
                    tmp = x.data

                if (self.debug):
                    print '    ',x.name, x.shape, x.rank, transpose, len(tmp), len(x.data)
                
                if (byComponent):
                    self.writeByComponent(tmp, x.rank)
                else:
                    self.FP.write(tmp,endRecord=1)
                return

    def asciiWrite(self):

        print >> self.fpwatch, 'Writing ascii '+ self.fileName,

        sys.stdout.flush()
        self.FP          = open(self.fileName, 'w')

	# First check all variables live on the chosen mesh
	if self.consistencyCheck(): 
	   return	

        # Write the header
        self.writeHeader()

        # Write the mesh
	if self.seq_no == 0:
           self.writeMesh()

	# Write a case file
	self.writeCase()

        # Write the variables
        self.writeVariables()

        # Write the flags
        self.writeFlags(self.bFlags)
        self.writeGroups(self.bFlags)
        
        # Close the file
        self.finalizeAndClose()

	print >> self.fpwatch, ' done'
        # return success
        return

    def writeHeader(self):

	#for the base writer just write a simple header

        print >> self.FP, 'ascii'
        

    def writeArray(self,rank,data,hdr,format,fp,nperline=7,sin='',
                   tprint=0,aprint=1,ibuflen=4096):

        # writes out a data array in a specified ascii format
        from time import time
        import array
        import sys

        if data==None:
           return
        
        if(tprint): t1 = time()

        if ( hdr == None): hdr = ''
        if ( sin == None): sin = ''
        if(not self.triedascii):

            try:
                Estring     = 'WRITEutils import error'
                import asciiwriter
            except:
                print >> self.fpwatch
                print >> self.fpwatch, '  _________________________________'
                print >> self.fpwatch, '  ---------------------------------'
                print >> self.fpwatch
                print >> self.fpwatch, '  ERROR: '+Estring
                print >> self.fpwatch
                print >> self.fpwatch, '  Using python based ascii writer for '
                print >> self.fpwatch, '  all ascii array writes'
                print >> self.fpwatch, '  ---------------------------------'
                print >> self.fpwatch
                self.triedascii = 1

            myShape = Numeric.array(Numeric.shape(data),'i')
            
            "Accomodate Numpy/Numeric differences"
            try:    # Numpy
               dtyp   = data.dtype.char
               myData = data.flat
               myData = Numeric.reshape(myData,myShape)
            except: # Numeric
               dtyp   = data.typecode()
               myData = Numeric.array(data)

            try:
                Estring = 'Array writing error'
                asciiwriter.writeArrayASCII(dtyp, rank, myShape, myData, hdr+sin, format, nperline, fp, 0)
                if(tprint):
                    t2 = time()
                    print >> self.fpwatch, '___time taken is: ',t2-t1
                return
            except:
                print >> self.fpwatch
                print >> self.fpwatch, '  _________________________________'
                print >> self.fpwatch, '  ---------------------------------'
                print >> self.fpwatch
                print >> self.fpwatch, '  ERROR: %s' %(Estring)
                print >> self.fpwatch
                print >> self.fpwatch, '  Using python based ascii writer for '
                print >> self.fpwatch, '  this array type (%s) only.'%(dtyp)
                print >> self.fpwatch, '  ---------------------------------'
                print >> self.fpwatch

        try:


            if ( rank <= 0 ):
                print >> fp,data

            elif (rank == 1):
                s      = hdr + sin
                iCount = 0
                for v in data:
                    s      += format%v
                    iCount += 1
                    if ( iCount%nperline == 0 and iCount < len(data)-1):
                        s  += '\n'
                        if(aprint or len(s) > ibuflen):
                            fp.write(s)
                            s  = ''
                s += '\n'
                if(aprint):
                    fp.write(s)
            else:
                i = 0
                s=''

                if(0 and (hdr == 'hex 8')):
                    # hot fix for ascii writing of cell connectivity
                    s = ''
                    fmt = hdr + format*8
                    for v in data:
                        s += fmt%(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7])
                        if(len(s)>ibuflen):
                            fp.write(s)
                            s = ''
                    fp.write(s)
                    s = ''
                else:
                    "Accomodate Numpy/Numeric syntax difference"
                    try:
                       dtyp = data.dtype.char
                    except:
                       dtyp = data.typecode()
                    dumData = Numeric.array(data,dtyp)
                    for v in dumData:
                        i=i+1
                        # Recursively call one rank lower
                        s=self.writeArray(rank-1,v,hdr,format,fp,nperline,sin,
                                          tprint=0,aprint=0)
                        if (aprint or len(s)>ibuflen):
                            fp.write(s)
                            s = ''
                    if (aprint or len(s)>ibuflen):
                        fp.write(s)
                        s = ''

            if(not aprint): return s
            if(tprint):
                t2 = time()
                print >> self.fpwatch, '___time taken is: ',t2-t1

        except:
            print >> self.fpwatch
            print >> self.fpwatch
            print >> self.fpwatch, 'error in writearray'
            print >> self.fpwatch
            print >> self.fpwatch
            sys.exit(5)

    def writeCompsArray(self,data,comps,format,fp,nperline=7):

        # writes out a rank=2 data array in a specified ascii format
        # with header = [comps] containing each component field name
        import array
        s      = ''
        cCount = 0
        for v in data:
            if (len(comps) > cCount):
                self.writeArray(1,v,None,format,fp,nperline,comps[cCount])
                print >> fp
                cCount += 1
            
    def writeReducedFieldAndRegion(self,face,field,indices,posns,hdr,fp,
                                   coordsfmt = '%12.5e',fieldfmt = '%20.5e'):

        #writes out a field defined only on the region
        #along with the region's indices and posns

        
        if field.rank > 2:
            errstring  = "\n In %s rank = %d but face = %d \n" \
                         %(self.__class__.__name__,field.rank,face)
            assert face > 0 and face < 7, errstring
            myField    = Numeric.array(field.data[:,face-1,:])
            myShape    = []
            myShape.append(field.shape[0])
            myShape.append(field.shape[2])
            myShape    = Numeric.array(myShape,'i')
            myRank     = 2
        else:
            myField = Numeric.array(field.data)
            myShape = Numeric.array(Numeric.shape(field.data),'i')
            myRank  = field.rank

        "Accomodate Numpy/Numeric differences"
        try:    # Numpy
           dtyp   = myField.dtype.char
           myData = myField.flat
           myData = Numeric.reshape(myData,myShape)
        except: # Numeric
           dtyp   = myField.typecode()
           myData = Numeric.array(myField)

        lenindices = len(indices)

        try:
            import asciiwriter
            cformat    = '%s       %s %s %s       '%(self.iformat,coordsfmt,coordsfmt,coordsfmt)            
            asciiwriter.writeQuery(myRank, myShape, face, myData, indices, lenindices, posns, hdr, fp, cformat, fieldfmt, self.debug)
        except:
            print >> self.fpwatch
            print >> self.fpwatch, 'Ascii C writer failed in %s' %(self.__class__.__name__)
            print >> self.fpwatch, 'Writing out reduced field and region arrays in Python.',\
                  'This will be slow!'
            print >> self.fpwatch

            stt         = hdr
            thiscformat = '%s,%s,%s'%(coordsfmt,coordsfmt,coordsfmt)
            
            if (field.rank <=1):

                cnt = 0
                for indice in indices:
                    xpos = posns[indice-1][0]
                    ypos = posns[indice-1][1]
                    zpos = posns[indice-1][2]
                        
                    pos  = thiscformat %(xpos,ypos,zpos)

                    thisdata = fieldfmt%myData[cnt]
                    stt     += ' %8i %+30s %s' \
                           %(indice,pos,thisdata)
                    stt += '\n'
                    cnt += 1

            elif(field.rank == 2):

                cnt = 0
                for indice in indices:
                    xpos = posns[indice-1][0]
                    ypos = posns[indice-1][1]
                    zpos = posns[indice-1][2]
                        
                    pos  = thiscformat %(xpos,ypos,zpos)

                    stt += ' %8i %+30s     ' %(indice,pos)
                    for e in myData[cnt]:
                       thisdata = fieldfmt%e
                       stt     += thisdata
                    stt += '\n'

                    cnt += 1
                    
            elif(field.rank == 3):

                cnt = 0 
                for indice in indices:
                    xpos = posns[indice-1][0]
                    ypos = posns[indice-1][1]
                    zpos = posns[indice-1][2]
                        
                    pos  = thiscformat %(xpos,ypos,zpos)
                        
                    stt += ' %8i %+30s     ' %(indice,pos)
                    count = 0
                    for e in myData[cnt]:
                       count += 1
                       if count == face:
                          for f in e:
                             thisdata = fieldfmt%f
                             stt     += thisdata
                    stt += '\n'

                    cnt += 1
                    
            s = stt +'\n'

            fp.write(s)


    def writeFieldAndRegion(self,face,field,indices,posns,hdr,fp,
                            coordsfmt = '%12.5e',fieldfmt = '%20.5e'):
                                    
        #writes out field on a region along with the region's indices and posns

        import string
        
        if field.rank > 2:
            errstring  = "\n In %s rank = %d but face = %d \n" \
                         %(self.__class__.__name__,field.rank,face)
            assert face > 0 and face < 7, errstring
            myField    = Numeric.array(field.data[:,face-1,:])
            myShape    = []
            myShape.append(field.shape[0])
            myShape.append(field.shape[2])
            myShape    = Numeric.array(myShape,'i')
            myRank     = 2
        else:
            myField = Numeric.array(field.data)
            myShape = Numeric.array(Numeric.shape(field.data),'i')
            myRank  = field.rank

        "Accomodate Numpy/Numeric differences"
        try:    # Numpy
           myData = myField.flat
           myData = Numeric.reshape(myData,myShape)
        except: # Numeric
           myData = Numeric.array(myField)

        lenindices = len(indices)

        try:
           import asciiwriter
           cformat    = '%s       %s %s %s       '%(self.iformat,coordsfmt,coordsfmt,coordsfmt)
           asciiwriter.writeQuery(myRank, myShape, face, myData, indices, lenindices, posns, hdr, fp, cformat, fieldfmt, self.debug)
        except:
            print >> self.fpwatch
            print >> self.fpwatch, 'Ascii C writer failed in %s' %(self.__class__.__name__)
            print >> self.fpwatch, 'Writing out field and region arrays in Python. This will be slow!'
            print >> self.fpwatch

            stt         = hdr
            thiscformat = '%s,%s,%s'%(coordsfmt,coordsfmt,coordsfmt)
            
            if (field.rank <=1):

                for indice in indices:
                    xpos = posns[indice-1][0]
                    ypos = posns[indice-1][1]
                    zpos = posns[indice-1][2]
                        
                    pos  = thiscformat %(xpos,ypos,zpos)
                    
                    if(indice <= len(field.data) and indice > 0):

                        thisdata = fieldfmt%field.data[indice-1]
                        stt     += ' %8i %+30s %s' %(indice,pos,thisdata)
                        stt     += '\n'

            elif(field.rank == 2):

                for indice in indices:
                    xpos = posns[indice-1][0]
                    ypos = posns[indice-1][1]
                    zpos = posns[indice-1][2]
                        
                    pos  = thiscformat %(xpos,ypos,zpos)

                    if(indice <= len(field.data) and indice > 0): 
                        stt += ' %8i %+30s     ' %(indice,pos)
                        for e in field.data[indice-1]:
                            thisdata = fieldfmt%e
                            stt     += thisdata
                        stt += '\n'

            elif(field.rank == 3):

                for indice in indices:
                    xpos = posns[indice-1][0]
                    ypos = posns[indice-1][1]
                    zpos = posns[indice-1][2]
                        
                    pos  = thiscformat %(xpos,ypos,zpos)
                        
                    if(indice <= len(field.data) and indice > 0): 
                        stt += ' %8i %+30s     ' %(indice,pos)
                        count = 0
                        for e in field.data[indice-1]:
                            count += 1
                            if count == face:
                                for f in e:
                                    thisdata = fieldfmt%f
                                    stt     += thisdata
                        stt += '\n'
                    
            s = stt +'\n'

            fp.write(s)

    def writeStatistics(self,storage,column=15):

        #creates a string containing cycle, time, variable
        #statistics associated with storage object

        cycles = []
        times  = []
        probes = []        

        for step in storage.tlist:
            cycles.append(step.cycle)
            times.append(step.time)

        for probe in storage.plist:
            probes.append(probe.name)

        stats  = '\n'
        tmp    = 'Available output cycles : ' + str(cycles)
        stats += self.rjustln(column, tmp)
        tmp    = 'Available output times  : ' + str(times)
        stats += self.rjustln(column, tmp)
        stats += '\n'
        tmp    = 'Available probes        : ' + str(probes)
        stats += self.rjustln(column, tmp)
        stats += '\n'

        tmp    = 'Available variables     : '
        stats += self.rjustln(column, tmp)
        stats += '\n'
        
        maxl    = 0
        varlist = []

        for mesh in storage.mlist:
           for meshspace in mesh.mslist:
              if meshspace.name == 'VERTEX':
                 for var in meshspace.vlist:
                    if var.name == 'COORDS':
                       varlist.append(var)
              elif meshspace.name == 'CELL':
                 for var in meshspace.vlist:
                    if var.name == 'CENTROIDS' or var.name == 'VOLUMES': 
                       varlist.append(var)
                       maxl = len(var.name)

        if len(storage.tlist):        
           for var in storage.tlist[0].vlist:
              if var.shape[0] > 1:
                 maxl    = max(maxl,len(var.name))
                 varlist.append(var)

        if len(storage.aborts.nlrlist):
           for nlr in storage.aborts.nlrlist:
              for var in nlr.vlist:
                 if var.shape[0] > 1:
                    maxl    = max(maxl,len(var.name))
                    varlist.append(var)

        if len(storage.aborts.lrlist):
           for lr in storage.aborts.lrlist:
              for var in lr.vlist:
                 if var.shape[0] > 1:
                    maxl    = max(maxl,len(var.name))
                    varlist.append(var)        

        format  = '%s'
        colen   = 10 + 2*maxl - len('NICKNAME')
        format += '%%+%ds' %(colen)
        colen   = 2*maxl 
        format += '%%+%ds' %(colen)
        format += '%20s'
        format += '%25s'
        tmp     = format%('NAME', 'NICKNAME', 'MESH', 'SHAPE', 'MESHSPACE')
        stats += self.rjustln(column, tmp)
        stats += '\n'

        for var in varlist:
           if var.shape[0] > 1:
              colen   = 5 + 2*maxl - len(var.name)
              format  = '%s'
              format += '%%+%ds' %(colen)
              colen   = 5 + 2*maxl 
              format += '%%+%ds' %(colen)
              format += '%20s %20s'
              tmp     = format%(var.name,var.nickname,str(var.mesh),
                                str(var.shape),str(var.meshspace))
              stats  += self.rjustln(column, tmp)

        return stats

    def consistencyCheck(self):

        #check all chosen variables should live on the chosen mesh
        #if they don't then we have a consistency error and we must warn the user
        meshName  = self.mesh.name
        for x in self.vars:
            if(len(x.mesh) and meshName != x.mesh):
                print >> self.fpwatch
                print >> self.fpwatch, 'Inconsistent choice of variables and mesh !'
                print >> self.fpwatch, x.name
                print >> self.fpwatch, 'Input mesh:'
                print >> self.fpwatch, meshName
                print >> self.fpwatch, 'Mesh associated with this variable :'
                print >> self.fpwatch, x.mesh
                print >> self.fpwatch
                return(1)

        # success, return 0
        return(0)

    def writeFlags(self, flags):

        return

    def writeGroups(self, flags):

        return

    def writeCase(self):

        return
    
    def writeMesh(self):

        return
    
    def writeVariables(self):

        return
    
    def finalizeAndClose(self):
        
	#for the base writer just close the file

        self.FP.close()


    #now define auxiliary routines needed by visualisation writers

    def writeThisVar(self,x):

        #check to decide if we should write this variable x out
        """
        LJCox Note:
           This form is computationally more efficient (than the
           previous), since when ANY false condition is found, the
           work is done.  In the other form additional tests were
           required at each 'tmp = tmp and' entity.
           It is also easier to read.
        """

        tmp  =     (not x.name.startswith('M_')) \
               and ( x.name      is not 'MATSLOT') \
               and ('_STRAIN_'      not in x.name) \
               and ('_STRESS_'      not in x.name) \
               and ('Z_VF'          not in x.name) \
               and ('RHS'           not in x.name) \
               and ('MU'            not in x.name) \
               and ('SIGMA'         not in x.name) \
               and ('COIL'          not in x.name) \
               and ('sens_function' not in x.name) \
               and ('FREQ'          not in x.name) \
               and ('UHFS'          not in x.name) \
               and ('OLD'           not in x.name) \
               and ('Old'           not in x.name) 

        return tmp


    def rjustln(self,column,s):

        tmp = s.rjust(column + len(s))
        tmp += '\n'

        return tmp



 

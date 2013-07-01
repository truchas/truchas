import fortransupport
try:
    import numpy.oldnumeric as Numeric
except ImportError:
    import Numeric
except:
    raise

if __name__=='__main__':
    # Set sys.path for component testing mode
    import os, sys
    thisdir    = os.path.dirname(__file__)
    parserdir  = os.path.abspath(thisdir+'../')
    sys.path.append(parserdir)


class T_Binary_File:
    """

    This is the support class for reading binary files written by
    Truchas.  Some assumptions are made about the file.  First is that
    it is a truchas ts file written in binary format.  This means that
    the first record in the file is an integer 1, the second entry is
    a real_kind 1.0, and the third entry in the file is a logical
    true.  This provides me all the information I need to decode my
    size information.

    Provides:
      array1=a.get(typecode, shape=[1], offset=last_offset)
        This returns an array of type 'typecode' with shape specified
        by shape starting in the file with an offset last_offset.  In
        case the offset you want is not the last offset, 

    Here is how you would use this module:
      a=T_ts_File('filename')
      array1=a.get('i',[1], offset=200)
             will return an array of integers of single dimension
             having a single entry, basically a rank 0 at an offset
             of 200 bytes from the beginning of the file
      array2=a.get('d',[2,3])
             will return an array of doubles having dimension 2x3 from
             the current location of the file
                            
    """
    def __init__(self, filename, f95HeaderSize=4, debug=0):
        
        self.Init          = 0
        self.filename      = filename
        self.offset        = 0
        self.f             = None
        self.debug         = debug
        self.f95HeaderSize = f95HeaderSize #/4

        if filename == None:
            return None

        if(self.debug):
            print 'opening file using mode "rb"'
            print filename

        # First find out about byte swapping
        import array
        
        self.f = open(filename, 'rb')
        b      = array.array('i')
        
        if(self.debug):
            print 'created array'
            print 'b:', b
            
        b.fromfile(self.f,1) #self.f95HeaderSize)
        
        if(self.debug):
            print 'from file'
            print self.f
            print 'b: ', b
            
        if(b[0] > 256):
            # our first record is an integer 1, so
            # record size better be very small!
            self.byteSwapfp = 1
        else:
            self.byteSwapfp = 0
            
        if(self.debug):
            print '\nclosing file'
            print self.f

        self.f.close()

        # Now open the file as a F95Binary file
        self.f = fortransupport.F95Binary(self.f95HeaderSize,self.debug)
        
        if(self.debug):
            print '\nopening F95BinaryFile'
            
        self.f.open(filename,'r',self.byteSwapfp)

        # Integer size, confirm byteswapping
        r            = self.f.getRecords()[0]
        self.intSize = len(r)
        i            = Numeric.fromstring(r,'i',1)

        if(self.debug):
            print '\nfrom file' 
            print filename
            print 'integer size is :', self.intSize
            print 'fromstring is   :',i
        
        # The following is to check whether Numeric needs byte swapping
        try:
            ib = i.byteswap()[0]        
        except:
            ib = i.byteswapped()[0]

        if (self.debug):
            print '\ncheck byte swapping:'
            print 'ib is           :',ib
            
        if(ib == 1):
            self.byteSwap = 1
        elif(i[0] == 1):
            self.byteSwap = 0
        else:
            import sys
            print "\n\nNot a TS File \n"
            sys.exit(2)
            
        if(self.debug):
            print 'integer i is    : ', i[0]
            print 'record is       : ', r
            print 'byteswapping    : ', self.byteSwap
            
        if(not self.intSize == 4 ):
            print "\nintSize is    : ", self.intSize
            print "not 4, weirdness may happen.\n"
            
        self.offset += self.intSize

        if(self.debug):
            print 'got int size    : ', self.intSize 
            print 'got offset      : ', self.offset

        # real_kind size, confirm byteswapping
        r               = self.f.getRecords()[0]
        self.doubleSize = len(r)
        d               = Numeric.fromstring(r,'d',1)[0]

        if(self.debug):
            print '\nfrom file'
            print filename
            print 'double size is  :', self.doubleSize
            print 'fromstring is   :',d

        if(self.byteSwap):
            try:
                d = d.byteswap()[0]
            except:
                d = d.byteswapped()[0]

        if(self.debug):
            print 'double d is     : ', d
            print 'record is       : ', r
            print 'byteswapping    : ', self.byteSwap

        if(not (d == 1.0)):
            import sys
            print "\n\nPANIC: Byteswapping seems to be different for ints and double\n"
            sys.exit(2)
            
        if(not self.doubleSize == 8 ):
            print "\ndoubleSize is    : ", self.doubleSize
            print "not 8, weirdness may happen.\n"

        self.offset += self.doubleSize

        if(self.debug):
            print 'got double size : ', self.doubleSize 
            print 'got offset      : ', self.offset
        
        # logical size, confirm byteswapping
        r               = self.f.getRecords()[0]
        self.logSize    = len(r)

        if(self.debug):
            print '\nfrom file'
            print filename
            print 'logical size is :', self.logSize

        if(self.debug):
            print 'record is       : ', r
            print 'byteswapping    : ', self.byteSwap

        if(not self.logSize == 4 ):
            print "\nlogicalSize is    : ", self.logSize
            print "not 4, weirdness may happen.\n"
           
        self.offset += self.logSize

        if(self.debug):
            print 'got logical size: ', self.logSize 
            print 'got offset      : ', self.offset        

        if(self.debug):
            print
            print '*'*70
            print '\n    Sizes are: int: %d, real_kind: %d, logical: %d'%(self.intSize, self.doubleSize, self.logSize)
            print 'byteswapping for F95Binary: %d, byteswapping for regular data: %d\n'%(self.byteSwapfp, self.byteSwap)
            print '*'*70
            print
            
        self.Init = 1

    def __reopen(self):
        # reopen the file and reset the offset
        self.close();
        self.f.open(self.filename, 'r', self.byteSwapfp)
        self.offset = 0
        
    def close(self):
        # reopen the file and reset the offset
        if(self.f  != None):
            self.f.close();
        self.offset = 0
        if(self.debug):
            print 'file closed '
            print self.filename

    def getValues(self, typecode, shape, offset=-1):
        # synonym for self.get
        return(self.get(typecode, shape, offset))

    def get(self, variable=None, type=None, shape=None, offset=-1):
        
        if(self.debug):
            print 'type is            :',type
        if(self.debug):
            print 'shape is           :', shape, len(shape)
            
        if(not self.Init ):
            return(None)

        if(offset < 0):
            myoffset = self.offset
        else:
            myoffset = offset

        t = Numeric.array(0,type,copy=1)

        if(myoffset < self.offset):
            self.__reopen();

        if(self.debug):
            print 'discarding: ', myoffset - self.offset, 'bytes, (', myoffset, ',', self.offset,').'

        self.f.discardBytes(myoffset-self.offset)
        
        # count the number of bytes needed
        iCount = 1
        for i in range(0,len(shape)):
            iCount = iCount * shape[i]
            
        try:
            iBytes = iCount*t.itemsize
        except:
            iBytes = iCount*t.itemsize()
            
        if(self.debug):
            print 'iBytes gotten      :',iBytes
            
        self.rr = self.f.getBytes(iBytes)
        
        if(self.debug):
            print 'record from file is:___%s___'%self.rr, type, iCount, shape
            print Numeric.fromstring(self.rr, type, iCount)
            
        t = Numeric.reshape(Numeric.fromstring(self.rr, type, iCount), shape)
        
        if(self.debug):
            print 'Array gotten       :', [t]
            
        if(self.byteSwap):
            try:
                t = t.byteswap()
            except:
                t = t.byteswapped()

        self.offset = myoffset+iBytes

        return([t])
    
    def getInfo(self):
        """
        Returns information on records in file
        """
        if(self.f == None):
            return None

        # Note that we ignore the first three records for truchas files
        return self.f.getInfo()[3:]
    
    def __del__(self):
        
        self.close()
        self.f      = None
        self.offset = 0

             
class F95Binary:
    """
    F95Binary is a class for reading from and (eventually) writing
    into FORTRAN binary files.
    
    Author: Sriram Swaminarayan (sriram@lanl.gov)

    Since I use the array module, you are only allowed the following
    types:
    
        Type code   C Type             Minimum size in bytes 
        'c'         character          1 
        'b'         signed integer     1 
        'B'         unsigned integer   1 
        'u'         Unicode character  2 
        'h'         signed integer     2 
        'H'         unsigned integer   2 
        'i'         signed integer     4 
        'I'         unsigned integer   4 
        'f'         floating point     4 
        'd'         floating point     8 

    Note that gets are cross-record gets.  In the future I shall
    implement a strict mode where only one record is parsed at a time,
    but you can fake it currently by doing a self.getRecords(1),
    saving the result and processing your parameters within a single
    record.

    **Strict mode is currently unimplemented**
    
    """
    # Private type conversion dictionaries
    __c         = {1:'c', 2:'u'}
    __f         = {4:'f', 8:'d'}
    __iSigned   = {1:'b', 2:'h', 4:'i'}
    __iUnsigned = {1:'B', 2:'H', 4:'I'}

    def __init__(self, f95HeaderSize=4, debug=0):
        
        import array

        self.fp            = None
        self.Init          = 0
        self.mode          = 'r'
        self.name          = 'unknown'
        self.debug         = debug
        self.buffer        = array.array('c')
        self.byteSwap      = 0
        self.info          = []
        self.bsize         = array.array('i')
        self.f95HeaderSize = f95HeaderSize

    def __del__(self):

        if(self.Init):
            self.close()

    def getArray(self, type='', count=1):
        """

        return an array of type with count items read in from the
        file associated with the type.

        Allowed types are:
        Type code   C Type             Minimum size in bytes 
        'c'         character          1 
        'b'         signed integer     1 
        'B'         unsigned integer   1 
        'u'         Unicode character  2 
        'h'         signed integer     2 
        'H'         unsigned integer   2 
        'i'         signed integer     4
        'I'         unsigned integer   4 
        'f'         floating point     4 
        'd'         floating point     8 

        Default count is 1.
        
        """
        if ( not (type in ('c','b','B','u','h','H','i','I','f','d'))):
            print "Type is wrong.  See 'help(fortransupport)' for allowed types"
            return(None)

        retval = array.array(type)

        try:
            retval.fromstring(self.getBytes(count, retval.itemsize))
        except:
            retval.fromstring(self.getBytes(count, retval.itemsize()))

        if(self.byteSwap):
            try:
                retval = retval.byteswap()
            except:
                retval = retval.byteswapped()
        return(retval)

    def getInts(self, count=1, size=4, signed=1):
        """
        This subroutine gets count ints of size size from the file.
        The variable signed controls the conversion for
        signed/unsigned integers.

        Byte swapping is controlled by the value set for the F95Binary
        type when you opened the file.
        
        Note that gets are cross-record gets unless Strict() is specified.
        """
        retval = None
        
        if not self.Init:
            print "Must call open() first"
            
        elif not (size==1 or size==2 or size==4):
            print "size must be one of 1, 2, or 4 for ints"
            
        else:
            if(signed > 0):
                retval=array.array(self.__iSigned[size])
            else:
                retval=array.array(self.__iUnsigned[size])
                
            retval.fromstring(self.getBytes(count, size))

            if(self.byteSwap):
                try:
                    retval = retval.byteswap()
                except:
                    retval = retval.byteswapped()

        return(retval)
            
    def getReals(self, count=1, size=8):
        """
        This subroutine gets count reals of size size from the file.
        The default is double precision or size=8.  You can send in
        size=4 for single precision reals

        Byte swapping is controlled by the value set for the F95Binary
        type when you opened the file.
        
        Note that gets are cross-record gets unless Strict() is specified.
        """
        
        retval = None
        
        if not self.Init:
            print "Must call open() first"
            
        elif not (size==8 or size==4):
            print "size must be one of 4 or 8 for reals"
            
        else:
            retval = array.array(self.__f[size])

            retval.fromstring(self.getBytes(count, size))
            
            if(self.byteSwap):
                try:
                    retval = retval.byteswap()
                except:
                    retval = retval.byteswapped()

        return(retval)
            
    def getBytes(self, count=1, size=1):
        """
        This subroutine gets count bytes of size 'size' from the file.
        Note that gets are cross-record gets unless Strict() is specified.
        """

        import array
        retval = None

        if not self.Init:
            print "Must call open() first before"
        
        iRet     = 0
        iRetSize = count*size
        retval   = array.array('c')
        
        if not len(self.buffer) :
            self.buffer = self.getRecords(1)[0]
            
        while(iRet < iRetSize ):
            if(len(self.buffer) >= iRetSize-iRet ):
                retval      = retval+self.buffer[0:iRetSize]
                tmp         = self.buffer[iRetSize:]
                self.buffer = tmp
                iRet        = iRetSize
            else:
                retval      = retval+self.buffer
                iRet        = iRet+len(self.buffer)
                self.buffer = self.getRecords(1)[0]

        return(retval)


    def discardBytes(self, count=1, size=1):

        self.discardBytesOld(count,size)

        #self.discardBytesNew(count,size)

    def discardBytesOld(self, count=1, size=1):
        """
        This subroutine discards count bytes of size 'size' from the file.
        Note that discards are cross-record discards unless Strict() is specified.
        """

        if ( not self.Init ):
            print "Must call open() first before discarding bytes"
            return

        iDiscard     = 0
        iDiscardSize = count*size
        
        if not len(self.buffer):
            self.buffer = self.getRecords(1)[0]
            
        while iDiscard < iDiscardSize:
            
            if len(self.buffer) >= iDiscardSize-iDiscard:
                tmp         = self.buffer[iDiscardSize:]
                self.buffer = tmp
                iDiscard    = iDiscardSize
            else:
                iDiscard   += len(self.buffer)
                self.buffer =self.getRecords(1)[0]
        return

    """
    def discardBytesNew(self, count=1, size=1):
        #This subroutine discards count bytes of size size from the file.
        #Note that discards are cross-record discards unless Strict() is specified.
        import array
        if ( not self.Init ):
            print "Must call open() first before discarding bytes"
            return(None)
        iDiscard = 0
        iDiscardSize=count*size
        if(len(self.buffer) >= iDiscardSize-iDiscard ):
            # within the same record
            self.buffer=self.buffer[iDiscardSize:]
            return

        while(iDiscard < iDiscardSize ):
            isz = self.nextRecordSize()
            if(isz < 0): break            
            if(isz < iDiscardSize-iDiscard ):
                self.skipRecord()
            iDiscard=iDiscard+isz
        self.buffer=(self.getRecords(1)[0])[iDiscardSize-iDiscard:]
        return        
    """
    
    def open(self, name='unknown', mode='r', byteSwap=0):

        """
        name: name of file to be opened
        mode = 'r' or 'w' and determines whether the file is to be
               opened for reading or writing.  
        byteSwap  is a logical parameter which determines whether
               data will be byteswapped or not.
        """

        import array
        
        if ( mode != 'r' and mode != 'w'):
            print "Sorry, mode must be either 'r' or 'w'"
            return

        self.name = name
        self.mode = mode+'b'

        if byteSwap:
            self.byteSwap = 1
        else:
            self.byteSwap = 0

        self.fp = open(self.name, self.mode)
        # useful for writes
        if self.mode == 'wb':
            self.buffer = ''
        
        self.recl   = array.array('i',[1])
        self.Init   = 1
        return

    def close(self):
        """
        Close the file pointer opened by above open call.
        """
        if self.fp != None and self.mode=='wb':

            if self.buffer:
                self.write(endRecord=1)
                
            self.fp.close()
            self.__init__()

    def write(self, data=None, endRecord=0, debug=0):
        """
        Writes a strings and numeric arrays to file.
        if endRecord is true, the record is flushed to disk
        """
        if self.fp == None:
            print '\nUninitialized / unopened F95Binary file\n'
            return
        
        import array
        
        if self.mode != 'wb':
            print
            print '\nSorry, you can only write to a binary file using this routine.'
            print
            return

        if data != None:

            if type(data) == '':
                self.buffer += data
                
            else:
                if type(data) != type(self.recl):
                    # try to convert a numeric array
                    try:
                        # Numpy/Numeric kludge. Array does not know "S"
                        # which is "ASCII string" versus "character"
                        try:
                            tc = data.dtype.char
                            if (tc == 'S'):
                                tc = 'c'
                            d = array.array(tc, data.tostring())
                        except:
                            d = array.array(data.typecode(),data.tostring())
                    except:
                        print
                        print '\nSorry, you can only send arrays to the write() method'
                        return
                else:
                    d = data

                if(self.byteSwap):
                    try:
                        d = d.byteswap()
                    except:
                        d = d.byteswapped()
                self.buffer += d.tostring()

                # Superstition demands that I convert d back
                if(self.byteSwap):
                    try:
                        d = d.byteswap()
                    except:
                        d = d.byteswapped()
                # free data
                d = None
        #end if(data)
        
        if(endRecord):
            self.recl = array.array('i',[len(self.buffer)])
            if(self.byteSwap):
                try:
                    self.recl.byteswap()
                except:
                    self.recl.byteswapped()

            self.fp.write(self.recl.tostring())
            self.fp.write(self.buffer)
            self.fp.write(self.recl.tostring())
            self.buffer = ''

        return

    def getRecords(self,count=1):
        """
        Get count records from the file.  The length of the records is
        not fixed unless the file is composed of fixed length records.

        Typically you should not be using this for your data
        extraction, You should be using the getInts() , getFloats(),
        getBytes(), routines for data extraction.
        
        """

        records = []
        for x in range(0,count):
            records += [self.__getSingleRecord()]

        if(self.debug):
            print 'getRecords returning          :',records

        return(records)
    
    def __getSingleRecord(self):
        """
        This is a private subroutine.
        
        This subroutine will get the next record from file.
        Note that we do not byteSwap the result.  Only the size of the
        current record is byteSwapped.
        """
        import array
        bSize  = array.array('i')
        bSize2 = array.array('i')

        bSize.fromfile(self.fp,1) #self.f95HeaderSize)

        if(self.byteSwap):
            try:
                bSize = bSize.byteswap()
            except:
                bSize = bSize.byteswappped()

        if(self.debug):
            print '\n__getSingleRecord recordsize: ',bSize, bSize[0]

        record = array.array('c')
        record.fromfile(self.fp,bSize[0])
        bSize2.fromfile(self.fp,1) #self.f95HeaderSize)

        if(self.byteSwap):
            try:
                bSize2.byteswap()
            except:
                bSize2.byteswapped()

        if(self.debug):
            print '__getSingleRecord returning   :',record
            print 'recordsize                    :',bSize[0],bSize2[0]

        return(record)

    def nextRecordSize(self, size=1,savepos=1):
        """
        Returns the size of the next record.
        File position is left untouched.

        Returns -1 on failure
        """
        if self.fp  == None:
            return -1

        try:
            if savepos:
                pos = self.fp.tell()

            self.bsize.fromfile(self.fp,1) #self.f95HeaderSize)

            if(self.byteSwap):
                try:
                    self.bsize.byteswap()
                except:
                    self.bsize.byteswapped()
            isz = self.bsize[0]
            
            if(self.debug):
                print 'nextRecordSize            :', isz

            if(savepos):
                self.fp.seek(-4,1)
            self.bsize.pop()
            
            return isz
        except:
            return -1


    def skipRecord(self,count=1,bytes=0):
        """
        skips count records
        returns 0 on success
        """
        import array

        if self.fp  == None:
            return -1

        try:
            if bytes >0 :
                self.fp.seek(bytes,1)
            else:
                self.bsize.fromfile(self.fp,1) #self.f95HeaderSize)
                if(self.byteSwap):
                    try:
                        self.bsize.byteswap()
                    except:
                        self.bsize.byteswapped()

                if(self.debug):
                    print 'skipRecord              :',self.bsize, self.bsize[0]

                #self.fp.seek(self.bsize[0]+4,1)
                self.fp.seek(self.bsize[0]+self.f95HeaderSize,1)
                self.bsize.pop()
            return 0
        except:
            return -1

    def getInfoNew(self):
        """
        Returns information on records in file
        """
        import array
        import sys

        if len(self.info):
            return self.info
        
        if self.fp == None:
            return array.array('i')

        fpos = self.fp.tell() # save original file position
        self.fp.seek(0)       # goto beginning of file

        self.info = array.array('i')
        try:
            pos = 0
            while (1):
                isz = self.nextRecordSize(savepos=0)
                if isz < 0 :
                    break

                self.info.append(isz)
                self.skipRecord(bytes=isz+4)

                if pos%5000 == 0:
                    sys.stdout.flush()
                
                pos += 1

        except:
            if(self.debug):
                print "file"
                print self.name
                print "info len is                  : ",pos

        self.fp.seek(fpos) # goto original file location

        return self.info

    def getInfo(self):
        """
        Returns information on records in file
        """
        import sys, array

        if len(self.info):
            return self.info

        if self.fp == None:
            return array.array('i')

        self.info = array.array('i')
        tmpInfo   = array.array('i')
        fpos      = self.fp.tell() # save original file position
        self.fp.seek(0)            # goto beginning of file

        try:
            pos = 0
            while(1):
                tmpInfo.fromfile(self.fp,1) #self.f95HeaderSize)

                if(self.byteSwap):
                    try:
                        tmpInfo = tmpInfo.byteswap()
                    except:
                        tmpInfo = tmpInfo.byteswapped()

                #self.fp.seek(tmpInfo[0]+4,1)
                self.fp.seek(tmpInfo[0]+self.f95HeaderSize,1)
                self.info.append(tmpInfo[0])

                if( pos%5000 == 0 ):
                    sys.stdout.flush()
                pos += 1

                if self.debug:
                    print "tmpinfo                    :",tmpInfo[0]
                    print "info is                    :",self.info
                    print "info len is                :",pos

                tmpInfo.pop()

        except:
            if(self.debug):
                print "file"
                print self.name
                print "info len is                  :",pos

        self.fp.seek(fpos) # goto original file location

        return self.info

    
# End class F95Binary

if __name__ == '__main__':

    import subprocess #popen2
    
    'for testing this component'

    from PYTHONutils import uTestOpts
    binfile    = 'header.bin'
    #opts      = uTestOpts('fd',defaults={'f':binfile})
    dir        = '../../../scripts/test_TBrookParse/samples/phasechange_mixed_output/' 
    dfltfl     = dir + 'phasechange_mixed.DefaultMesh.00001.bin'
    opts       = uTestOpts('fd', defaults={'f':dfltfl})
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    try:
        #P = popen2.Popen4('g77 -o testheader createheader.f')
        #P.tochild.close()
	cmd = ['g77','-o','testheader','createheader.f']
	rtn = subprocess.call(cmd)
        #P = popen2.Popen4('testheader')
        #P.tochild.close()
        rtn = subprocess.call('testheader')
    except:
        print "---> Unable to create %s file " %(binfile)
        print "\nTry different binary file as input\n"
    try:
        file = T_Binary_File(opt.f,debug=1)
        file.getInfo()
        file.get(type='c',shape=[1],offset=16)
        file.get(type='i',shape=[1],offset=19)
        file.get(type='i',shape=[100],offset=23)
    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise

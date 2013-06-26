"""

 dataSource

 -----------------------------------------------------------------------------
  Purpose:
  
    Instantiates a dataSource object given a file name and type.

    It provides anonymous get and getInfo methodss to the outside
    world to get specified data from the file. The actual methods used
    depend on the <filetype> argument used in object creation.

    The class is instantiated without opening the data file.  Only
    when the first get is called will the file be opened.

    If the file does not exist, a warning message will be printed
    during instantiation.

  Expected Usage:
      This class will be provided to the variables in the srcPointer
      array as the value for the 'tfile' key.  The class will be
      instantiated in the XMLgetTimestep() just before the 'tfile' key
      is set in XMLgetTimestep.self.pointer{}.  When the first call is
      made to the get() function, the file will be opened.  It is
      considered an error if the file does not exist and a message will
      be printed to the stderr and None will be returned for the data.
      Calling close() will close the file pointer.  

  Public Interface(s):
     ds = dataSource(filename,filetype)
     ds.getInfo(*a)
     ds.get(*a)
     ds.close()

  Contains:
    class dataSource
        __init__(self, filename, filetype, debug=0)
        __getDefault(self, *a)
        __getXML(self,*a)
        __getInfoXML(self,*a)
        __getF95b(self,name, type, shape, offset)
        __getInfoF95b(self,name, type, shape, offset)
        close(self)

    ---> NO Unit Test Block
  
  Author(s): 
 -----------------------------------------------------------------------------
"""

import os, sys
class dataSource:
    """
    Expected Usage:
      This class will be provided to the variables in the srcPointer
      array as the value for the 'tfile' key.  The class will be
      instantiated in the XMLgetTimestep() just before the 'tfile' key
      is set in XMLgetTimestep.self.pointer{}.  When the first call is
      made to the get() function, the file will be opened.  It is
      considered an error if the file does not exist and a message will
      be printed to the stderr and None will be returned for the data.
      Calling close() will close the file pointer.  
    """
    
    def __init__(self, filename, filetype, debug=0):
        """
          filename: name of file where data sits
          filetype: type of file, i.e. 'xml', 'XML','f95b',
                    'fortranbinary', ('hdf', 'udm'), etc 
          debug: flag that controls printing of debug information
        Note that no default is provided for filetype.
        If file does not exist, a warning message will be printed to
        screen once.
        """
        from FORTRANutils import fortransupport
        self.filename = filename
        self.filetype = filetype
        self.file     = None
        self.debug    = debug

        if(self.debug): print '  __dataSource: file=', filename, ' type = ', filetype

        if ( os.path.exists(filename) ):
            self.status = 0
        else:
            sys.stderr.write('\n\n  __dataSource WARNING: '+filename+' does not exist\n\n')
            self.status = -1

        if ( (self.filetype == 'f95b') or (self.filetype == 'fortranbinary')):
            self.get = self.__getF95b
        elif ( (self.filetype == 'xml') or (self.filetype == 'XML')):
            self.get = self.__getXML
        else:
            self.get = self.__getDefault
            sys.stderr.write('\n\n  __dataSource WARNING: Unknown type: '+filetype+'\n\n')
            self.status = -1

        if ( (self.filetype == 'f95b') or (self.filetype == 'fortranbinary')):
            self.getInfo = self.__getInfoF95b
        elif ( (self.filetype == 'xml') or (self.filetype == 'XML')):
            self.getInfo = self.__getInfoXML
        else:
            self.getInfo = self.__getDefault
            sys.stderr.write('\n\n  __dataSource WARNING: Unknown type: '+filetype+'\n\n')
            self.status = -1

        return

    def __getDefault(self, *a):
        """
        This is the default method.  It does nothing and returns None
        for values.  This happens only if the file did not exist on
        initialization or a bad type was given.
        
        """
        if(self.debug):
            print '  __dataSource getDefault(): status is: ', self.status
        return None

    def __getXML(self,*a):
        print """
        getXML is not implemented yet
        """
        return None

    def __getInfoXML(self,*a):
        
        print """
        getInfoXML is not implemented yet
        """
        return None
    
    def __getF95b(self,name, type, shape, offset):

        if(self.debug):
            self.prefix = '  __dataSource.getF95b: '
            print self.prefix, 'args are: ',name,type,shape,offset
            print self.prefix, 'status is: ', self.status
        from FORTRANutils import fortransupport


        if (self.status < 0 ): return None
        if(self.file == None):
            self.file = fortransupport.T_Binary_File(self.filename)
            if(self.file):
                self.status =  1
            else:
                self.status = -1
            
        if(self.file):
            return(self.file.get(name,type,shape,offset))
        else:
            return None

    def __getInfoF95b(self,name, type, shape, offset):

        if(self.debug):
            self.prefix = '  __dataSource.getInfoF95b: '
            print self.prefix, 'args are: ',name,type,shape,offset
            print self.prefix, 'status is: ', self.status
        from FORTRANutils import fortransupport


        if (self.status < 0 ): return None
        if(self.file == None):
            self.file = fortransupport.T_Binary_File(self.filename)
            if(self.file):
                self.status =  1
            else:
                self.status = -1
            
        if(self.file):
            return(self.file.getInfo())
        else:
            return None
        

    def close(self):
        if (self.status == 1 ):
            if(self.debug): print 'closing ',self.filename
            self.file = None
        else:
            if(self.debug): print 'file '+self.filename+'is already closed'

        self.status = 0
        return

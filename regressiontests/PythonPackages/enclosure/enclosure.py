"""
This script will allow users to open viewfactor files
and return data in a viewfactor data type
"""
import sys


class enclosure:
    init    = 0
    name    = 'none'
    AType   = 0
    type    = 'b'
    comment = 'no comment'
    cell    = []
    face    = []
    counts  = []
    offset  = []
    entryid = []
    entries = []

    def info(self, v=0):
        """
        Print information about the enclosure.
        How much information is controlled by value of v.
        Higher is more!
        """
        # Basic information is always printed
        print ''
        print '  File:    ',self.name
        print '  Type:    ',self.type
        print '  Comment: ',self.comment
        print '  Nfaces:  ',self.nfaces_g
        print '  NEntries:',self.rcounts_total

        # if v >= 1 then we print cell/face ids and counts
        if ( v >= 1):
            print     '         CELL   FACE    ENTRIES',
            if(v>=2):
                print '     VF Sum',
            if(v>=3):
                print '     VFs',
            print ' '
            for i in range(self.nfaces_g):
                print '  %10d %6d %10d'%(self.cell[i], self.face[i], self.counts[i]),
                if(v>=2):
                    print '  %12.9f'%self.querysum(i),
                if (v >= 3):
                    for x in self.entries[self.offset[i]:self.offset[i]+self.counts[i]]:
                        print x,
                print ' '
        print ''

        
    def __init__(self,name='none',type='b',debug=0):
        """
        Open a file and read the viewfactor data.
        If there is an error this will return None
        """
        if(debug):
            print 'Trying file = ',name, ' of type ', type

        if ( type == 'a' or type == 'ascii'):
            type = 'a'
            AType = 1
        elif (type == 'b' or type == 'binary'):
            type = 'b'
            AType = 0
        else:
            print ''
            print 'Wrong type in enclosure()'
            print ''
            return None

        if ( AType) :
            # We have to read data from ascii file
            try:
                fileInput  = open(name, 'r')         # Open file
            except:
                print ''
                print 'error opening file ',name
                print ''
                return None
         
            self.comment = fileInput.readline()
            if(debug):
                print '\n\ncomment is: '
                print self.comment.rstrip()
                print ' '

                print 'getting file info'

            [self.nfaces_g, self.rcounts_total]  = self.getInts(fileInput.readline())
            if(debug):
                print 'nfaces_g = ',self.nfaces_g
                print 'rcounts_total = ', self.rcounts_total

            if(debug): print '  cell'
            self.cell = self.getInts(fileInput.readline()) 

            if(debug): print '  face'
            self.face = self.getInts(fileInput.readline()) 

            if(debug): print '  count'
            self.counts = self.getInts(fileInput.readline()) 

            if(debug): print '  offset'
            itmp = 0
            for n in self.counts:
                self.offset = self.offset + [itmp]
                itmp = itmp + n

            if(debug): print '  entryID'
            self.entryid = self.getInts(fileInput.readline()) 

            if(debug): print '  entry'
            self.entries = self.getFloats(fileInput.readline()) 

            if(debug):
                print 'done getting file info'
                print
                print 'total faces = ',nfaces_g
                print 'total entries = ', rcounts_total
            fileInput.close()
            fileInput=None

        else:
            # This is a binary file.
            # Use the fortransupport library
            try:
                from fortransupport import F95Binary
            except:
                print '  '
                print '  Unable to import fortransupport module.'
                print '  Ensure that PYTHONPATH environment variable '
                print '  includes truchas/tools/PythonPackages/TBrookSupport'
                print '  '
                raise 'unable to import fortransupport error'
            f = F95Binary()
            f.open(name)
            if ( not f ):
                print ''
                print 'error opening binary file ',name
                print ''
                return None
            self.comment = ((f.getRecords(1))[0]).tostring().rstrip()
            if(debug): print self.comment


            [self.nfaces_g, self.rcounts_total]  = f.getInts(2)
            if(debug):
                print 'nfaces_g = ',self.nfaces_g
                print 'rcounts_total = ', self.rcounts_total

            if(debug): print '  cell'
            self.cell = f.getInts(self.nfaces_g) 

            if(debug): print '  face'
            self.face = f.getInts(self.nfaces_g) 

            if(debug): print '  count'
            self.counts = f.getInts(self.nfaces_g) 

            if(debug): print '  offset'
            itmp = 0
            i = 0
            self.offset=range(self.rcounts_total)
            for n in self.counts:
                if ( i > 99030):print i, itmp
                self.offset[i] = itmp
                itmp = itmp + n
                i = i + 1
                

            if(debug): print '  entryID'
            self.entryid = f.getInts(self.rcounts_total)

            if(debug): print '  entry'
            self.entries = f.getReals(self.rcounts_total) 

            if(debug):
                print 'done getting file info'
                print
                print 'total faces = ',self.nfaces_g
                print 'total entries = ', self.rcounts_total
            f.close()
            f=None

   

        self.init = 1

    def querysum(self,i):
        """
        Returns sum of VFs for face entry i
        """
        
        if (i>=len(self.cell) or i < 0):
            print '\n  ***************\n'
            print '  Invalid request.  '
            print '  ***************\n'
            return None
        vfsum = 0.0
        for x in self.entries[self.offset[i]:self.offset[i]+self.counts[i]]:
            vfsum = vfsum + x
        return(vfsum)

    
    def getInts(self,t):
        """
        Split a string and return a list of integers
        """
        return([int(x) for x in t.split()])
        
    def getFloats(self,t):
        """
        Split a string and return a list of floats
        """
        return([float(x) for x in t.split()])

def pArray(a, label, format):
        """
        print array a preceeded by label using format
        """
        print label,
        for x in a: print format%x,
        print 

                        

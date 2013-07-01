

Estring = 'System Import Error'

try:
    Estring = 'os/sys Import Error'
    import os
    import sys

    Estring = 'Numeric Import Error'
    try:
        import numpy.oldnumeric as Numeric
    except ImportError:
        import Numeric
    except:
        raise
    
    Estring = 'Path Setting Error'
    pathdir = os.path.abspath(os.path.dirname(__file__))
    if(pathdir == '/'):
        pathdir = './'
    pathdir = pathdir + '/..'                             
    sys.path   = sys.path + [pathdir +'/..']
    sys.path   = sys.path + [pathdir +'/../TBrookParser']
    sys.path   = sys.path + [pathdir +'/../enclosure']
    sys.path   = sys.path + [pathdir +'/POSTPROCESSORutils'] 


    Estring = 'XMLutils import error'
    from XMLutils               import BrookSupport

    Estring = 'POSTPROCESSORutils import error'
    from POSTPROCESSORutils        import createFileObject,usubs
    from getPostProcessorObject    import getPostProcessorObject
    from GMVwriteStorageObject     import GMVwriteStorageObject
    from ENSIGHTwriteStorageObject import ENSIGHTwriteStorageObject

except:
    print "\n"
    print "\n"
    print "\n"
    print "  ERROR: (in file "+__file__+")"
    print "    " + Estring
    print "    The pathmagic script will now exit"
    print
    print
    sys.exit(1)


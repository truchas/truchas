"""
get_aborted_residuals.py
Written by Sharen Cummins

This module reads in a TBrook file, and selects any linear and non-linear
residual information written out when the Truchas simulation is aborted. It plots the
results in a binary GMV file.

GMV files for the non-linear residual are of the form
    gmv_simulation_name_nonlinres.iteration
and are placed in the same directory as the TBrook file.

GMV files for the linear residual are of the form
    gmv_simulation_name_linres.iteration
and are placed in the same directory as the TBrook file.

The script is executed by entering:

python /path_to_script/get_aborted_residuals.py TBrook_file_name

This script uses the scriptUtils module written by Sriram Swaminarayan
to access the TBrook file data.

"""

import os, sys, getopt
try:
    # Do path magic
    truchasdir     = os.path.abspath(os.path.dirname(sys.argv[0]))
    if(len(truchasdir)):
        truchasdir = truchasdir+'/../'
    else:
        truchasdir = '../'
    sys.path.append(truchasdir)
    from PythonPackages import *
    sys.path.pop()
except:
    print """
    
    Unable to import standard modules Please ensure that
    <truchasdir>/tools is in your PYTHONPATH variable.

    It could also be that you don't have Numeric installed.

    """
    sys.exit(1)


argv        = sys.argv[1:]
opts, rest  = getopt.getopt(argv, 'hf:e:')
if (not len(rest)):
     print
     print ' No TBrook.xml file specified on command line'
     print ' Command line should be of the form '
     print 
     print '      python path_to_script/get_aborted_residuals.py file.TBrook.xml ' 
     print 
     sys.exit()
    
files = rest[0]

# load file and primary mesh
f = suFile(files)                  # try f.help() to get info on the suFile structure
m = suMesh(f.meshes[0],f)          # try m.help() to get info on the suMesh structure

#print files

input = usubs.input()

#plot any nonlinear residuals..

for nlr in f.storage.aborts.nlrlist:
    f.storage.getValues(nlr.vlist)
    idx = f.storage.aborts.nlrlist.index(nlr)
    
    # Write GMV file(s)
    filename  = files[0:len(files)-11]                   # Strip the .TBrook.xml suffix from the file name
    fileout   ='%s_%s_gmv.%06d'%(filename,'nonlinres',idx)   # The output file name
    outputfmt = 'binary'                                 # The output format
    vars      = nlr.vlist                                # Variables to be put into GMV file
    seq_no    = 1                                        # The sequence number
    flags     = {}

    # Write our variables to a binary GMV file
    GMVwriteStorageObject(fileout, outputfmt, f.meshes[0],
                          0.0, 0, vars, seq_no, flags, debug=0)

#plot any linear residuals..

for lr in f.storage.aborts.lrlist:
    f.storage.getValues(lr.vlist)
    idx = f.storage.aborts.lrlist.index(lr)
    
    # Write GMV file(s)
    filename  = files[0:len(files)-11]                   # Strip the .TBrook.xml suffix from the file name
    fileout   ='%s_%s_gmv.%06d'%(filename,'linres',idx)      # The output file name
    outputfmt = 'binary'                                 # The output format
    vars      = lr.vlist                                 # Variables to be put into GMV file
    seq_no    = 1                                        # The sequence number
    flags     = {}

    # Write our variables to a binary GMV file
    GMVwriteStorageObject(fileout, outputfmt, f.meshes[0],
                          0.0, 0, vars, seq_no, flags, debug=0)            


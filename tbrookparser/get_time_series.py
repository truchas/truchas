"""
get_time_series.py
Written by Sharen Cummins

This module writes out a given variable and its index for a given set of cycles [Cyclea, Cycleb].
If the variable is a scalar, a simple text file is produced of the form

        cycle   |   t     |   VariableName (at Index = I)
       -------------------------------------------------   
          5     |   1.1   |     V_I
          6     |   1.2   |     ...

If the variable is a vector, the text file will be of the form

       cycle    |      t     |   VariableName (at Index = I)
      ------------- ----------------------------------------   
         5      |    1.1   |     X_I  Y_I  Z_I
         6      |    1.2   |     .....

If the variable does not live on a meshspace (such as sensitivity function values) the file is of the form

       cycle    |      t     |   VariableName 
      ------------- ----------------------------------------   
         5      |    1.1   |     X_I  Y_I  Z_I...(i.e number of elements = shape of the variable)
         6      |    1.2   | 

The resulting file has a default name of the form

         VariableName_Index_Cyclea_to_Cycleb.txt

and is placed in the same directory as the TBrook file.

The script is executed by entering:

python /path_to_script/get_time_series.py.

This script uses the scriptUtils module written by Sriram Swaminarayan
to access the TBrook file data.

"""
import os, sys, string

try:
    # Do path magic
    truchasdir = os.path.abspath(os.path.dirname(sys.argv[0]))
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

def_file     = os.path.abspath(truchasdir) + '/scripts/test_TBrookParse/samples/phasechange_mixed_output/phasechange_mixed.TBrook.xml'
inFile       = def_file
#check if there is a TBrook.xml file in the current working directory
for file in os.listdir(os.getcwd()):
    if file[len(file)-11:len(file)] == '.TBrook.xml':
        inFile  = file
def_exofile     = ''

from POSTPROCESSORutils     import createFileObject,usubs
from getPostProcessorObject import getPostProcessorObject

# Set up our input buffer
input      = usubs.input('')

# Now create a file object
X          = createFileObject()

ppobject   = getPostProcessorObject(def_file,def_exofile,input,debug=0)
iExit      = None
try:
    inFile = ppobject.input.usetc('Name of file to be read:', inFile)
    if inFile not in ppobject.files:
        X  = ppobject.load(inFile)
except:
    print
    print 'Error in loading file. Try loading a different file'
    print
    sys.exit(1)

while (iExit != 1):

    #Now get the variable and the index (if the variable lives on a meshspace)
    varformat    = 3
    this_storage = X.storage

    count           = 0
    varlist         = '\n'
    vlist           = []
    for var in this_storage.tlist[0].vlist:
        if var.rank > 0:
            s           = '%+30s' %(str(var.nickname))
            s           = s + ','
            varlist     = varlist + s
            vlist.append(var.nickname)
            count       = count + 1
            if (count%varformat == 0):
                varlist = varlist + '\n'
    varlist     = varlist[0:len(varlist)-1]
    prompt      = 'Time series for which variable? (%s)' %(varlist)
    varname     = ppobject.input.usetc(prompt,str(this_storage.tlist[0].vlist[1].nickname))

    if (len(varname)<1):
        print
        print 'No variable specified..try again'
        print
        while (len(varname)<1): 
            print
            print '*************************************************'
            print 'No variable specified... try again'
            print '*************************************************'
            print
            varname     = ppobject.input.usetc(prompt,str(this_storage.tlist[0].vlist[1].nickname))
    elif not (varname in vlist):
        while not (varname in vlist): 
            print
            print '*************************************************************'
            print 'This variable is not available in the XML output... try again'
            print '*************************************************************'
            print
            varname     = ppobject.input.usetc(prompt,str(this_storage.tlist[0].vlist[1].nickname))
    else:
        novar               = 1
        face                = 0
        for k in this_storage.tlist[0].vlist:
            if(k.name == varname or k.nickname == varname):
                thevar       = k
                novar        = novar*0
                if k.rank == 3:
                    print 'This variable is face-centred vector with %i components' %(k.shape[len(k.shape)-1])
                    facelist = range(k.shape[len(k.shape)-2]+1)
                    facelist = str(facelist[1:len(facelist)])
                    face     = ppobject.input.usetc('Specify which face you would like to analyze out of ' + facelist, '1')
                    face     = int(face)

    #Check the variable lives on a meshspace..if so inquire as to which index to analyze
    #otherwise just print out the variable

    nomeshspace = 0
    if thevar.meshspace == 'CELL':
        for meshspace in this_storage.mlist[0].mslist:
            if meshspace.name == 'CELL':
                if meshspace.size != thevar.shape[0]:
                    nomeshspace = 1

    if not nomeshspace:
        for meshspace in this_storage.mlist[0].mslist:
            if thevar.meshspace == meshspace.name:
                ninds     = int(meshspace.size)
                nindsh    = int(0.5*ninds)
                prompt    = 'Index to be analyzed in the time series [1 -> %s]' %(str(ninds))
                index     = ppobject.input.useti(prompt,str(nindsh))
    else:
        index = 'Not relevant'

    #Now inquire about which cycles to print out...

    lowcycle   = this_storage.tlist[0].cycle
    highcycle  = this_storage.tlist[len(this_storage.tlist)-1].cycle
        
    lower      = ppobject.input.useti('Low Cycle number?',lowcycle)
    upper      = ppobject.input.useti('High Cycle number?',highcycle)

    idx2       = []
    cnt        = 0
    for i in this_storage.timesteps.clist[thevar.name]:
        if (i <= upper and i >= lower):
            idx2 = idx2 + [cnt]
        cnt = cnt + 1

    #Header..

    s = '#' + '-' * 115 +'\n'
    s +='#' + '%+10s %+20s %+30s (Index = %s)' %('CYCLE', 'TIME', varname, str(index))
    s += '\n' + '#' + '-' * 115 + '\n'

    for i in idx2:
        
        for var in this_storage.tlist[i].vlist:
            if (var.name == varname or var.nickname == varname):
                thevar = var

        this_storage.getValues([thevar])

        tmp = ''
        if thevar.rank == 1:
            if nomeshspace :
                #e.g sensitivity variables or em variables - i.e variables that don't live on a meshspace
                for j in thevar.data:
                    tmp += '%20.3E' %(j)
            else:
                tmp    = '%20.3E' %(thevar.data[index-1])
        if thevar.rank == 2:
            for j in thevar.data[index-1]:
                tmp += '%20.3E ' %(j)
        if thevar.rank == 3:
            cnt = 0
            for j in thevar.data[index-1]:
                cnt += 1
                if cnt == face:
                    for f in j:
                        tmp += '%20.3E ' %(f)
    
        s += '%10d %20.3f %+30s '%(this_storage.tlist[i].cycle, this_storage.tlist[i].time, tmp)
        s += '\n'

    #Now inquire about the name of the txt output file 
    prompt      = 'Filename to output results to?'
    if nomeshspace:
        defname     = '%s_%s_to_%s.txt' %(thevar.nickname,str(lower),str(upper))
    else:
        defname     = '%s_%s_%s_to_%s.txt' %(thevar.nickname,str(index),str(lower),str(upper))
    fname       = ppobject.input.usetc(prompt, defname)
    fp          = open(fname,'w')
    print >> fp, s
    fp.close()

    prompt      = 'Quit now (q) or not (no)?'
    chckqt      = ppobject.input.usetc(prompt,'q')
    if chckqt == 'q':
        iExit = 1

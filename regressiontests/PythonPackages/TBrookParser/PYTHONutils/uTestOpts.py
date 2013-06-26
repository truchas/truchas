"""

 uTestOpts

 -----------------------------------------------------------------------------
  Purpose:

    This package provides a simple wrapper of OptionParser to give
    command-line argument support for unit testing.  The user passes
    a string of option characters (e.g. 'fd')

    Defaults are provided for 'f' (file: a TBrook.xml file) and 'd'
    (debug: False)

    The defaults, actions (how to store the info), and types
    (argument types) can also be set, but are defaulted to None,
    'store' and 'string', respectively

    This class provides a header method that writes out the python
    file being tested and the (TBrook) input file being used.

    It also provides a help option that gives hints to its use.

    Usage example:
       A = uTestOpts('fd')
       (o,a) = A.parse_args()
       print o.f, o.d

  Public Interface(s):

    opt = uTestOpts(flags)

  Contains:
    class uTestOpts(OptionParser)
        __init__(self, flags, defaults={}, actions={}, types={}, dir=None)
        header(self, unit, file)
        help(self)

    Unit Test Block

  Version:
    $Id: uTestOpts.py,v 1.2 2007/01/25 21:59:23 scummins Exp $

  Author(s): Larry Cox (ljcox@lanl.gov)
 -----------------------------------------------------------------------------
"""

from optparse import OptionParser

class uTestOpts(OptionParser):
    "Simple interface to OptParser for typical unit testing"

    def __init__(self, flags, defaults={}, actions={}, types={}, dir=None):

        "Initialize requested options setting defaults/actions (if any)"

        OptionParser.__init__(self)

        # default source directory and TBrook.xml file
        if dir == None:
            # look up three levels, then down for the samples directory
            # the argument 'dir' overrides this location
            self.__dir  = '../../../scripts/test_TBrookParse/samples/'
        else:
            self.__dir = dir

        self.__file = self.__dir + \
                      'phasechange_mixed_output/phasechange_mixed.TBrook.xml'

        for k in flags:
            (act, df, typ) = ('store', None, None)

            if actions.has_key(k):  act = actions[k]
            elif k == 'd':          act = 'store_true'

            if defaults.has_key(k): df  = defaults[k]
            else:
                if   k == 'f':      df = self.__file
                elif k == 'd':      df = False

            if types.has_key(k):    typ = types[k]

            if typ:
                self.add_option('-'+k, dest=k,
                                default=df, action=act, type=typ)
            else:
                self.add_option('-'+k, dest=k,
                                default=df, action=act)
                
    def header(self, unit, file=None):
        " print unit test header "
        import os
        print "*",os.path.basename(unit),"Component Test ..."
        if file:
            print "* Using input file:", os.path.basename(file)
            print "*   from directory:", os.path.dirname(file)
        print "*********************\n"

    def help(self):
       print  \
       """
       Simple interface to OptParser for typical unit testing
       Establishes one command-line option for each character provided.
       Usage Example:
        from PYTHONutils import uTestOpts
        flags = 'fd'
        p = uTestOpts('fd')
        (opt,arg) = p.parse_args()

        print 'File is', opt.f
        if opt.d: # if '-d' was provided
        ...
       Note that internal defaults are provided for 'f' and 'd' only
       """

if __name__=='__main__':

    p = uTestOpts('fdt',
                  actions  = {},
                  defaults = {'t': 1},
                  types    = {'t':'int'})
    (o,a) = p.parse_args()

    # Print test header with tested filename and input filename
    #   Note, must occur after call to p.parse_args()
    p.header(__file__,o.f)

    # List established options and their received or defaulted values
    print 'f:', o.f
    print 'd:', o.d
    print 't:', o.t


#!/usr/bin/python

####
# Copyright (C) 2000 Yianilos Laboratories
#    Author: Daniel Gindikin <dan@netrics.com>
#
# This source is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation, version 2.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You can retrieve a copy of the GNU Library General Public License
# from http://www.gnu.org/.  For a copy via US Mail, write to the
#
#     Free Software Foundation, Inc.
#     59 Temple Place - Suite 330,
#     Boston, MA  02111-1307
#     USA
#
####

"""
    This module provides functions to help you process the command line

    Common Definitions:
    ------------------
	Just some definitions to avoid any confusion.

	argument: every entry in sys.argv is an argument.

	option: an argument that is either prefixed by '-', or '--', or
	    has a value associated with it (i.e. "opt=val" or "opt:val").

	value: the value associated with an option.

	long option: an argument that is prefixed by '--'.

	short option: an argument that is prefixed by '-'. Short options
	    can be combined together into sets. For example, if sys.argv
	    contains the entry '-abcd', then all of the following calls
	    will return 1:
			cmdline.receivedOption('-a')
			cmdline.receivedOption('-d')
			cmdline.receivedOption('-bc')
			cmdline.receivedOption('-abcd')
	    Short options can have values associated with them just like
	    long options, except if they are in a set. If there is an
	    argument on the command line '-abcd=3', and you call
			cmdline.getIntegerValueOf('-a')
	    this module will complain and exit, since it would have no way
	    of knowing whether the 3 was associated with '-a' or with '-b'.
			cmdline.getIntegerValueOf('-abcd')
	    will return 3.
	    CAUTION: it is up to you to make sure there are no conflicts
	             between short options. For example expecting both
		     '-a' and '-abcd' is a mistake, since this module
		     will have no easy way of knowing whether an argument
		     '-abcd' on the command line is equivalent to the short
		     option '-abcd', or, for example, to the two short
		     option '-a' and '-bcd'.
	
	parameter: an argument which is neither an option, nor a value for
	    an option.


    The Interface:
    -------------
	cmdline.receivedOption(opt)
	    opt -- either an option, or a list of options, e.g. '-h --help'
	    return
		1 if the option is present, 0 otherwise

	cmdline.getIntegerValueOf(opt, [default])
	    opt -- either an option, or a list of options, e.g. '-p --port'
	    default -- the default value to return if there is no option opt
		on the command line, or if the option is there, but it does
		not have a value associated with it.
	    return
		the value passed to the option opt. If the value is of the
		wrong type ('--opt=3.3'), an error message will be printed
		and the program will exit. If there is no option opt on the
		command line, or if the option is there, but it does not have
		a value, and a default was passed in, the default will
		be returned. Otherwise an error message will be printed
		and the program will exit.

	cmdline.getFloatValueOf(opt, [default])
	cmdline.getStringValueOf(opt, [default])
	    same as cmdline.getIntegerValueOf(), but for different types

	cmdline.getRemaining([description of what you expect])
	    Get the first argument, that is not an option, that has not
	    been used so far. Useful if you are expecting some parameters
	    to be passed in without associated options. You should call
	    this after you are done processing all the options. If you pass
	    in the description of what you expect, and there are no more
	    unused arguments, an error message will be printed and the
	    program will exit. Otherwise None will be returned, and its up
	    to you to deal with it.

	cmdline.printWarnings()
	    This prints out warnings about unused arguments. You should call
	    this after you are done processing the command line. It will
	    return 0 if there were no warnings printed, 1 otherwise. If you
	    are a risk averse type individual, you should probably exit
	    if 1 is returned.

    Example Usage:
    -------------
	import cmdline

	def usage():
	    import sys
	    print 'Usage: %s --bias float [-p|--port [port]] [-t|--trace [-l|--lvl|--level [level]]] hostname' % sys.argv[0]
	    raise SystemExit

	if cmdline.receivedOption('-h --help'):
	    usage()

	# 5555 is the default value returned if there is no
	# option '-p' or '--port' on the command line
        # or if they are there, but they do not have a value
	# associated with them
	port = cmdline.getIntegerValueOf('-p --port', 5555)
	print 'the port is', port

	# if there is no option '--bias' on the command line
	# or if it does not have a float argument, the error
	# information will be printed, and the program will exit
	bias = cmdline.getFloatValueOf('--bias')
	print 'the bias is', bias

	if cmdline.receivedOption('-t --trace'):
	    print 'turning on trace'
	    lvl = cmdline.getIntegerValueOf('-l --lvl --level', 5)
	    print '\tlevel', lvl

	# hostname passed in without any flags
	hostname = cmdline.getRemaining('hostname')
	print 'The hostname is', hostname

	# print warnings about unused arguments and options
	cmdline.printWarnings()
"""

import sys
import re
from string import find, atoi, atof, split

__version__ = '1.0'

# some globals
# XXX mayhaps bite the bullet and have an internal class?
gOpts = []
gValues = []
gNotused = []
gCmdLine = []

def hasValue( arg ):
    if ('=' in arg) or (':' in arg):
	return 1
    else:
	return 0

def isOption( arg ):
    "see if arg is an option"
    if isLongOption(arg) or isShortOption(arg):
	return 1
    # it may have a value, if does, give it the benefit of the doubt
    # and consider it an option
    if hasValue(arg):
	return 1
    return 0

def isLongOption( arg ):
    "an option of the form '--bla' is considered a long option"
    if arg[0:2] == '--':
	return 1
    else:
	return 0

def isShortOption( arg ):
    "an option of the form '-bla' is considered a short option"
    if (arg[0] == '-') and (arg[1] <> '-'):
	return 1
    else:
	return 0

def processArg( arg ):
    """
    returns the tuple (option, value), for example
	processArg('--foo:bar') -> ('--foo', 'bar')
	processArg('--foo=bar') -> ('--foo', 'bar')
	processArg('--foo')	-> ('--foo', None)
	processArg('-f:bar')	-> ('-f', 'bar')
	processArg('-f=bar')	-> ('-f', 'bar')
	processArg('-f')	-> ('-f', None)
    """
    splitInd = find(arg, ':')
    if splitInd <> -1:
	return arg[0:splitInd], arg[splitInd+1:]

    splitInd = find(arg, '=')
    if splitInd <> -1:
	return arg[0:splitInd], arg[splitInd+1:]

    return arg, None

def removePrefix( arg ):
    if arg[0:2] == '--':
	return arg[2:]
    if arg[0] == '-':
	return arg[1:]

def indexOfOption( arg ):
    i = 1
    while i < len( gOpts ):
	opt = gOpts[i]
	if opt == arg:
	    return i
	if opt == removePrefix(arg):
	    if gValues[i]:
		return i
	if isShortOption(opt) and isShortOption(arg):
	    if find(opt, arg[1:]) <> -1:
		return i
	i = i+1
    return -1

def receivedOption_aux( arg ):
    "this deals only with single options"
    ret = indexOfOption( arg )
    if ret == -1:
	return 0
    else:
	try:
	    gNotused.remove(ret)
	except ValueError:
	# what to do here?
	# either arg was already processed as a value for another option
	# in which case there is an error on the command line, or somebody
	# simply called receivedOption() twice, which should not generate
	# an error. Have to return 1 because of sets of short options
	    pass
	return ret

def getValueOf_aux( opt ):
    """ a helper function for getValueOf, it only handles
	single options (getValueOf() handles lists of options) """
    ind = indexOfOption( opt )
    if ind == -1:
	raise ValueError(opt+' is not present in the command line')
    # is opt is a short option, i.e. '-hk', and is in a set with other
    # short options, i.e. '-hkj', then the user can't pass it a value
    # because there is not way to know whether the value is for '-hk'
    # or for '-j'
    if isShortOption( opt ):
	if gOpts[ind] <> opt and gOpts[ind] <> removePrefix(opt):
	    printProblemAndExit( [ind], 'an argument is expected for option "%s",\n\tseparate it from the set "%s" and pass it the argument' % (opt, gOpts[ind]))
    try:
	gNotused.remove(ind)
    except ValueError:
	# this is fine, this index might have been
	# removed by a call to receivedOption()
	pass

    if gValues[ind] == None:
	if (ind+1) >= len( gCmdLine ):
	    return None
	else:
	    if not isOption( gCmdLine[ind+1] ):
		try:
		    gNotused.remove(ind+1)
		except ValueError:
		    # this means that gCmdLine[ind+1] was already processed
		    # and is not in fact a value for the option opt, so:
		    return None
		return gCmdLine[ind+1]
	    else:
		return None
    else:
	return gValues[ind]

def getValueOf( opts, receivedDefault, desc, conversionFunc = None ):
    listOfOpts = split(opts)
    for opt in listOfOpts:
	try:
	    val = getValueOf_aux(opt)
	except ValueError:
	    # we have not received option o, continue
	    continue
	# return the default if we didn't get a value for the option
	# if there is one
	if (val == None) and not receivedDefault:
	    ind = indexOfOption( opt )
	    printProblemAndExit( [ind], 'option \"%s\" requires %s argument' % (gCmdLine[ind], desc) )
	    # and this code is no more of this earth
	if val <> None:
	    if conversionFunc:
		try:
		    val = conversionFunc( val )
		except ValueError:
		    ind = indexOfOption( opt )
		    printProblemAndExit( [ind], 'option "%s" requires %s argument, "%s" is not %s' % (opt, desc, val, desc) )
		    # and this code is no more of this earth
	return val

    # if we are here, then we have not received any
    # options matching the list opt passed in
    if receivedDefault:
	# let the caller handle this
	return None
    else:
	if ' ' in opts:
	    msg = re.sub(' +', '" or "', opts)
	else:
	    msg = opts
	printProblemAndExit( [], 'option "%s" with %s argument is expected' % (msg, desc) )

def printProblem( underline, msg, type='error' ):
    if type == 'error':
	print 'ERROR in command line:'
    elif type == 'warning':
	print 'WARNING about command line:'
    else:
	raise ValueError, 'error type must be either "error" or "warning"'

    # first duplicate the command line
    sys.stdout.write('\t')
    for arg in gCmdLine:
	sys.stdout.write(arg+' ')
    sys.stdout.write('\n')
    # now underline
    sys.stdout.write('\t')
    for ind in range( len( gCmdLine ) ):
	if ind in underline:
	    glyph = '^'
	else:
	    glyph = ' '
	if ind <> 0:
	    sys.stdout.write(' ')
	for i in range( len( gCmdLine[ind] ) ):
	    sys.stdout.write( glyph )
    sys.stdout.write('\n')
    # and print the message
    print  '\t' + msg + '\n'


def printProblemAndExit( underline, msg ):
    printProblem( underline, msg )
    raise SystemExit

def init(cmdline):
    global gOpts, gValues, gNotused, gCmdLine
    gCmdLine = cmdline
    gOpts = []; gValues = []; gNotused = []
    for arg in gCmdLine:
	opt, value = processArg( arg )
	gOpts.append( opt )
	gValues.append( value )
    gNotused = range(0, len( gCmdLine ) )
    gNotused.remove(0)

############# The Public Inteface to the Module ##############

def receivedOption( args ):
    argList = split(args)
    for arg in argList:
	if receivedOption_aux( arg ):
	    return 1
    return 0

# need to have *arg rather than default=None, so that the
# caller may pass None as the default, and we may may correctly
# distinguish that from the case when there is no default
def getStringValueOf( opts, *arg ):
    receivedDefault = 0
    if len(arg) > 0:
	receivedDefault = 1

    val = getValueOf( opts, receivedDefault, 'a string' )

    if val == None and receivedDefault:
	return arg[0]
	
    return val

def getIntegerValueOf( opts, *arg ):
    receivedDefault = 0
    if len(arg) > 0:
	receivedDefault = 1

    val = getValueOf( opts, receivedDefault, 'an integer', atoi )

    if val == None and receivedDefault:
	return arg[0]
	
    return val

def getFloatValueOf( opts, *arg ):
    receivedDefault = 0
    if len(arg) > 0:
	receivedDefault = 1

    val = getValueOf( opts, receivedDefault, 'a float', atof )

    if val == None and receivedDefault:
	return arg[0]
	
    return val

def getRemaining(desc=None):
    """
    Get the first argument, that is not an option, that has not
    been used so far
    """
    for i in range(len(gNotused)):
	arg = gCmdLine[gNotused[i]]
	if not isOption(arg):
	    gNotused.pop(i)
	    return arg

    # if there is a description of what was expected
    # then this argument was mandatory
    if desc:
	printProblemAndExit([], 'you must specify the %s' % desc)
    else:
	return None

def printWarnings():
    if gNotused <> []:
	printProblem( gNotused, 'The underlined arguments were ignored', 'warning' )
	return 1
    else:
	return 0

###################### The Tests ###########################

def testLongOptions():
    print 'Running the long options test...',
    cmdline = 'proggy --string=bla --int:5 --float 3.3 wrd=6'
    init( split(cmdline) )
    str = getStringValueOf('-s --string')
    int = getIntegerValueOf('-i --int')
    flt = getFloatValueOf('-f --float')
    wrd = getFloatValueOf('-w --wrd')

    if ([str, int, flt, wrd] <> ['bla', 5, 3.3, 6.0]):
	raise RuntimeError('long options test failed, values do not match')

    if printWarnings():
	raise RuntimeError('long options test failed, no warnings should have been printed')

    print 'done'

def testShortOptions():
    print 'Running the short options test...'
    cmdline = 'proggy -s=bla -i:5 unused options -f 3.3 w=6'
    init( split(cmdline) )
    str = getStringValueOf('-s --string')
    int = getIntegerValueOf('-i --int')
    flt = getFloatValueOf('-f --float')
    wrd = getFloatValueOf('-w --wrd')

    if ([str, int, flt, wrd] <> ['bla', 5, 3.3, 6.0]):
	raise RuntimeError('short options test failed, values do not match')

    if not printWarnings():
	raise RuntimeError('short options test failed, warnings should have been printed')

    print 'done'

def testTypeChecking():
    print 'Running the type checking test...'
    cmdline = 'proggy int=3.3 int2=foo float=foo'
    init( split(cmdline) )

    try: getIntegerValueOf('--int')
    except SystemExit: pass
    else: raise RuntimeError('integer type checking is broken')

    try: getIntegerValueOf('--int2')
    except SystemExit: pass
    else: raise RuntimeError('integer type checking is broken')

    try: getFloatValueOf('--float')
    except SystemExit: pass
    else: raise RuntimeError('float type checking is broken')

    print 'done'

def testShortOptionSets():
    print 'Running the set of short options test...'
    cmdline = 'proggy -abc -ghi=bar'
    init( split(cmdline) )

    if not receivedOption('-a'):
	raise SystemExit
    if not receivedOption('-bc'):
	raise SystemExit
    if receivedOption('--abc'):
	raise SystemExit

    try: getStringValueOf('-gh')
    except SystemExit: pass
    else: raise RuntimeError('incorrectly getting value for elements of sets of short options')

    print 'done'

def testDefaults():
    print 'Running the defaults test...'
    cmdline = 'proggy present=4 --should-have-arg --noarg1 --noarg2'
    init( split(cmdline) )

    prt = getIntegerValueOf('--present', 3)
    narg1 = getIntegerValueOf('--noarg1', 5)
    narg2 = getFloatValueOf('--noarg2', None)
    abst = getStringValueOf('-a --absent', 'foo')
    try: getIntegerValueOf('--should-have-arg')
    except SystemExit: pass
    else: raise RuntimeError('defaults test failed, should have complained about an absence of a mandatory argument')
    try: getFloatValueOf('-nd --absent-no-default')
    except SystemExit: pass
    else: raise RuntimeError('defaults test failed, should have complained about an absent option')
    # need this to get complete code coverage in getValueOf
    try: getStringValueOf('--absent-no-default-single-option')
    except SystemExit: pass
    else: raise RuntimeError('defaults test failed, should have complained about an absent option')

    if [prt, narg1, narg2, abst] <> [4, 5, None, 'foo']:
	raise RuntimeError('defaults test failed, values do not match')
    
    if printWarnings():
	raise RuntimeError('defaults test failed, no warnings should have been printed')

    print 'done'

def testParameters():
    print 'Running the parameter test...'
    cmdline = 'proggy -foo=3.3 param1 --bar bar param2'
    init( split(cmdline) )

    getFloatValueOf('-foo')
    getStringValueOf('--bar')

    prm1 = getRemaining()
    prm2 = getRemaining()
    prm3 = getRemaining()

    try: getRemaining('necessary parameter')
    except SystemExit: pass
    else: raise RuntimeError('parameter test failed, should have complained because of an absent necessary parameter')

    if [prm1, prm2, prm3] <> ['param1', 'param2', None]:
	raise RuntimeError('parameter test failed, values do not match')

    if printWarnings():
	raise RuntimeError('parameter test failed, no warnings should have been printed')

    print 'done'

def testCompleteCoverage():
    """a few contrived situations that did not really fit into any other
       categories, put here to get complete code coverage"""
    print 'Running the rest of the tests to get complete code coverage...',
    cmdline ='proggy --opt bla=3 --foo 5'
    init( split(cmdline) )
    
    if int(getRemaining()) <> 5:
	raise RuntimeError('code coverage test failed, values do not match')
    if not receivedOption('--opt'):
	raise RuntimeError('code coverage test failed, option is actually present')
    if getIntegerValueOf('--opt', 3) <> 3:
	raise RuntimeError('code coverage test failed, values are different')
    if getIntegerValueOf('--foo', None):
	raise RuntimeError('code coverage test failed, got back a phantom value')

    print 'done'

def test():
    "test the module, it breaks by raising a RuntimeError"
    try:               testLongOptions()
    except SystemExit: raise RuntimeError('long options test failed')

    try:               testShortOptions()
    except SystemExit: raise RuntimeError('short options test failed')

    try:               testTypeChecking()
    except SystemExit: raise RuntimeError('type checking test failed')

    try:               testShortOptionSets()
    except SystemExit: raise RuntimeError('short option sets test failed')

    try:               testDefaults()
    except SystemExit: raise RuntimeError('defaults test failed')

    try:               testParameters()
    except SystemExit: raise RuntimeError('parameters test failed')

    try:		testCompleteCoverage()
    except SystemExit: raise RuntimeError('code coverage test failed')

################### Initialization #########################

init(sys.argv)

############################################################

if __name__ == '__main__':
    if receivedOption('-v --version'):
	print __version__
	raise SystemExit

    if receivedOption('-t --test'):
	print '\n================================================================'
	print 'About to test this module.'
	print 'A lot of stuff will be printed to the screen.'
	print 'All the tests pass unless the program dies with some exception'
	print 'even though "ERROR" or "WARNING" may be printed several times.'
	print '================================================================\n'

	if not receivedOption('-n --no-warning'):
	    print 'Press any key to begin the tests'
	    sys.stdin.read(1)

	try:
	    test()
	except RuntimeError, msg:
	    print 'RuntimeError:', msg
	    print '\n================================================================'
	    print 'FAILED'
	    print '================================================================\n'
	else:
	    print '\n================================================================'
	    print 'all tests PASSED'
	    print '================================================================\n'

    else:
	print 'Run with "-t" or "--test" to test this module'

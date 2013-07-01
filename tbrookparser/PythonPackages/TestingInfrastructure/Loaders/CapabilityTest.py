import os, sys

if __name__=='__main__':
    print "\n for component test in %s \n" % __file__
    thisdir    = os.path.abspath(os.path.dirname(__file__))
    testingdir = thisdir + '/../'
    sys.path.append(testingdir)
    parserdir  = thisdir + '/../../TBrookParser'
    sys.path.append(parserdir)
    
from TestCases   import TruchasCapabilityTest

class CapabilityTest(TruchasCapabilityTest):
    "Base template for developer to build on"

    def testSomethingNow(self):
	"single line comment describing the test"


if __name__=='__main__':
    "put component test here"

    print " this component test is empty\n"

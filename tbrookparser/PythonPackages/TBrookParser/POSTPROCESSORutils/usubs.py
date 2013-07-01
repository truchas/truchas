"""

 usubs

 -----------------------------------------------------------------------------
  Purpose:
  
     usubs is a set of utility routines for getting user input.
     Typically the user is given a command and the results are
     interpreted on a word-by word basis. Strings enclosed in double
     or single quotes are interpreted as literals.

     Remember, more than a set of subroutines, usubs is an input
     philosophy.  To really buy into it, you must do _ALL_ of your
     input through usubs.  Mixing usubs calls with regular input calls
     can have disastrous results.

  Public Interface(s):
    i = input([inBuffer,noNewInput,debug])
    i.usetc([prompt,default])
    i.usetf([prompt,default])
    i.useti([prompt,default])
     
  Usage Examples:
     Some examples of how you would use this:
     
     
     >>> import usubs
     >>> input = usubs.input()
     >>>
     >>> i = input.uset('Please give a value for i' [, default]);
         ... j = input.uset('Please give a value for j', 0)
        Please give a value for i (None): 24 36
     >>> print i, j
      24 36
     >>>
     
     Note that the prompt was not displayed for j because the user 'typed
     ahead' and provided an answer because they knew something was
     forthcoming.
       
     This is a partial implementation which provides:
     
       uset(prompt, default)  General gset routine
       useti(prompt, default) Returns an integer
       usetf(prompt, default) Returns a float/real
       usetc(prompt, default) Returns a string (synonym for uset())
       
       @filename              execute all commands in file.mac
     
     All commands take a default value which is returned if the user
     types a return or if the user types a '-' the default value is
     returned.

  Contains:
    class input
        __init__(self, inBuffer=None, noNewInput=0, debug=0)
        nextWord(self, debug=0)
        ugetLine(self)
        condition(self,string)
        uset(self, prompt, default=None)
        usetc(self, prompt=None, default=None)
        usetf(self, prompt=None, default=None)
        useti(self, prompt=None, default=None)

    ---> NO Unit Test Block
  
  History:
     Usubs is a clone of a set of routines that were originally developed
     by Art Voter in T12 in f77 which I had converted to C.  Turns out that
     the python implementation is straightforward because python will do
     the loops for you.  So, I haven't bothered to program in the loops
     
  Author(s):  Sriram Swaminarayan (sriram@lanl.gov)
 -----------------------------------------------------------------------------
"""


import sys
from copy import copy
import re

class input:
    buffer   = []
    dash     = re.compile('^-+$')
    maxTries = 4                

    def __init__(self, inBuffer=None, noNewInput=0, debug=0):
        self.debug = debug
        if(inBuffer != None):
            self.buffer = self.buffer + self.condition(inBuffer)
        if(self.debug):
            print self.buffer
        self.noNewInput = noNewInput
        
    def nextWord(self, debug=0):
        if(not len(self.buffer)):
            if(not self.noNewInput):
                self.ugetLine()
            else:
                if(self.debug): print 'returning nothing because no new input!'
                return('')
        
        if(self.debug): print 'buffer length is: ', len(self.buffer)
        if(not len(self.buffer)):
            return('')

        word = self.buffer.pop(0)
        if(self.debug): print 'word is: ', word
        while(word[0] == '@'):
            try:
                if(debug or self.debug): print 'getting data from file', word[1:]
                file=open(word[1:],'r')
                lines = file.readlines()
                lines.reverse()
                if(debug or self.debug): print 'lines are: ', lines
                for line in lines:
                    self.buffer = self.condition(line) + self.buffer
                if(debug or self.debug): print 'buffer is : ', self.buffer
                
                file.close()
            except:
                # do nothing
                print '\n   Unable to get commands from file ',word[1:],'\n'
                if(debug or self.debug):
                    print 'rejected word'
            if(len(self.buffer)):
                word=self.buffer.pop(0)
            else:
                return('')

        # Replace any dashes with one less
        m = self.dash.search(word)
        if(m):
            word = word[1:]

        return(word)

    def ugetLine(self):
        print self.prompt,
        sys.stdout.flush()
        a = sys.stdin.readline()
        self.buffer = self.buffer + self.condition(a)

    def condition(self,string):
        """
         Here is where we do our post-processing of input.
         Basically we want to crack symbols, identify
         string input, etc.  Of course, other than double
         quoted strings, none of this has been implemented
         yet, but I am working on it!
         """
        
        retval = []
        s = string
        if(s == None):
            return(retval)

        # first double quoted strings
        while(len(s)):
            i1 = s.find('"')
            i2 = s.find('"',i1+1)
            if(i1 < 0 or i2 < 0):
                retval = retval + s.split()
                return(retval)
            if(i1>0):retval = retval + s[:i1].split()
            retval = retval + [s[i1+1:i2]]
            s = s[i2+1:]


        return(retval)    

    def uset(self, prompt, default=None):
        """
        The general uset routine.  This routine does not care what
        gets returned.
        """
        #print "uset prompt:", prompt
        
        if(type(prompt) != type('')):
            print """
            Sorry,  prompt must be a string and is required
            """
            raise TypeError

        abak = copy(default) # Backup our default value

        a = abak
        self.prompt = ' ' + prompt
        if(a != None):
            self.prompt = self.prompt + '(' + str(a) + ') '
        else:
            self.prompt = self.prompt + '() '
        tmp = self.nextWord()
        if(tmp == None or tmp == ''):
            tmp = a
        return(tmp)
        
    def usetc(self, prompt=None, default=None):
        """
        Return a character string array.
        Returns default if user types a '-'
        or a carraige return

        Note that any string that contains
        spaces must be enclosed in double quotes

        In some senses, this is just a synonym for uset
        """
        #print "usetc:", prompt
        prompt = str(prompt)

        a = self.uset(prompt, default)

        return(a)

    def usetf(self, prompt=None, default=None):
        """
        Return a float.
        Does a check to ensure that the user
        enters a floating point number
        """

        i = 0
        abak = copy(default) # Backup our default value

        while(i<self.maxTries):
            tmp = self.uset(prompt,default)
            try:
                a = float(tmp)
                i = self.maxTries # preload failure
            except:
                # Print warning
                print
                print "  WARNING: Invalid Entry.  Please enter a floating point number "
                print 
                # reload the default
                a = abak
                i = i+1

        return(a)

    def useti(self, prompt=None, default=None):
        """
        Return an integer
        
        Does a check to ensure that the user
        enters an integer number
        
        Floating point numbers are rounded down
        """
        
        i = 0
        abak = copy(default) # Backup our default value

        a = abak
        while(i<self.maxTries):
            tmp = self.uset(prompt,default)
            try:
                a = float(tmp)
                a = int(a)
                i = self.maxTries # preload failure
            except:
                # Print warning
                print
                print "  WARNING: Invalid Entry.  Please enter an integer!!"
                print 
                # reload the default
                a = abak
                i = i+1
                
        return(a)
    

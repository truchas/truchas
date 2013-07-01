"""

 getXMLString

 -----------------------------------------------------------------------------
  Purpose:

    Make an XML file into a string
    based on the criteria for TBrookParsing

    The approach is to apply regular expressions to a single
    string, eliminating that which is not needed by TBrook

    Checking of the formed string for tag-coherency is provided as an
    option (on by default).  It can be turned off to allow for faster
    processing if the user knows the XML is well formed.

  Public Interface(s):

    x = XMLstr(file) # instantiates object and checks it (by default)

    x = XMLstr(file, checkIt=False) # instantiates object but doesn't
                                    # check its tags for coherency
    x.fstr           # the resulting XML string itself
    x.checked        # True if the result has been checked for tag coherency
    x.str()          # prints the XML string

   The following methods are also available.  They are invoked on
   initialization if checkIt==True (the default).  They can be invoked
   manually, if desired, to (re)verify a string or to retrieve the
   tags themselves.

    (o,e) = x.getTags()  # retrieves lists of opening and closing tags
                         # NOTE: the two tag lists will not be identical
                         #       if there are ANY nested tags.
    n     = x.checkTags()# if n > 0, there are open tags (to close)
                         # if n==-1, the tags are not coherent or fixable
    if n: x.fixTags()    # closes unclosed tags

  Contains:
    class XMLstr
        __init__(self,file, checkIt=True, debug=0)
        __formString(self, debug=0)
        __processCycles(self, str)
        getTags(self)
        checkTags(self)
        fixTags(self)
        str(self)

    Unit Test Block

  Version:
    $ID$

  Author(s): Larry Cox (ljcox@lanl.gov)
 -----------------------------------------------------------------------------
"""

if __name__=='__main__':
    "For component testing, set sys.path"
    import os, sys
    thisdir   = os.path.dirname(__file__)
    parsedir  = os.path.abspath(thisdir+'../')
    sys.path.append(parsedir)

class XMLstr:
    def __init__(self,file,checkIt=True,debug=0):
        "Instantiate the object"
        self.file    = file    # file name
        self.fstr    = None    # initialize string to Nada
        self.checked = False   # has tag coherency been checked?
        self.debug   = debug   # internal debug flag
        self.__etags = []      # holder for extra tags, if found on checkTags
        self.__formString()    # form the real string
        if checkIt and self.checkTags():
            self.fixTags()     # closes unclosed tags in __etags

    def __formString(self):
        """
        Suck in the whole file into one big string and purge it
        of unwanted elements, then store the string in the object
        """
        import sys, re
        fp   = open(self.file,'r')
        fstr = ''.join(fp.readlines())
        fp.close()

        "Replace spurious XML characters"
        fstr = fstr.replace(' < ',' lt ')
        fstr = fstr.replace(' > ',' gt ')
        fstr = fstr.replace('&','and')

        "Remove all XML comments"
        myRe = re.compile('<\!\[CDATA.*?\]>',re.DOTALL|re.MULTILINE)
        fstr = myRe.sub('',fstr)
        fstr = re.sub('<\?xml.*?\?>','',fstr)

        "Process CYCLE elements"
        fstr = self.__processCycles(fstr)
                    
        "Remove legacy text"
        myRe = re.compile('^\[[^<]*',re.MULTILINE)
        fstr = myRe.sub('',fstr)

        "Special case for legacy text interrupted by preserved <CYCLE>s"
        myRe = re.compile('</CYCLE>[^<]*',re.MULTILINE)
        fstr = myRe.sub('</CYCLE>\n',fstr)

        "Compress whitespace"
        fstr = re.sub('^[^<]*','',fstr) # Toss chars before the first tag
        myRe = re.compile('\n(\s*\n)+',re.MULTILINE)
        fstr = myRe.sub('\n',fstr)      # Toss embedded blank lines
        fstr = re.sub('  +',' ',fstr)   # Trim multiple spaces to one

        "Handle interrupted tags: trim from end back to last '>'"
        lc = fstr.rfind('>')
        if lc < len(fstr)-1:
            if self.debug: print "dropping '", fstr[lc+1:],"'"
            fstr = fstr[:lc+1] + '\n'
        
        self.fstr = fstr

    def __processCycles(self,fstr):
        """
        Keeps only CYCLEs that contain  LINEAR_RESIDUAL, 
          NONLIN_RESIDUAL or JOULEHEAT information (tags)
        Assumptions: LINEAR/NONLIN will be in last cycle only
                     JOULEHEAT may be in any cycle, but not in all
        """
        import re
        cycleRe = re.compile('<CYCLE.*?</CYCLE>',re.DOTALL|re.MULTILINE)
        cycles = cycleRe.findall(fstr) # A list of CYCLE elements as strings 
        if cycles:                     # There may not be any
            cstr = ''.join(cycles)     # one big CYCLE string
            cpos = cycleRe.search(fstr).start()
            if self.debug:
                print "Processing CYCLES..."
                print "\tThere is(are) %d CYCLE(s)." %(len(cycles))
                cc = 0
                for c in cycles: cc += len(c)
                print "\tTotal character length: %d" %cc

            jHeat = re.search('<JOULEHEAT',cstr)
            if jHeat:
                if self.debug:print "\tDetected JOULEHEAT..."
                "Iterate over cycles"
                jHeat  = ''
                jc = 0
                for c in cycles:
                    if 'JOULEHEAT' in c:
                        jc +=1
                        jHeat += c + '\n'
                if self.debug:
                    print "\t\tThere are %d cycles." %(len(cycles))
                    print "\t\tKeeping %d with JOULEHAEAT: total length %d" \
                          %(jc,len(jHeat))

            resid = re.search('<LINEAR_RESIDUAL|<NONLIN_RESIDUAL',cstr)
            if resid:
                if self.debug: print "\tDetecting RESIDUALs..."
                "Preserve last cycle"
                resid = cycles[-1] + '\n'
                if 'JOULEHEAT' in resid:
                    resid = None
                elif self.debug:
                    print "\t\tKeeping the last cycle of length %d" \
                          %(len(resid))
            else: resid = None

            if self.debug:
                print "\tDeleting all cycles..."
                print "\t\tlength before: %d" %(len(fstr))
            fstr = cycleRe.sub('',fstr)
            if self.debug: print "\t\tlength after:  %d" %(len(fstr))

            if jHeat:
                if self.debug: print "\tReinserting the JOULEHEAT cycle(s)"
                fstr = fstr[:cpos] + jHeat + fstr[cpos:]
                cpos = cpos + len(jHeat)

            if resid:
                if self.debug: print "\tReinserting the RESIDUALs cycle"
                fstr = fstr[:cpos] + resid + fstr[cpos:]

        return fstr

    def getTags(self):
        " Return lists of opening and closing XML tags"
        if not self.fstr:
            print "No String to check. formString() first."
            return (None,None)
        else:
            import re
            oIter = re.finditer('<(\w+).*?>',self.fstr) # opening tags
            eIter = re.finditer('</(\w+)>',self.fstr)   # closing tags
            "Trim each item to just the TAG keyword"
            olist = []
            for o in oIter:
                olist.append(o.group(1))
            elist = []
            for e in eIter:
                elist.append(e.group(1))
            return (olist,elist)

    def checkTags(self):
        if self.debug: print "Checking tags coherency..."
        "Look for unclosed tags"
        if not self.fstr:
            print "No String to check. formString() first."
            return None
        elif self.checked:
            return len(self.__etags)
        else:
            (olist,elist) = self.getTags()
            olist.reverse() # reverse to make a stack
            nest = ['TAIL'] # stack for out-of-order (nested) tags
            for i in range(len(elist)):
                "Each iteration must consume one each of open/close tags"
                e = elist[i] # next closing tag
                while e:
                    "descend olist and nest until we find e"
                    if e == nest[-1]:
                        e = None   # consumes e
                        nest.pop() # consumes last nested opener
                    else:
                        try:
                            o = olist.pop()
                            if e == o:
                                e = None       # consumes e and, implicitly, o
                            else:
                                nest.append(o) # consumes neither e nor o
                        except:
                            print "  !!! Unrepairable XML tag incoherency encountered!!!"
                            print "      Likely caused by a spurious end tag",\
                                        "or a missing opening tag"
                            print "\tOpening tag list emptied prematurely\n",\
                                  "\t\tNested depth %d\n\t\tElist length %d" \
                                  %(len(nest)-1,len(elist))
                            print "\tCurrent End Tag: </%s>\n\tNested tag list:" %e
                            for n in nest[1:]:
                                print "\t\t<%s>" %n
                            return -1
            "Store extra tags"
            if len(nest) > 1: self.__etags = nest[1:] # nested ones
            while len(olist) > 0:
                self.__etags.append(olist.pop())      # extra ones
            self.checked = True
            return len(self.__etags)

    def fixTags(self):
        "Close remaining tags"
        if not self.fstr:
            print "No String to check. formString() first."
        elif len(self.__etags) > 0:
            if self.debug:
                print "  Found %d leftover opening tags." %len(self.__etags)
            while len(self.__etags)>0:
                t = self.__etags.pop()
                if self.debug: print '\t<%s>' %(t)
                self.fstr += "</%s>\n" %(t)

    def str(self):
        print self.fstr

if __name__ == '__main__':

    from PYTHONutils import uTestOpts
    p = uTestOpts('fdpct',
                  defaults={'p':False, # print XML string generated
                            'c':False, # autocheck tag coherency
                            't':False},# print tag list
                  actions ={'p':'store_true',
                            'c':'store_true',
                            't':'store_true'})
    (opt,a) = p.parse_args()
    p.header(__file__,opt.f)

    try:
        x = XMLstr(opt.f,debug=opt.d,checkIt=opt.c)
        (ol,el) = x.getTags()
        print "Number of tags:\n\topening = %d\n\tclosing = %d" \
              %(len(ol),len(el))
        if opt.t:
            print "Tag List:"
            for i in range(max(len(ol),len(el))):
                (o,e) = (' ',' ')
                if  i < len(ol): o = ol[i]
                if  i < len(el): e = el[i]
                print "%3d:\t%25s\t%25s" %(i,o,e)
        if opt.p:
            x.str()
        print "Tags Checked:",x.checked
        print "XML string length = %d" %(len(x.fstr))

        if not x.checked:
            n = x.checkTags()
            if n > 0:
                print "\tTags not complete: %d extra tags" %(n)
                x.fixTags()
            elif n != -1:
                print "\tTags checkout okay."

    except:
        print "---> Test failed in some aspect <---"
        print "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise

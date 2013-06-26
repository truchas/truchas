"""

 Parsing

 -----------------------------------------------------------------------------
  Purpose:
  
    To define a very general document parser, this could be UDM, HDF or
    XML the restriction is that we assume the document contains
    elements, with each element containing element attributes and
    element data

    This module also creates a specific Parser subclass (MiniDom).
  
  Public Interface(s):

    For all parsers
      p = Parser
      p.getDoc(file)
      p.getElemAttribute(element, AttName)
      p.getElemParent(element)
      p.getElement(TheNode, TagName)
      p.getElements(TheNode, TagName)
      p.getValue(element)

    Additionally for MiniDom parsers
      m = MiniDom
      m.getXMLFileAsString(filename)
      m.typeConverter(typeCode)
  
  Contains:
    class Parser
        getDoc(self,file,debug=0)
        getElements(self,TheNode,TagName)
        getElement(self,TheNode,TagName)
        getValue(self,element)
        getElemAttribute(self,element,AttName)
        getElemParent(self,element)

    class MiniDom(Parser)
        getDoc(self,file)
        getXMLFileAsString(self,filename)
        getElement(self,TheNode,TagName)
        getElements(self,TheNode,TagName)
        getValue(self,element)
        getElemParent(self,element)
        getElemAttribute(self,element,AttName)
        typeConverter(self, typeCode)

    Unit Test Block
  
  Author(s): Sharen Cummins (scummins@lanl.gov)
 -----------------------------------------------------------------------------
"""
import sys
if __name__=='__main__':
    # Set sys.path for component testing mode
    import os
    thisdir   = os.path.dirname(__file__)
    parsedir  = os.path.abspath(thisdir+'../')
    sys.path.append(parsedir)

class Parser:

    "A very general document parser"
    
    def getDoc(self,file):
        print
        print 'In getDoc'
        print 'Obtains the highest node in the file'
        print 'All other nodes are children of this node'
        print
        
    def getElements(self,TheNode,TagName):
        print
        print 'In getElements'
        print 'Obtains a list of nodes with node name', \
              '"TagName" that live directly underneath '
        print '(i.e are children of) the node "TheNode"'
        print

    def getElement(self,TheNode,TagName):
        print
        print 'In getElement'
        print 'Obtains the node with node name "TagName" that lives',\
              'directly underneath (i.e is a child of) the node "TheNode"'
        print
        
    def getValue(self,element):
        print
        print 'In getValue'
        print 'Obtains the data associated with the node "element"'
        print 

    def getElemAttribute(self,element,AttName):
        print
        print 'In getElemAttribute'
        print 'Obtains the attribute with the name "AttName"',\
              'associated with the node "element" '
        print

    def getElemParent(self,element):
        print
        print 'In getElemParent'
        print 'Obtains the nearest parent node that contains',\
              'the node "element"'
        print 

class MiniDom(Parser):

    "A particular Parser subclass (MiniDom)" 

    def getDoc(self,file,fpwatch=sys.stdout,debug=0):

        import xml
        import xml.dom.minidom

        """
        extract only relevant (cycle) data from XML file
        
          ensure all tags in the xml file match,
          if not then repair the xml file (???)
        """
        
        try:
            import sys
            xmlstr = self.getXMLFileAsString(file)
            T      = xml.dom.minidom.parseString(xmlstr)
        except:
            print >> fpwatch, '\n In MiniDom parser, reading in entire file \n'
            T      = xml.dom.minidom.parse(file)

        return T   

    def getXMLFileAsString(self,filename):

        from XMLutils import XMLstr
        s = XMLstr(filename)
        return s.fstr
   
    def getElement(self,TheNode,TagName):

        elem = TheNode.getElementsByTagName(TagName)[0]
        elem.normalize()
        return elem

    def getElements(self,TheNode,TagName):

        elems = TheNode.getElementsByTagName(TagName)
        for elem in elems:
            elem.normalize()
        return elems

    def getValue(self,element):

        if element.firstChild != None:
            S = element.firstChild.nodeValue
            if (    type(S.encode()) == str
                and len(S) > 1):
                S = S[1:len(S)-1]
        else:
            S = ''
        return S

    def getElemParent(self,element):

        return element.parentNode

    def getElemAttribute(self,element,AttName):

        return element.getAttribute(AttName)

    def typeConverter(self, typeCode):
        
         "Given a Numeric type code, return the matching Python type converter."
         conversionDict = {'d': float, 'f':float, 'i': int, 'c':str}

         return conversionDict[typeCode]


if __name__=='__main__':

    'test "Parser"-ing an XML file'

    from PYTHONutils import uTestOpts
    opts = uTestOpts('fds',defaults={'s':False},actions={'s': 'store_true'})
    (opt,args) = opts.parse_args()
    opts.header(__file__,opt.f)

    fpwatch = sys.stdout
    try:
        parse = MiniDom

        if opt.s:
            "Just check getXMLFileAsString"
            xmlstr = parse().getXMLFileAsString(opt.f)
            print >> fpwatch, 'length of XML string: %d' %(len(xmlstr))
            print >> fpwatch, 'xml-'*18
            print >> fpwatch, xmlstr
            print >> fpwatch, 'xml-'*18

        else:
            theDoc    = parse().getDoc(file=opt.f, debug=opt.d, fpwatch=fpwatch)
            variable  = parse().getElement(theDoc,'PROGRAMSPECIFICATIONS')
            
            thisspec  = 'CODE'
            specnode  = parse().getElement(variable,'CODE')
            specvalue = parse().getValue(specnode)
            
            print >> fpwatch,'the %s simulation specification is:' %(thisspec)
            print >> fpwatch, specvalue
            
            print >> fpwatch, 'meshes in this simulation are:'
            meshelems = parse().getElements(theDoc,'MESH')
            for elem in meshelems:
                filenode = parse().getElement(elem,'FILE')
                filefmt  = parse().getElemAttribute(filenode,'FORMAT')
                file     = str(parse().getValue(filenode))
                s        = file + ' which is in a ' + filefmt + ' format'
                print >> fpwatch, s

            print >> fpwatch, 'timesteps in this simulation are:'
            timelems = parse().getElements(theDoc,'TIMESTEP')
            for elem in timelems:
                filenode = parse().getElement(elem,'FILE')
                filefmt  = parse().getElemAttribute(filenode,'FORMAT')
                file     = str(parse().getValue(filenode))
                s        = file + ' which is in a ' + filefmt + ' format'
                print >> fpwatch, s
    except:
        print >> fpwatch, "---> Test failed in some aspect <---"
        print >> fpwatch, "\nNature-Of-Error:", sys.exc_info()[0],"\n"
        if opt.d: raise

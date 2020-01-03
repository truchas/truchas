#===============================================================================
#
# This file is part of Truchas. 3-Clause BSD license; see the LICENSE file.
#
#===============================================================================

from docutils import nodes
from docutils.parsers.rst import Directive, directives
from docutils.statemachine import StringList

from sphinx.locale import _
from sphinx.util.docutils import SphinxDirective

# Class name for all 'namelist parameter' elements in the document
NAMELIST_PARAMETER_CLASS = 'nml-param'

class namelist_parameter(nodes.block_quote):
    pass

def visit_nmle_node(self, node):
    self.visit_block_quote(node)

def depart_nmle_node(self, node):
    self.depart_block_quote(node)


class NamelistParameterDirective(SphinxDirective):
    """
    A custom directive that describes a namelist parameter.
    
    The directive takes three mandatory options:
        
        :type:    - the Fortran type of the namelist parameter
        :domain:  - the set of valid values for the parameter
        :default: - the default value if the parameter is not specified
    """
    
    has_content = False
    option_spec = {
        'type': directives.unchanged_required,
        'domain': directives.unchanged_required,
        'default': directives.unchanged_required,
    }
    
    def run(self):
        # Check all options were passed
        for option in self.option_spec.keys():
          if not option in self.options:
            raise ValueError("option '{}' required but none supplied".format(option))
        
        # Construct text tokens to be parsed
        text = []
        text.append('**Type:** ``{}``'.format(self.options['type'].upper()))
        text.append('')
        text.append('**Domain:** {}'.format(self.options['domain']))
        text.append('')
        text.append('**Default:** {}'.format(self.options['default']))
        
        content = StringList(text)
        
        # Create namelist parameter and add the parsed content
        nml_node = namelist_parameter(content, classes=[NAMELIST_PARAMETER_CLASS])
        self.state.nested_parse(content, self.content_offset, nml_node)

        return [nml_node]


def setup(app):
    app.add_node(namelist_parameter,
                 html=(visit_nmle_node, depart_nmle_node),
                 latex=(visit_nmle_node, depart_nmle_node),
                 text=(visit_nmle_node, depart_nmle_node))

    app.add_directive('namelist_parameter', NamelistParameterDirective)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }

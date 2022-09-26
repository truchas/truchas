# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import sphinx_rtd_theme
import sphinx.ext.imgmath

# -- Project information -----------------------------------------------------

project = 'Truchas Reference Manual'
copyright = '2021, LANL'
author = 'LANL'

# The full version, including alpha/beta/rc tags
release = '1.0.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_rtd_theme"
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'


# html_show_sourcelink = False



html_context = {
  'display_github': False,
  'github_repo': 'truchas',
  
  #'gitlab_user': 'naren.ragav',
  #'gitlab_version': 'master/docs/',
}


numfig = True 
numfig_format = {'figure': 'Figure %s'}
math_numfig = True
numfig_secnum_depth = 2
math_eqref_format = "Eq.{number}"
imgmath_latex_preamble = r'\usepackage{array}'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

def setup(app):
   app.add_css_file('css/custom.css')

extensions = ['sphinxcontrib.bibtex']
bibtex_bibfiles = ['refs.bib']

latex_engine = 'pdflatex'
latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '10pt',
    }

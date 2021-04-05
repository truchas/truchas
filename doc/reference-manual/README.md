# Truchas Reference Manual Source

This document explains how to build and contribute to the Truchas Reference
Manual written using Sphinx.

## Prerequisites

Building the Reference Manual requires Python and the following Python
packages:

* sphinx (https://www.sphinx-doc.org/en/master/usage/installation.html)
* sphinx-rtd-theme (https://sphinx-rtd-theme.readthedocs.io/)
* sphinxcontrib-bibtex (https://readthedocs.org/projects/sphinxcontrib-bibtex/)

These python packages can all be installed with pip:

    $ pip install sphinx sphinx-rtd-theme sphinxcontrib-bibtex

## How To Build

### Manually
To compile the Reference Manual directly, use the `sphinx-build` command:

    $ mkdir build
    $ sphinx-build doc/reference-manual build/html

The generated HTML files are located in `build/html`

### Using CMake
To build the Reference Manual as part of the CMake build, add the argument
`-D BUILD_HTML=ON` to the initial `cmake` command line:

    $ mkdir build
    $ cd build
    $ cmake -C ../config/<config_file> \
            -D TRUCHAS_TPL_DIR=<truchas_tpl_dir> \
            -D BUILD_HTML=ON ..
    $ make

The generated HTML files are located in `build/doc/reference-manual/html`.
See [BUILDING.md](BUILDING.md) for more information on building Truchas.

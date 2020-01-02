# Truchas Documentation

## Requirements
Compiling the Truchas documentation requires the following Python packages.

* sphinx
* sphinx-jsonschema
* sphinxcontrib-bibtex
* sphinx_bootstrap_theme

All requirements can be installed with pip.

    $ pip install sphinx sphinx-jsonschema sphinxcontrib-bibtex sphinx_bootstrap_theme

For convenience, we recommend installing Sphinx's autobuild extension to
automatically build the documentation whenever a file is changed.

    $ pip install sphinx-autobuild


## Building

### cmake
To compile the Truchas documentation as part of the cmake build, add the `-D BUILD_HTML=ON`
definition flag to the initial cmake command:

    $ mkdir build
    $ cd build
    $ cmake -C ../config/<config_file> \
            -D TRUCHAS_TPL_DIR=<truchas_tpl_dir> \
            -D BUILD_HTML=ON ..
    $ make
    $ make install

See [BUILDING.md](BUILDING.md) for more information on building Truchas.

### Directly
To compile the Truchas documentation directly, use the `sphinx-build` command

    $ mkdir build_doc
    $ sphinx-build doc/sphinx build_doc/html

You can also use `sphinx-autobuild` to automatically recompile the documentation whenever a file is
changed.

    $ mkdir build_doc
    $ sphinx-autobuild doc/sphinx build_doc/html

This starts a lightweight local server that can be accessed from any web browser.

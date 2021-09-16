# Truchas Documentation

## Requirements
Compiling the Truchas documentation requires the following Python packages.

* sphinx
* sphinx-jsonschema
* sphinxcontrib-bibtex
* sphinx-rtd-theme

All requirements can be installed with pip.

    $ pip install sphinx sphinx-jsonschema sphinxcontrib-bibtex sphinx-rtd-theme

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

### Manually
To compile the Truchas documentation directly, use the `sphinx-build` command

    $ mkdir build_doc
    $ sphinx-build doc/sphinx build_doc/html

You can also use `sphinx-autobuild` to automatically recompile the documentation whenever a file is
changed.

    $ mkdir build_doc
    $ sphinx-autobuild doc/sphinx build_doc/html

This starts a lightweight local server that can be accessed from any web browser.


## Style Guide
Refer to the Sphinx Documentation Style Guide (`doc/sphinx/style_guide.rst`) for tips and
best-practices for writing the Truchas Sphinx documentation.

The style guide is hidden from users (there are no visible links to it). However, it is still part
of the document tree and will be generated with the rest of the documentation. The HTML file is
located at `<build root>/style_guide.html`. When building with `sphinx-autobuild`, the style guide
page is located at [localhost:8000/style_guide.html](localhost:8000/style_guide.html).

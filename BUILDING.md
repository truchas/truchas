Building Truchas
------------------------------------------------------------------------------
### Quick Start Guide

#### Requirements
The Truchas build system assumes a UNIX-like environment. Current development
and testing is done on 64-bit Linux and Cray CLE platforms.
* Fortran and C/C++ compilers.  The compiler executables must be in your path.
  We use and test with the following compilers.
    - Intel Fortran and C/C++:
        - version 17.0.7 (anything 17.0.1 or later should be okay)
        - version 18.0.5 (anything 18.0.2 or later should be okay)
        - version 19.1.0. Version 19.0.x **will not** work due to compiler bugs.
    - NAG Fortran (with GNU C/C++):
        - version 6.1 (build 6149 or later preferred)
        - version 6.2 (build 6252 or later preferred)
        - most any version of GNU C/C++ should be okay
    - GFortran is *not* currently supported due to incomplete and/or flawed
      support for some Fortran 2003 features. There are GFortran configuration
      files for internal testing purposes, which you can try if you are feeling
      adventurous.
* Cmake version 3.5 or later; but not 3.6.0 and 3.6.1.
* Standard software development tools: make, patch, perl
* Zlib development library and header files
* Python, version 3.5 or later, along with the packages h5py (version 2.6.0 or
  later) and scipy (version 0.18.0 or later)
* MPI.  The C compiler wrapper (`mpicc`, for example) must be in your path.

Truchas requires some additional libraries, but these can be built by the
third party library build step described below.

#### Compiling
Compiling Truchas for the first time is usually a two-stage process.  The
first stage involves building and installing additional third party libraries
(TPL) that Truchas requires and which are not present on your system.  This
only needs to be done once.  A cmake superbuild project for this stage can be
found in the [truchas-tpl](https://gitlab.com/truchas/truchas-tpl) repository
on GitLab. This version of Truchas is tested against the "v10" bundle of TPLs;
do a `git checkout v10` after cloning the TPL repository. *v10 is required if
using the Intel 19.1 compiler*. See its README file for further instructions.

Once the required TPLs are installed, the procedure for building Truchas is
straightforward. You create a build directory, run cmake from that directory,
and then run make. What you choose for a build directory is irrelevant (other
than it must be different than the current directory). Here is an example:

    $ mkdir build
    $ cd build
    $ cmake -C ../config/intel-opt.cmake \
            -D TRUCHAS_TPL_DIR=<truchas_tpl_dir> ..
    $ make
    $ make install

* *Don't overlook the final `..` argument on the cmake command line!*
* The `-C` argument pre-loads the cmake cache with settings from the following
  file. The `config` subdirectory contains some examples. If none of those are
  suitable, create your own, or simply define the various variables directly
  on the cmake command line (using the `-D` flag).
* Set the `TRUCHAS_TPL_DIR` variable to the TPL installation directory you
  used in the first stage. It must be an absolute path.
* By default Truchas will be installed into the `install` subdirectory of the
  top-level source directory. Use the `-D CMAKE_INSTALL_PREFIX=<truchas_dir>`
  cmake argument to specify a different directory.

#### Testing
From the build directory run the command

    $ ctest

to run the regression test suite. On multi-core systems use the `-j<n>` option
to tell ctest how many processes it can run simultaneously; `-j8`, for example.
All tests should pass.

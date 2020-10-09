Building Truchas
===============================================================================
## Quick Start Guide

### Requirements
The Truchas build system assumes a UNIX-like environment. Current development
and testing is done on 64-bit Linux and Cray CLE platforms with both Intel and
NAG compilers and on MacOS with NAG Fortran (experimental).
* Fortran and C/C++ compilers.  The compiler executables must be in your path.
  We use and test with the following compilers.
    - Intel Fortran and C/C++:
        - version 18.0.5 (anything 18.0.2 or later should be okay)
        - version 19.1.0
        - Many versions are known to **not work** due to compiler bugs. These
          include 19.0.x (any x), 19.1.1, and 19.1.2.
    - NAG Fortran (with GNU C/C++):
        - version 6.2 (build 6252 or later preferred)
        - version 7.0 (build 7026 or later required)
        - most any version of GNU C/C++ should be okay
    - NAG Fortran on MacOS (with Apple Clang C/C++):
	    - NAG version 6.2 and 7.0 (build 7028 or later required)
		- Likely any version of Apple Clang should be fine
    - GFortran is *not* currently supported due to incomplete and/or flawed
      support for some Fortran 2003 features. There are GFortran configuration
      files for internal testing purposes, which you can try if you are feeling
      adventurous.
* Cmake version 3.16 or later
* Standard software development tools: make, patch, perl
* Zlib development library and header files
* Python, version 3.5 or later, along with the packages h5py (version 2.6.0 or
  later) and numpy (version 1.12.0 or later)
* MPI.  The C compiler wrapper (`mpicc`, for example) must be in your path.

Truchas requires some additional libraries, but these can be built by the
third party library build step described below.

### Compiling
Compiling Truchas for the first time is usually a two-stage process.  The
first stage involves building and installing additional third party libraries
(TPL) that Truchas requires and which are not present on your system.  This
only needs to be done once.  A cmake superbuild project for this stage can be
found in the [truchas-tpl](https://gitlab.com/truchas/truchas-tpl) repository
on GitLab. This version of Truchas is tested against the "v14" bundle of TPLs;
do a `git checkout v14` after cloning the TPL repository. See its README file
for further instructions.

Once the required TPLs are installed, the procedure for building Truchas is
straightforward. You create a build directory, run cmake from that directory,
and then run make. What you choose for a build directory is irrelevant (other
than it must be different than the current directory). Here is an example:

    $ mkdir build
    $ cd build
    $ cmake -C ../config/linux-intel.cmake \
            -D CMAKE_BUILD_TYPE=Release \
            -D TRUCHAS_TPL_DIR=<truchas_tpl_dir> ..
    $ make
    $ make install

* *Don't overlook the final `..` argument on the cmake command line!*
* The `-C` argument pre-loads the cmake cache with settings from the following
  file. The `config` subdirectory contains some examples. If none of those are
  suitable, create your own, or simply define the various variables directly
  on the cmake command line (using the `-D` flag).
* `Release`directs CMake to configure an optimized build of Truchas. Another
  option is `Debug` for an unoptimized build with lots of additional runtime
  checking.
* Set the `TRUCHAS_TPL_DIR` variable to the TPL installation directory you
  used in the first stage. It must be an absolute path.
* By default Truchas will be installed into the `install` subdirectory of the
  top-level source directory. Use the `-D CMAKE_INSTALL_PREFIX=<truchas_dir>`
  cmake argument to specify a different directory.

### Compiling on a Mac
The test suite is currently failing so mac support is still considered
experimental.  So far, this has only been tested using the NAG Fortran
compiler and the Apple Clang compilers.  This is a 3 step process

1. OpenMPI
2. TPL
3. Truchas

##### Notes for OpenMPI
There is an issue with failing to trigger the NAG preprocessor on the
case insensitive Mac filesystem.  This should be fixed for OpenMPI
versions > 4.0.3.  If compiling for version <= 4.0.3, the following
configuration step should suffice (from the top level source
directory):

	$ mkdir build
	$ cd build
	$ ../configure FC=nagfor FCFLAGS=-fpp --prefix=<my_openmpi_install_dir>

The `-fpp` option can be removed for later OpenMPI releases

##### Notes for TPL
Configuring cmake using the `config/linux-nag.cmake` is sufficent.  No other
special flags have been needed

##### Notes for Truchas
`PGSLib` is a subpackage of Truchas that depends heavily on the linux
shared library linking model.  The linux linking behavior can be
emulated with the `CMAKE_SHARED_LINKER_FLAGS` that are set in
`mac-nag-chk.cmake`.  With this in mind, the following allowed for
building truchas on a mac

	$ mkdir build
    $ cd build
    $ cmake -C ../config/mac-nag.cmake \
          -D CMAKE_BUILD_TYPE=Release \
	        -D TRUCHAS_TPL_DIR=<truchas_tpl_dir> ..
    $ make
    $ make install

### Testing
From the build directory run the command

    $ ctest

to run the regression test suite. On multi-core systems use the `-j<n>` option
to tell ctest how many processes it can run simultaneously; `-j8`, for example.
All tests should pass.

CMakeTestFortranCompiler.cmake is a modified version of the distributed
version that fixes an error (test program assumes implicit typing and so
fails when the "-u" compiler option is used).  This bug is fixed in 2.8.11.
This needs to be removed as soon as possible, and then the minimum version
level of cmake set to 2.8.11.

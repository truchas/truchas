The small utility program upgradevff.c can be used to upgrade an original
view factor file to work with the current version of Truchas and RadE tools.
What it does is add a new array "icount" to the netCDF file.  This new array
substitutes for the "ia" array which is no longer used.  See the comments at
the top of the source file.

### Compiling
Compiling upgradevff.c can be as easy as running the command

    $ gcc -std=c99 upgradevff.c -o upgradevff -lnetcdf

It only depends on the netcdf library and header files (any version).

### Usage
To upgrade the view factor file `vf.nc`, run the command

    $ upgradevff vf.nc

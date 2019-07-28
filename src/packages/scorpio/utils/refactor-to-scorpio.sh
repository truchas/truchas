#! /bin/bash

# Author: Sarat Sreepathi
# Date: March 05, 2013.
# This script is useful to port any applications using 
# earlier ASCEM-IO library calls to the new library interface.

# As executing some of these transformations twice may result in undesired changes, there are some checks 
# to see if running it is needed.

# MYFILES=*.[ch]
# MYFILES="*.c *.h *.F90 *.f90"
MYFILES=*.f90
sed -i 's/parallelIO/scorpio/g' $MYFILES; 
sed -i 's/ASCEMIO/SCORPIO/g' $MYFILES; 
sed -i 's/ASCEM-IO/SCORPIO/g' $MYFILES; 
sed -i 's/ascemio/scorpio/g' $MYFILES; 
sed -i 's/parallelio/scorpio/g' $MYFILES;
sed -i 's/PARALLELIO_LIB/SCORPIO/g' $MYFILES;

grep "\<PIO_INTEGER\>" $MYFILES
if [ $? -eq 0 ] ; then 
    echo "Fixing stuff"; 
	for i in PIO_INTEGER PIO_DOUBLE PIO_FLOAT PIO_LONG PIO_CHAR PIO_BYTE PIO_STRING
	do
	  # echo $i
	  sed -i "s:$i:${i//PIO/SCORPIO}:g" $MYFILES
	done
fi

grep "\<piof.h\>" $MYFILES
if [ $? -eq 0 ] ; then 
	sed -i 's/piof\.h/scorpiof\.h/g' $MYFILES;
fi

grep -E "\<FILE_READONLY\>|\<FILE_CREATE\>|\<FILE_READWRITE\>" $MYFILES
if [ $? -eq 0 ] ; then 
    echo "Fixing stuff"; 
	for i in NONUNIFORM_CONTIGUOUS_READ NONCONTIGUOUS_READ EVERYONE_ENTIRE_DATASET_READ NONUNIFORM_CONTIGUOUS_WRITE NONCONTIGUOUS_WRITE FILE_CREATE FILE_READONLY FILE_READWRITE UNIFORM_CONTIGUOUS_READ UNIFORM_CONTIGUOUS_WRITE 
	do
	  # echo $i
	  sed -i "s:\<$i\>:SCORPIO_$i:g" $MYFILES
	done
fi


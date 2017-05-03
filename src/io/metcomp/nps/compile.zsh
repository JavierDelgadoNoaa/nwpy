#!/usr/bin/env zsh
##
# Stand alone compilation script to build the library without distutils
##
set -aeu
DEST=../../lib/nps/
FFLAGS="-convert big_endian"
f2py -m -c --f90exec=`which ifort` --f90flags="$FFLAGS" write_nps_int write_nps_int.f90
mv write_nps_int.so $DEST

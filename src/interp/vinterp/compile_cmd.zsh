#!/usr/bin/env zsh
##
# Script to manually compile the interp_routines using f2py (as opposed
# to using distutils
##
set -aeu

DEST_DIR=../lib
# list of f2py opts: http://docs.scipy.org/doc/numpy-dev/f2py/usage.html#command-f2py
f2py -m -c --f90exec=`which ifort` interp_routines interp_routines.f90  
mv interp_routines.so $DEST_DIR

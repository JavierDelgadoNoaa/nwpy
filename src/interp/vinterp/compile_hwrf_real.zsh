#!/usr/bin/env zsh
##
# This script was used to compile the hwrf interpolation wrapper on NOAA's
# Theia, using the Intel Fortran Compiler
##


#DEST_DIR=../lib
DEST_DIR=../
PERFORMANCE_OPTS="-DF2PY_REPORT_ATEXIT -DF2PY_REPORT_ON_ARRAY_COPY=10"
DEBUG_OPTS=' --debug --f90flags="-g -traceback -O0 $PERFORMANCE_OPTS"'
[[ -z $1 ]] && echo -e "USAGE: $0 <module_name>.\n e.g. ``$0 foo'' will build "\
                       "foo.f90 and create foo.so in $DEST_DIR" && exit 1
MODULE=$1

set -aeu

# list of f2py opts: http://docs.scipy.org/doc/numpy-dev/f2py/usage.html#command-f2py

#f2py -m -c --f90exec=`which ifort` hwrf_real_interp hwrf_real_interp.f90
#f2py -m -c --f90exec=`which ifort` --debug hwrf_real_interp hwrf_real_interp.f90
#f2py -m -c --f90exec=`which ifort` --debug --f90flags="-g -traceback -O0" $MODULE $MODULE.f90
eval f2py -m -c --f90exec=`which ifort` $DEBUG_OPTS $MODULE $MODULE.f90
mv $MODULE.so $DEST_DIR

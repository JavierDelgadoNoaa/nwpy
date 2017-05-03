#!/usr/bin/env python
'''
Get the value of a given set of grib parameters from a particular
grib message in a given grib file.
By default, the 10-meter UGRD message is used, but this can be changed
in the script via 3 variables.
Note that this script is intended for extracting basic information
about the grid. For extracting large arrays (e.g. for plotting) see
the "plot_grib_fields.py" script.

OUTPUT: List of values for the requested parameters, separated by space

USAGE: get_grid_extents.py <file name> <list of parameters requested>
'''

import os
import sys
import pygrib


# 
# Configuration options
#
# grib message field from which to extract requested parameter
GRIB_FIELD_SHORT_NAME = '10u'
# type of level for the field from which to extract requested parameter
GRIB_LEVEL_TYPE = 'heightAboveGround'
# Level value ...
GRIB_LEVEL_VALUE = 10

#
# MAIN
#
if __name__ == '__main__':
    if len(sys.argv) > 2:
        file_name = sys.argv[1]
        params_requested = sys.argv[2:]
        if sys.argv[1].find("-h") > -1: 
            print __doc__
            sys.exit(0)
    else:
        #print 'USAGE: %s file params' %sys.argv[0]
        print __doc__
        sys.exit(1)
    if not os.path.exists(file_name):
        print 'File not found: ', file_name
        print __doc__
        sys.exit(1)

    grbindx = pygrib.index(file_name, 'shortName', 'typeOfLevel', 'level')
    selected_grbs = grbindx.select(shortName=GRIB_FIELD_SHORT_NAME, 
                                   typeOfLevel=GRIB_LEVEL_TYPE,
                                   level=GRIB_LEVEL_VALUE)
    if len(selected_grbs) == 0:
        print 'Did not find a u10 field in the given file'
        sys.exit(1)
    if len(selected_grbs) > 1:
        print 'More than one match, this is unexpected'
        sys.exit(1)
    grbmsg = selected_grbs[0]
    for param in params_requested:
        if grbmsg.has_key(param):
            sys.stdout.write("%s " %grbmsg[param])
        else:
            sys.stderr.write("The %s field of this GriB file does not contain the"
                             " field '%s'\n" %(GRIB_FIELD_SHORT_NAME, param))
            sys.stderr.write("The following fields are available:\n%s\n" %grbmsg.keys())
            sys.exit(2)
    sys.stdout.write("\n")

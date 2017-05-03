#!/usr/bin/env python
'''
This script allows you to plot one or more fields from a Grib file.
It can also be used to dump the contents of a grib file ala `wgrib -s`.

WARNING: If extents are relative to 360 in the Grib file,
         that is how they must be specified in Basemap instantiation
         Otherwise, contour/contourf (and probably other methods) will
         mask out all values thinking they are outside of the plotting area
         Currently this script uses gributils.basemap_from_grib(), 
         which uses get_grid_extents()
'''

import os
import sys
import logging as log
from optparse import OptionParser

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pygrib

import pycane.postproc.grib.util as gributils
from pycane.postproc.grib.objects import GribMessage
from nwpy.viz.map import bling

# basic constants - modify at your own risk 
CONTOUR = 3
FILLED_CONTOUR = 4
FIELD_PARAM_FORMAT = "(field shortName, typeOfLevel, level [, plotType]). "  \
                         "Paraentheses are optional. The plotType is optional, " \
                         "and its value can be either 'contour' (default) or 'filledContour'." \
                         "The shortName is the abbreviated field name."\
                         "e.g. \"(gh, isobaricInhPa, 450, contour) \" "
FIELDSET_HELP = "Specify fields to plot. Either a text file or a string"\
                " of fields may be given. If giving a text file, it should" \
                " have one field per line. If it is a string, each field"\
                " should be separated by a semicolon(;). The field should be"\
                " specified as follows: %s" %FIELD_PARAM_FORMAT 

# Field(s) to plot are specified using a 4-tupple consisting of the
#  ( shortName, level_type, level, plotType)
# - shortName is the abbreviated field name
# - level_type is hybrid, meanSea, isobaricInhPa, heightAboveGround, etc.
# - level is the value of the level
# - plotType is the kind of plot, which can be either CONTOUR or FILLED_CONTOUR
# indices in fields_to_plot tupples of different values
FIELD_NAME_IDX = 0
LEVEL_TYPE_IDX = 1
LEVEL_VALUE_IDX = 2
PLOT_TYPE_IDX = 3

#
# Defaults
#
# default log level, which can be overriden in the command line
DEFAULT_LOG_LEVEL = log.DEBUG
# padding at each edge of generated image, beyond the grib message extents
MAP_PADDING = 15 # degrees 

def _process_field_request_string(request_string):
    '''
    Given a string of 3 or 4 comma-separated tokens, optionally surrounded
    by parentheses, verify that the format is correct and if necessary
    set default value for the 4th token, which is optional and can be either
     'contour' or 'filledContour'
    RETURN - 4-tupple containing the verfied tokens
    '''
    toks = request_string.split(",")
    toks = [t.strip() for t in toks]
    if not len(toks) in (3,4):
        raise Exception("Unexpected input for --fieldset. Excpected format: %s"
                        %FIELD_PARAM_FORMAT)
    # remove leading/trailing parens
    if toks[0].startswith("("):
        toks[0] = toks[0][1:]
    if toks[-1].endswith(")"):
        toks[-1] = toks[-1][:-1]
    # set plot_type
    if len(toks) == 4:
        if toks[3] == 'contour':
            plot_type = CONTOUR
        elif toks[3] == 'filledContour':
            plot_type = FILLED_CONTOUR
        else:
            raise Exception("Plot type should be 'contour' or 'filledContour")
    else:
        plot_type = CONTOUR
    return  (toks[0],  toks[1],  int(toks[2]), plot_type)

def parse_cmdline():
    '''
    Parse command line arguments and options

    '''
    global fields_to_plot, input_file

    parser = OptionParser(usage="%prog [options] input_file. e.g. %prog -h")
    parser.add_option("-f", "--fieldset", dest="desired_fieldset",
                      help=FIELDSET_HELP)
    parser.add_option("-d", "--dump", dest="dump", 
                      action="store_true", default=False,
                      help="Dump inventory of the input grib file and exit."\
                           " Similar to wgrib -s output, but only show"\
                           " shortName, parameterName, levelType, level")
    parser.add_option("-l", "--log-level", dest="log_level", 
                      default=DEFAULT_LOG_LEVEL, type='int',
                      help='Specify the logging level. Should be an integer between' \
                           ' 0 and 50. 0 being very verbose, and 50 being the least')

    (options, args) = parser.parse_args()

    # Set the log level
    log.basicConfig(format='%(levelname)s : %(message)s', level=options.log_level)
    # set input file
    if len(args) < 1:
        parser.error("Input GriB file not specified")
    input_file = args[0]
    # just dump contents of input file, if --dump option passed
    if options.dump:
        gributils.inventory(input_file, 'shortName', 'parameterName', 'typeOfLevel', 'level')
        sys.exit(0)
    # process fieldset request 
    if options.desired_fieldset is None:
        raise Exception("Need to pass in fields to plot if not using --dump")
    if os.path.exists(options.desired_fieldset):
        log.debug("using input file for fields")
        with open(options.desired_fieldset) as fieldset_file:
            for line in fieldset_file.readlines():
                fields_to_plot.append(_process_field_request_string(line))
    else:    
        log.debug("Parsing fieldset given in command line")
        fieldsets = options.desired_fieldset.split(";")
        if len(fieldsets) == 0:
            raise Exception("Invalid input to --fieldset")
        for fs in fieldsets:
            fields_to_plot.append(_process_field_request_string(fs))
    if len(fields_to_plot) == 0:
        raise Exception("Either --fieldset not passed or format is wrong")
    log.debug("Will plot %i fields" %len(fields_to_plot))


#
# MAIN
#
if __name__ == '__main__':
    # fields to plot
    fields_to_plot = []
    # array to hold GribMessage objects that we are plotting
    grib_entries = []

    parse_cmdline() # logging will be set up here according to desired level

    # find desired grib message(s), create GribMessage from it, and add to 
    # grib_entries
    grbidx = pygrib.index(input_file, 'shortName', 'typeOfLevel', 'level')
    for fieldTupple in fields_to_plot:
        msg = grbidx.select(shortName = fieldTupple[FIELD_NAME_IDX],
                            typeOfLevel = fieldTupple[LEVEL_TYPE_IDX],
                            level = fieldTupple[LEVEL_VALUE_IDX])
        if len(msg) > 1 :
            log.info("Found more than one matching message for field = %s, " \
                     "level type = %s, level value = %i. Will use the first one" \
                     %(fieldTupple[FIELD_NAME_IDX], fieldTupple[LEVEL_TYPE_IDX],
                       fieldTupple[LEVEL_VALUE_IDX]))
        msg = msg[0]
        # Convert mslp values to make contour plots more readable
        if msg.shortName == 'msl':
            msg.values = msg.values / 100
        log.debug('Found matching grib entry:\t %s' %msg)
        lats, lons = msg.latlons()
        grib_entries.append( 
            GribMessage(msg.parameterName, msg.level, lats, lons, msg.values,
                        fieldTupple[PLOT_TYPE_IDX]))


    # plot the fields in grib_entries
    m = gributils.basemap_from_grib(msg, padding=MAP_PADDING, projection='cyl')
    for entry in grib_entries:
        if entry.plot_type == CONTOUR:
            entry.plot_contour_lines(m)
        elif entry.plot_type == FILLED_CONTOUR:
            entry.plot_contours(m)

    # final processing of map
    bling.decorate_map(m)
    #plt.show()
    plt.savefig('map.png')

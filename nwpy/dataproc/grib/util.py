#!/usr/bin/env python
'''
This module provides functions for extracting data from GriB files using PyGrib

Javier.Delgado@noaa.gov
'''

import os
import sys
import numpy
import pygrib
import logging as log
from mpl_toolkits.basemap import Basemap


def get_grid_extents_from_file(file_name, msgName='10u', 
                               msgLevelType='heightAboveGround',
                               msgLevel = 10):
    '''
    Get the lat and lon extents of a grib file. By default, use the "10u"
    message, i.e. the 10-meter UGRD. To get from a different message, adjust
    the `msgName', `msgLevelType', and `msgLevel'.
    The function accounts for scanning order of latitude points and for negative 
    values (the output will always show positive values)
    
    RETURNS
     a 4-tuple: (loLat, loLon, hiLat, hiLon)
    '''
    if not os.path.exists(file_name):
        print 'File not found: ', file_name
        sys.exit(1)
    
    grbindx = pygrib.index(file_name, 'shortName', 'typeOfLevel', 'level')
    selected_grbs = grbindx.select(shortName=msgName,
                                   typeOfLevel=msgLevelType, 
                                   level=msgLevel)
    if len(selected_grbs) == 0:
        print 'Did not find a u10 field in the given file'
        sys.exit(1)
    if len(selected_grbs) > 1:
        print 'More than one match, this is unexpected'
        sys.exit(1)
    grbmsg = selected_grbs[0]
    
    return get_grid_extents(grbmsg)

def get_grid_extents(grbmsg):    
    '''
    Get the lat and lon extents of a grib message.
    
    The function accounts for scanning order of latitude points and for negative 
    values (the output will always show positive values)
    
    Use the extents in the grib message fields if avaiable, 
    otherwise, estimate from latlons()
    
    RETURNS
     a 4-tuple: (loLat, loLon, hiLat, hiLon)
   '''
    use_latlons = False
    for k in ('latitudeOfFirstGridPointInDegrees',
              'longitudeOfFirstGridPointInDegrees',
              'latitudeOfLastGridPointInDegrees', 
              'longitudeOfLastGridPointInDegrees'):
        if not grbmsg.has_key(k):
            sys.stderr.write('Grib message does not contain necessary extent' \
                             'key. Using latlons()')
            use_latlons = True
    if not use_latlons:
        loLat = float(grbmsg['latitudeOfFirstGridPointInDegrees'])
        loLon = float(grbmsg['longitudeOfFirstGridPointInDegrees'])
        hiLat = float(grbmsg['latitudeOfLastGridPointInDegrees'])
        hiLon = float(grbmsg['longitudeOfLastGridPointInDegrees'])
        if int(grbmsg['jScansPositively']) == 0:
            tmp = loLat ; loLat = hiLat ; hiLat = tmp
    else:
        (lats,lons) = grbmsg.latlons() # latlons is [ 411 x [ 705 ] ] 
        loLat = numpy.amin(lats) # min(lats[-1])
        loLon = numpy.amin(lons) # max(lons[-1])
        hiLat = numpy.amax(lats) # max(lats[0]) 
        hiLon = numpy.amax(lons) # min(lons[-1])
    # Note : This is not working, but it is not a documented 
    # feature, so I'll leave it as is
    for attr in ('hiLat', 'hiLon', 'loLat', 'loLon'):
        if locals()[attr] < 0: locals()[attr] = 360 + locals()[attr]
    return (loLat, loLon, hiLat, hiLon)

def inventory(grib_file, *args):
    '''
    Dump the values of the GriB parameteres specified in *args for the given
    `grib_file'
    e.g. inventory('input.grib', 'parameterName', 'typeOfLevel', 'level')
    '''
    grbs = pygrib.open(grib_file)
    grbs.seek(0)
    outbuf = []
    for grb in grbs:
        outbuf.append("(")
        for i,param in enumerate(args):
            if i == len(args) - 1:
                outbuf.append("%s" %grb[param])
            else:
                outbuf.append("%s, " %grb[param])
        outbuf.append(")\n")
    print ''.join(outbuf)

def basemap_from_grib(grib_message, projection=None, resolution='l',
                      area_thresh=1000., northLat=None, southLat=None,
                      westLon=None, eastLon=None, padding=0):
    '''
    Create a BaseMap object using the given grib_message
    -If extents (northlat, southLat, westLon, eastLon) are not passed in, the 
     values from the grib files will be used.
    -padding specifies how much extra space to leave beyond the given/extracted
     extents, in degrees.
    The function has only been tested with cylindrical projection (i.e. cyl)      
    '''
    (lats, lons) = grib_message.latlons()
    (loLat, loLon, hiLat, hiLon) = get_grid_extents(grib_message)
    if northLat is None:
        northLat = hiLat
    if southLat is None:
        southLat = loLat
    if eastLon is None:
        eastLon = hiLon
    if westLon is None:
        westLon = loLon
    if projection is None:
        if grib_message['gridType'] == 'regular_ll':
            projection = 'cyl'
        else:
            log.warn("Setting projection to 'cyl' since nothing else has " \
                     "been tested")
            projection = 'cyl'
    log.debug("Creating basemap with the following extents (before padding of "\
              "%i), (westLon, eastLon) = (%f, %f); (southLat, northLat) = (%f,%f)"\
              %(padding, westLon, eastLon, southLat, northLat) )
    return Basemap(llcrnrlon=westLon-padding, llcrnrlat=southLat-padding,
                   urcrnrlon=eastLon+padding, urcrnrlat=northLat+padding, 
                   projection=projection, resolution=resolution,
                   area_thresh=area_thresh)

        
if __name__ == '__main__':
    # Test get_grid_extents_from_file() (which indirectly tests
    # get_grid_extents()
    if len(sys.argv) > 1:
        file_name = sys.argv[1]
    else:
        print 'USAGE: %s file' %sys.argv[0]
        sys.exit(1)
    #import pdb ; pdb.set_trace()
    (loLat, loLon, hiLat, hiLon) = get_grid_extents_from_file(file_name)
    print loLat, loLon, hiLat, hiLon

#!/usr/bin/env zsh
import os
import sys
import netCDF4 as nc4
import numpy as  np
import logging
from datetime import datetime as dtime
from datetime import timedelta as tdelta
    
# Inter mediate file version
VERSION = 5   
# Projection (only latlon is supported)
IPROJ = 0
# Radius of earth, in km
EARTH_RADIUS = 6367.47
# Reference point - change at your own risk
START_LOC = 'SWCORNER'
# Number assigned to missing vlaues in the nps_int format
NPS_INT_MISSING_VALUE = -1.E+30
# Maximum length (characters) of various NPS intermediate file attr names
MAX_NPS_FIELD_NAME_LENGTH = 9
MAX_NPS_UNITS_LENGTH = 26
MAX_NPS_DESCRIPTION_LENGTH = 46
MAX_NPS_OUTFILE_LENGTH = 300

def _default_log(log2stdout=logging.INFO, log2file=None):
    global _logger
    if _logger is None:
        _logger = logging.getLogger('nc2nps')
        _logger.setLevel(log2stdout)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        if log2file is not None:
            fh = logging.FileHandler('log.txt')
            fh.setLevel(log2file)
            fh.setFormatter(formatter)
            _logger.addHandler(fh)
        if log2stdout is not None:
            ch = logging.StreamHandler()
            ch.setLevel(log2stdout)
            ch.setFormatter(formatter)
            _logger.addHandler(ch)
    return _logger

def nc_to_nps_int(inFile, outFile, date, xfcst, fields, source=None, 
                  geos2wrf=False, log=None):
    '''
    Convert a netCDF4 file to NPS intermediate format. 
    @param inFile The input nc4 file name
    @param outFile The name of the File to output nps_int data to. It will be
                   __appended__ to.
    @param fields Array of 4-tupple consisting of: (0)variable name in the 
                            input file, (1) field name it should be given in 
                            the nps_int file, (2) units, (3) description
    @param date of the file, as datetime
    @param xfcst Forecast hour of the inFile
    @param source The model used to generate the file. This will help determine
                  the layout of the netCDF4 file. If not passed in, assume
                  the dimensions names are 'lev', 'lat', 'lon', 'time' and 
                  that lat and lon are 1-d. Known values are "g5nr" and "lis"
    @param geos2wrf True if output should be compatible with geos2wrf utilities.
                    This just changes the hdate to the format geos2wrf expects.
    '''               
    if log is None:
        log = _default_log()   
    flip_lats = False
    flip_lons = False
    # for each vertical level type in the netCDF file, map a standard
    # level ID (e.g. 'ps' for pressure) to its name in the netCDF file
    rootgrp_lev_types = {} # support multiple horizontal level types
    var2lev = {} # map 3d variables to lev_type
    if source == 'g5nr':
        (timeName,latName,lonName,rootgrp_lev_types['ps']) = ('time', 'lat', 'lon', 
                                                      'lev')
    elif source == 'lis':
        latName = 'north_south'
        lonName = 'east_west'
        rootgrp_lev_types['sm'] = 'SoilMoist_profiles'
        rootgrp_lev_types['st'] = 'SoilTemp_profiles'
        timeName = None
    else:
        (timeName,latName,lonName,rootgrp_lev_types['ps']) = ('time', 'lat', 'lon',
                                                      'lev')
    # Set basic attributes
    if geos2wrf:
        hdate = '{:%Y-%m-%d_%H}'.format(date)
    else:
        hdate = '{:%Y:%m:%d_%H:%M:%S}'.format(date)
        
    rootgrp = nc4.Dataset(inFile, 'r')
    
    # read the dimensions
    # hack! Estimate lat/lon for LIS
    # (TODO : Fix by flattening lat/lon to 1-d and accounting for 
    # the fact that lat/lon values are masked where there is no soil)
    # Actually, I don't think the nps_int file has a record of the lat/lon
    # values - it just uses the REF_LAT/REF_LON and DLAT/DLON, so we can
    # just use the attributes as already doing. The lat_var/lon_var are not
    # being used and the mask issue does not matter since we have the swCorner
    if source == 'lis':
        log.warn("Estimating lat/lon for LIS")
        swLat = rootgrp.getncattr("SOUTH_WEST_CORNER_LAT") 
        swLon = rootgrp.getncattr("SOUTH_WEST_CORNER_LON")
        dx = rootgrp.getncattr("DX") 
        dy = rootgrp.getncattr("DY")
        numLats = len(rootgrp.dimensions["north_south"])
        numLons = len(rootgrp.dimensions["east_west"])
        neLat = swLat + (numLats * dx) 
        neLon = swLon + (numLons * dy)
        lat_var = np.linspace(swLat, neLat, numLats)
        lon_var = np.linspace(swLon, neLon, numLons)
        # intermediate format wants west->east and south->north
        flip_lats = True
        flip_lons = True
    else:
        lat_var = rootgrp.variables[latName]
        lon_var = rootgrp.variables[lonName]
        if lat_arr[0] > lat_var[1]:
            log.info("Flipping latitude values to go South->North")
            flip_lats = True
            lat_var[:] = lat_var[::-1]
        if lon_var[0] > lon_var[1]:
            log.debug("Flipping longitude values to go West->East")
            flip_lons = True
            lon_var[:] = lon_var[::-1]
        deltalat = ( lat_var[1] - lat_var[0] )
        deltalon = ( lon_var[1] - lon_var[0] ) 
        dx = 110.0 * deltalat
        dy = 110.0 * deltalon
        
    
    # read the variables
    for (inName,outName,units,description) in fields:
        # sanity checks
        assert len(outName) <= MAX_NPS_FIELD_NAME_LENGTH
        assert len(units) <= MAX_NPS_UNITS_LENGTH
        assert len(description) <= MAX_NPS_DESCRIPTION_LENGTH
        assert len(outFile) <= MAX_NPS_OUTFILE_LENGTH
        log.debug("Processing geos5 variable '{}'".format(inName))
        var = rootgrp.variables[inName]
        # figure out if it is 2d or 3d
        is_3d = False
        for levType,levName in rootgrp_lev_types.iteritems():
            if levName in var.dimensions:
                is_3d = True
                log.debug("Treating variable '{}' as 3D".format(inName))
                # now know level type for this variable is `levType'
        # process
        if not is_3d:
            # NOTE : The slab should be a 2d variable with lon being the first
            # dimension (on the fortran side)
            dimNames = (timeName, latName, lonName, None)
            slab = get_2d_slab_from_var(var[:], dimNames, None, 
                                        flipLats=flip_lats,
                                        flipLons=flip_lons)
            xlvl = 200100.000000
                
            # set missing values - TODO this is SLOW, use Fortran
            slab[np.where(slab[:] == var.missing_value)] = NPS_INT_MISSING_VALUE
            write_2d_unformatted(outFile, VERSION, hdate, xfcst, source,
                                 outName, units, description, xlvl, numLons, 
                                 numLats, START_LOC, IPROJ, lat_var[0], 
                                 lon_var[0], dx, dy,
                                 deltalat, deltalon, EARTH_RADIUS, 0.0, 0.0,
                                 0.0, numLats, False, slab)
        else:
            dimNames = (timeName, latName, lonName, levName)
            log.info("For soil params, assuming we start at surface")
            curr_start_depth = 0.
            levIdx = var.dimensions.index(levName)
            for levCtr in range(1, var.shape[levIdx]+1):
                slab = get_2d_slab_from_var(var[:], dimNames, lev=levCtr, 
                                            flipLats=flip_lats,
                                            flipLons=flip_lons)
                # set missing values - TODO this is SLOW, use Fortran
                missingIdc = np.where(slab[:] == var.missing_value)
                slab[missingIdc] = NPS_INT_MISSING_VALUE
                # Set xlvl and outName (if necessary) according to levType
                if levType in ('sm', 'st'):
                    # soil moisture/temperature level - need to change 
                    # outName according to depth range
                    # This only works for LIS, AFAIK
                    xlvl = 200100.000000
                    thicknesses = rootgrp.getncattr('SOIL_LAYER_THICKNESSES')
                    thicknesses = [ v.round() for v in thicknesses ]
                    if thicknesses != [10.0, 30.0, 60.0, 100.0]:
                        log.warn("Unexpected thicknesses: {},{},{},{}"
                                .format(thicknesses))
                    curr_end_depth = curr_start_depth + thicknesses[levCtr]
                    pfx = levType.upper()
                    outName = nps_utils.get_nps_soil_field_name(
                        pfx, curr_start_depth, curr_end_depth  )
                    curr_start_depth = curr_end_depth
                elif levType == 'ps':
                    # pressure level meteorological variable 
                    #xlvl = rootgrp_lev_types[levType].levIdx
                    log.warn("Just putting indices for 'lev' ala NPS.")
                    xlvl = levCtr
                else:
                    raise Exception("Unknown height/level dimension type")
                
                write_2d_unformatted(outFile, VERSION, hdate, xfcst, source,
                        outName, units, desc, xlvl, numLons, numLats,
                        IPROJ, START_LOC, lat_var[0], lon_var[0], dx, dy,
                        deltalat, deltalon, EARTH_RADIUS, 0.0, 0.0,
                        0.0, numLats, False, slab[:])

        
def get_2d_slab_from_var(var, dimNames, lev=None, flipLats=False, 
                           flipLons=False):
    '''
    Get a 2D slab of data from NC variable `var' for level `lev'.
    `var' may have a scalar timeName dimension, in which case it will be 
    sliced out. But timeName must be the first dimension in this case.
    This function will return a 
    2-D slab of data. The data will be ordered in (lat,lon), since the 
    write_2d_unformatted routine expects the slab in (lon,lat) in 
    Fortran ordering (i.e. column major).
    @param var The netCDF variable to get the slab from
    @param dimNames 3-tupple of (timeName or None, latname, lonName) with the
                    names given to the dimension variables in var
    @param lev index of level number to extract. Ignored for 2d vars
    @param flipLats True if latitude values need to be reversed
    @param flipLons True if longitude values are to be reversed.
    '''
    msg_unsupported_fmt = "Unsupported netCDF variable dimension ordering"
    (timeName, latName, lonName, levName) = dimNames
    var_dims = var.dimensions 
    latIdx = var_dims.index(latName)
    lonIdx = var_dims.index(lonName)
    assert abs(latIdx - lonIdx) == 1 # must be adjacent
    if timeName is not None and timeName in var_dims:
        if var_dims.index(timeName) == 0:
            has_time_dim = True
        else:
            msg = "{}. {} dimension must be first".format(msg_unsupported_fmt, 
                                                    timeName) 
            raise Exception(msg)
    if levName is not None and levName in var_dims:
        log.debug("'{}' will be treated as 3-d variable".format(var.long_name))
        if has_time_dim:
            if var_dims.index(levName) == 1:
                slab2d = var[0,lev,:,:]    
            elif var_dims.index(levName) == len(var_dims) - 1:
                slab2d = var[0,:,:,lev]
            else:
                raise Exception(msg_unsupported_fmt)                
        else: # no time dimension
            if var_dims.index(levName) == 0:
                slab2d = var[lev,:,:]
            elif var_dims.index(levName) == len(var_dims) - 1:
                slab2d = var[:,:,lev]
            else:
                raise Exception(msg_unsupported_fmt)
    else:
        log.debug("'{}' will be treaded as 2-d variable".format(var.long_name))
        if has_time_dim:
            slab2d = var[0,:]
        else:
            slab2d = var[:]
    if latIdx > lonIdx:
        log.info("Lon comes before Lat. Will transpose")
        slab2d = np.transpose(slab2d)
    if flipLats:
        varOut = varOut[::-1,:]
    if flipLons:
        varOut = varOut[:,::-1]
    
    return slab2d    
 
                
def g5nr_to_nps_int(filename):
    nc_to_nps_int(filename, source='g5nr')
    field = (inName,outName,units,description)
   
def main():
    start_date = dtime(year=2006, month=9, day=10, hour=0)
    duration = tdelta(hours=0)
    frequency = tdelta(hours=3)
    log = _default_log(log2stdout=logging.DEBUG)
    
    # general
    date = start_date
    xfcst = 0
    outFile = 'FILE:foobar'
    
    #file-specific - MET
    inFile = 'g5nr_test_2_metfields_hires.20060910_0000z.nc4'
    #Field = (inName,outName,units,description)
    tField = ('T','TT','K','air temperature')
    slpField = ('SLP', 'PMSL', 'Pa', 'sea level pressure')
    fields = [tField, slpField]
    nc_to_nps_int(inFile, outFile, date, xfcst, fields, source='g5nr', 
                  geos2wrf=False)
        
    # file-specific - Soil
    inFile = 'LIS_HIST'
    stField = ('SoilTemp_tavg', 'ST', 'K', 'soil temperature')
    fields = [stField]
    nc_to_nps_int(inFile, outFile, date, xfcst, fields, source='g5nr', 
                  geos2wrf=False)
    
    # Think: idemponence, since this will likely be very time consuming
    
    # TODO  ? 
    # 0. Create fortran interface
    #
    # 1. Should the xfcst be the actual g5nr forecast hour or should it be 
    # relative to when we are initializing the regional forecast? With WRF we
    # usually use the GFS forecast that starts at the same time, so I don't
    # know if there is an assumption here.
    #
    # 2. Verify date against netCDF file.
    #  For LIS: begin_time (HHMMSS) and begin_date (ymd)
    
"""
This module contains various functions to faciliate use of and interaction 
with NPS intermediate format files.
Key functionalities:
 -Functions for converting netCDF data to NPS-intermediate format
ASSUMPTIONS
 - input (netCDF) data is on a latlon projection
"""
import os
import sys
import logging
from datetime import datetime as dtime
from datetime import timedelta as tdelta

import netCDF4 as nc4
import numpy as  np
from cfunits import Units

from write_nps_int import write_2d_slab
import nps_utils

# Projection (only latlon is supported)
IPROJ = 0
# Radius of earth, in km
EARTH_RADIUS = 6367.47
# Reference point - change at your own risk
START_LOC = 'SWCORNER'
# Number assigned to missing vlaues in the nps_int format
#NPS_INT_MISSING_VALUE = -1.E+30
NPS_INT_MISSING_VALUE = 0.0
# Maximum length (characters) of various NPS intermediate file attr names
MAX_NPS_FIELD_NAME_LENGTH = 9
MAX_NPS_UNITS_LENGTH = 26
MAX_NPS_DESCRIPTION_LENGTH = 46
MAX_NPS_OUTFILE_LENGTH = 300

# Assumed soil thicknesses
ASSUMED_SOIL_THICKNESSES = [10.0, 30.0, 60.0, 100.0]                  

# Set internal vars
_logger = None

# GLOBAL VARS
DEFAULT_EXPECTED_UNITS = \
  {
    "TT"     :    "K"  ,
    "RH"     :    "%"  ,
    "SPECHUMD"     :    "kg kg-1"  , # This is what it is in GFS files
    "UU"     :    "m s-1"  ,
    "VV"     :    "m s-1"  ,
    "HGT"     :    "m"  ,
    "PRESSURE"     :    "Pa"  ,
    "PINT"     :    "Pa"  ,
    "PSFC"     :    "Pa"  ,
    "PMSL"     :    "Pa"  ,
    "CLWMR"    :    "kg kg-1",
    "SM000010"     :    "proprtn", # Note that this is incorrectly specified
    "SM010040"     :    "proprtn", # as kg/m2 in the Vtable, although the
    "SM040100"     :    "proprtn", # parameters correspond to SOILW (volumetric
    "SM100200"     :    "proprtn", # soil moisture) in current GriB spec.
    "SM010200"     :    "proprtn",
    "ST000010"     :    "K"  ,
    "ST010040"     :    "K"  ,
    "ST040100"     :    "K"  ,
    "ST100200"     :    "K"  ,
    "ST010200"     :    "K"  ,
    "ST000010"     :    "K"  ,
    "ST010040"     :    "K"  ,
    "ST040100"     :    "K"  ,
    "ST100200"     :    "K"  ,
    "ST010200"     :    "K"  ,
    "SEAICE"     :    "proprtn"  ,
    "LANDSEA"     :    "proprtn"  ,
    "SOILHGT"     :    "m"  ,
    "SKINTEMP"     :    "K"  ,
    "SNOW"     :    "kg m-2"  ,
    "LAND_LIS" : "proprtn",  # LANDSEA field created from LIS data
  }

def create_seaice(frLakeVar, frOceanVar, outFile):
    '''
    Create an NPS intermediate format file from netCDF variables
    containing fraction of lake and ocean.
    @param frLakeVar netCDF4.Variable containing fraction of lake
    @param frOceanVar netCDF4.Variable containing fraction of ocean
    @param outFile String containing path to store output at
    '''

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

def get_int_file_name(source, date, includeMMSS=False):
    """
    Returns the intermediate file name, in typical format, given the variable
    source and date
    @param source Source the variable came from (e.g. 'GFS' or just 'FILE' to 
                   use default
    @param date datetime.datetime object of the date
    @param includeMMSS If True, include the minutes and seconds as well 
    @return a string <source>:YYYY-MM-DD_HH or <source>:YYYY-MM-DD_HH:MM:SS
    """
    if includeMMSS:
        return "{}:{:%Y-%m-%d_%H:%M:%S}".format(source, date)
    else:
        return "{}:{:%Y-%m-%d_%H}".format(source, date)

def write_slab_cyl(slab, outFile, hdate, xfcst, source,
                   outName, units, description, xlvl,  
                   startLat, startLon, deltalat, deltalon,
                   earthWindRelative=False, altOutFile=None ):
    '''
    Write 2-D slab of data in cylindrical equidistant coordeinates
    to intermediate file.
    @param slab The 2-d slab of data arranged in (lon,lat) ordering. f2py
                will internally transpose it when calling Fortran
    @param outFile String of output file path. It will be _appended_ to
    @param hdate String of forecast data in the NPS format, 
           e.g. 2006:09:10_00:00:00
    @param xfcst Floating point forecast hour
    @param source The source the data comes from (e.g. 'G5NR')
    @param outName name to give the field/variable in the output file
    @param units K, m/s, etc.
    @param description Up to 46 characters describing the field
    @param xlvl The level number to pass to the intermediate file
    @param startLat Starting latitude 
    @param startLon Starting longitude
    @param deltalat Grid spacing, in degrees
    @param deltalon Grid spacing, in degrees
    @param altOutFile If given, also write the slab to this file
    '''
    # sanity checks
    assert len(outName) <= MAX_NPS_FIELD_NAME_LENGTH
    assert len(units) <= MAX_NPS_UNITS_LENGTH
    assert len(description) <= MAX_NPS_DESCRIPTION_LENGTH
    assert len(outFile) <= MAX_NPS_OUTFILE_LENGTH
    #(nx,ny) = slab.shape
    #print 'jza (v) ny = ', ny, 'nx = ', nx
    iproj = 0
    # thse are not used for latlon grids in version 5
    # of the intermediate file format.
    truelat1 = truelat2 = xlonc = dx = dy =numLats = numLons = 0.0
    version = 5
    # pad strings to be expected size
    param_lengths = { 'hdate':24, 'source':32, 'outName':9, 'units':25, 
                      'description':46 } # TODO - use this
    #hdate = hdate.rjust(24)
    #source = source.rjust(32)
    #units = units.rjust(25)
    #description = description.rjust(46)
    #outName = outName.rjust(9)
    start_loc = START_LOC.rjust(8)
    xlvl = float(xlvl)
    print 'el_slab.shape: ', slab.shape
    write_2d_slab(outFile, version, hdate, xfcst, source,
                  outName, units, description, xlvl,  
                  iproj, start_loc, startLat, startLon, dx, dy,
                  deltalat, deltalon, EARTH_RADIUS, truelat1, truelat2,
                  xlonc, numLats, earthWindRelative, slab)
    if altOutFile is not None:
        try:
            os.mkdir(os.path.dirname(altOutFile))
        except OSError:
            pass
	    write_2d_slab(altOutFile, version, hdate, xfcst, source,
                  outName, units, description, xlvl,  
                  iproj, start_loc, startLat, startLon, dx, dy,
                  deltalat, deltalon, EARTH_RADIUS, truelat1, truelat2,
                  xlonc, numLats, earthWindRelative, slab)

def _get_alt_out_file_path(createIndividualFiles, outDir, varName, lev, date):
    """
    Get name of alternative output file. 
    If `createIndividualFiles' is False, return None
    @param createIndividualFiles Do we want to create the individual files?
    @param outDir Path to put files in
    @param varName Name of the variable in NPS
    @param lev The level number/id
    @param date Datetime object representing the date
    """
    if createIndividualFiles:
        yyyymmddHHMM = '{:%Y%m%d%H%M}'.format(date)
        fil = "{}_{}_{}.nps".format(varName, lev, yyyymmddHHMM)
        altOutFile = os.path.join(outDir, 'npsIntByField', fil)
    else:
        altOutFile = None
    return altOutFile

def __get_expected_units(file_path):
    """
    Read a file containing key, value pairs separated by ":" and return a 
    dictionary
    """
    d = {}
    with open(file_path) as f:
        for line in f.readlines():
            toks = line.split(":")
            d[toks[0].strip()] = toks[1].strip()

def __verify_units(expectedUnitsFile, var, srcVarName, destVarName, givenUnits, log):
    """
    Verify units in input against what is expected by NPS.
    If the nc4 variable `var' does not have a units attribute,
    use givenUnits.
    srcVarName is the name of the variable in the source dataset
    destVarName is the name to be used for input to Metgrid and NemsInterp.
    Return 2-tuple: (in_units, expected_units). in_units will be either
           the Variable's "units" attribute or givenUnits, with spaces
           stripped and possibly transformed to a value that cfunits understands.
           expected_units will be the expected units, with similar transformations.
    """
    if expectedUnitsFile is None:
        expected_units = DEFAULT_EXPECTED_UNITS[destVarName]
    else:
        expected_units = __get_expected_units(expectedUnitsFile)[destVarName]
    try:
        dsUnits = getattr(var, "units")
        if dsUnits != givenUnits:
            log.warn("Units read from dataset are '{inp}', but specified"
                     " units are '{out}'. Will use dataset values, since"
                     " that is what we are copying"
                     .format(inp=givenUnits, out=dsUnits))
            in_units = dsUnits
        else:
            in_units = givenUnits
    except AttributeError:
        log.info("Variable '{v}' does not have a 'units' attribute, "
                 "will use the passed-in value '{u}'"
                 .format(v=srcVarName, u=givenUnits))
        in_units = givenUnits
    try:
        # strip out spacing since it's not consistent
       # expected_units = expected_units.replace(" ","") 
       # in_units = in_units.replace(" ","")
        
        # Replace values used by GriB/NPS that do not have mapping in 
        if expected_units in [1, "1", "proprtn", "proportion", "proportn"]:
            expected_units = "1"
        if in_units in [1, "1", "proprtn", "proportion", "proportn"]:
            in_units = "1"
        
        #if in_units != expected_units:
        # strip out spacing since it's not consistent
        if in_units.replace(" ","") != expected_units.replace(" ",""):
            log.info("Variable {v}: Dataset/given units = {inp}, but  "
                     "expected units = {exp}. Conversion required"
                     .format(v=srcVarName, inp=givenUnits, exp=expected_units))
            # NOTE: conversion will be done in get_2d_slab_from_var                      

    except KeyError:
        log.critical("The expected_units dict does not have a mapping "
                     "for variable '{0}'".format(destVarName))
        sys.exit(44)
    
    return in_units, expected_units

def nc_to_nps_int(inFile, outFile, date, xfcst, fields, source=None, 
                  geos2wrf=False, log=None, createIndividualFiles=False,
                  expectedUnitsFile=None):
    """
    Convert a netCDF4 file to NPS intermediate format.
    ASIDE: This function was verified by comparing outputs generated to those
           generated using Gus Alaka's unpack_netcdf program. The generated 
           files were identical for TT. Note that some basic changes to the
           code were necessary for this (e.g. levels had to be read in reverse
           order and missing value was different), so we cannot expect a 
           an identical match with the current implementation.
            -Performance note:
                python nc2nps.py  272.64s user 38.81s system 91% cpu 5:39.54 total
                ./unpack_netcdf.exe  118.90s user 9.51s system 83% cpu 2:33.32 total
                # removing the code that looks for masked values did not help much:
                python nc2nps.py  266.17s user 25.14s system 93% cpu 5:11.67 total

    @param inFile The input nc4 file name
    @param outFile The name of the File to output nps_int data to. It will be
                   __appended__ to.
    @param fields Array of field data for the fields to be processed. 
                    Each element is a 4-tuple (inName, outName, units, desc):
                            (0) variable name in ``inFile'' 
                            (1) name it should be given in the nps_int file
                            (2) units
                            (3) description
                    The "units" will only be used if the ``inFile'' does not have a 
                    "units" attribute for the variable "inName". The units will
                    be compared to the expected units for "outName". The expected
                    units are specified via argument ``expectedUnits''.
                    *NOTE: It's very important to ensure the units are correct
                     here (i.e. what Metgrid and/or NemsInterp/Real) since those
                     programs do not make careful verification of the units. If 
                     the units to not match what is expected, this function will
                     convert them using cfunits
    @param date of the file, as datetime
    @param xfcst Forecast hour of the inFile
    @param source The model used to generate the file. This will help determine
                  the layout of the netCDF4 file. If not passed in, assume
                  the dimensions names are 'lev', 'lat', 'lon', 'time' and 
                  that lat and lon are 1-d. Known values are "g5nr" and "lis"
    @param geos2wrf True if output should be compatible with geos2wrf utilities.
                    This just changes the hdate to the format geos2wrf expects.
    @param createIndividualFiles If True, create separate individual nps_int 
            files (e.g. for debugging). The files will be created in the same
            directory as outFile, using the following naming convention:
            <outPath>/npsIntByField/<field>_<xlvl>_<yyyymmddHHMM>.nps
    @param expectedUnitsFile Path to file containing expected units. The file
           should have a line for each variable. Each line should consist of
           "varName : units", where ``varName'' is the `outVarName'', i.e. the
           second element of the tupples of the given ``fields'' array.
           If None, use DEFAULT_EXPECTED_UNITS, which has the mappings used
           in the Vtable for NPS as of 8/2016.
    """               
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
    log.debug("Reading file {}".format(inFile))    
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
        deltalon = rootgrp.getncattr("DX") 
        deltalat = rootgrp.getncattr("DY")
        numLats = len(rootgrp.dimensions["north_south"])
        numLons = len(rootgrp.dimensions["east_west"])
        neLat = swLat + (numLats * deltalat) 
        neLon = swLon + (numLons * deltalon)
        lat_var = np.linspace(swLat, neLat, numLats)
        lon_var = np.linspace(swLon, neLon, numLons)
        # intermediate format wants west->east and south->north
        flip_lats = True
        flip_lons = True
        dx = 110.0 * deltalon
        dy = 110.0 * deltalat
    else:
        lat_var = rootgrp.variables[latName]
        lon_var = rootgrp.variables[lonName]
        if lat_var[0] > lat_var[1]:
            log.info("Flipping latitude values to go South->North")
            flip_lats = True
            lat_var[:] = lat_var[::-1]
        if lon_var[0] > lon_var[1]:
            log.debug("Flipping longitude values to go West->East")
            flip_lons = True
            lon_var[:] = lon_var[::-1]
        deltalat = ( lat_var[1] - lat_var[0] )
        deltalon = ( lon_var[1] - lon_var[0] ) 
        dx = 110.0 * deltalon
        dy = 110.0 * deltalat
        
    
    # read the variables
    for (inName,outName,inUnits,description) in fields:
        log.debug("Processing {} variable '{}'".format(source, inName))
        #var = rootgrp.variables[inName]
        # figure out if it is 2d or 3d
        # hack - only do this for met fields since the variable name
        # passed in for LSM variables is not the actual variable name 
        # and we know that they are 3d
        if inName in ('SM', 'SoilMoist_tavg'):
            is_3d = True
            levType = 'sm'
            levName = rootgrp_lev_types['sm']
            log.warn("Reading 'SoilMoist_tavg' instead of passed in {}".format(inName))
            var = rootgrp.variables['SoilMoist_tavg']
            varForUnitsHack = "SM010200" # hack: Need somthing that's in expected_units
        elif inName in ('ST', 'SoilTemp_tavg'):
            is_3d = True
            levType = 'st'
            levName = rootgrp_lev_types['st']
            log.warn("Reading 'SoilTemp_tavg' instead of passed in {}".format(inName))
            var = rootgrp.variables['SoilTemp_tavg']
            #import pdb ; pdb.set_trace()
            varForUnitsHack = "ST010200" # hack: need something that's in expected_units
        else:
            is_3d = False # changed below if 3d
            try:
                var = rootgrp.variables[inName]
            except KeyError:
                log.critical("Variable {var} is not in dataset {inFile}"
                             .format(var=inName, inFile=inFile))
                sys.exit(1)
            for levType,levName in rootgrp_lev_types.iteritems():
                if levName in var.dimensions:
                    is_3d = True
                    log.debug("Treating variable '{}' as 3D".format(inName))
                    # now know level type for this variable is `levType'
            varForUnitsHack = outName

        (inUnits, out_units) = __verify_units(expectedUnitsFile, var, 
                                              #inName, outName, inUnits, log)
                                              inName, varForUnitsHack, inUnits, log)

        # process
        if not is_3d:
            # NOTE : The slab should be a 2d variable with lon being the first
            # dimension (on the fortran side)
            dimNames = (timeName, latName, lonName, None)
            slab = get_2d_slab_from_var(var, dimNames, None, 
                                        inUnits=inUnits, outUnits=out_units,
                                        flipLats=flip_lats,
                                        flipLons=flip_lons, log=log)
            xlvl = 200100.000000
                
            # set missing values - TODO this is SLOW, use Fortran
            try:
                slab[np.where(slab[:] == var.missing_value)] = NPS_INT_MISSING_VALUE
            except AttributeError:
                log.warn("Variable '{0}' does not have a 'missing_value' "
                         "attribute; unable to set the NPS_INT_MISSING_VALUE"
                         .format(inName))

            altOutFile = _get_alt_out_file_path(createIndividualFiles, 
                                                os.path.dirname(outFile),
                                                outName, 200100, date)
            #import pdb ; pdb.set_trace()
            write_slab_cyl(slab, outFile, hdate, xfcst, source, outName,  
                           out_units, description, xlvl, lat_var[0], lon_var[0], 
                           deltalat, deltalon, altOutFile=altOutFile)
        else: 
            # 3d field
            dimNames = (timeName, latName, lonName, levName)
            log.info("For soil params, assuming we start at surface")
            curr_start_depth = 0.
            levIdx = var.dimensions.index(levName)
            #for levCtr in range(1, var.shape[levIdx]+1):
            #for levCtr in range(var.shape[levIdx]-1, -1, -1):
            for levCtr in range(var.shape[levIdx]):
                slab = get_2d_slab_from_var(var, dimNames, lev=levCtr, 
                                            flipLats=flip_lats,
                                            inUnits=inUnits, outUnits=out_units,
                                            flipLons=flip_lons, log=log)
                # set missing values - This is a bit SLOW, but not a bottleneck
                # TODO : Works for LIS. Ensure this works for g5nr data too.
                #import pdb ; pdb.set_trace()
                if isinstance(slab, np.ma.masked_array):
                    missingIdc = np.where(slab.mask == True)
                else:
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
                    if thicknesses != ASSUMED_SOIL_THICKNESSES:
                        log.warn("Unexpected thicknesses: {},{},{},{}"
                                .format(thicknesses))
                    curr_end_depth = curr_start_depth + thicknesses[levCtr]
                    pfx = levType.upper()
                    log.info("Overriding variable name for soil moist./temp.")
                    outName = nps_utils.get_nps_soil_field_name(
                        pfx, int(curr_start_depth), int(curr_end_depth)  )
                    log.info("Overriding description for soil moist./temp.")
                    description = nps_utils.get_nps_soil_field_description(
                        pfx, int(curr_start_depth), int(curr_end_depth)   )
                    curr_start_depth = curr_end_depth
                elif levType == 'ps':
                    # pressure level meteorological variable 
                    #xlvl = rootgrp_lev_types[levType].levIdx
                    msg = "Just putting indices for 'lev' ala NPS."
                    if not msg in __already_logged:
                        log.warn(msg)
                        __already_logged.append(msg)
                    xlvl = levCtr + 1 # fortran
                else:
                    raise Exception("Unknown height/level dimension type")
                
                altOutFile = _get_alt_out_file_path(createIndividualFiles, 
                                                    os.path.dirname(outFile),
                                                    outName, xlvl, date)
                
                write_slab_cyl(slab, outFile, hdate, xfcst, source,
                        outName, out_units, description, xlvl, 
                        lat_var[0], lon_var[0], deltalat, deltalon, 
                        altOutFile=altOutFile)

        
def get_2d_slab_from_var(var, dimNames, lev=None, inUnits=None, outUnits=None,
                         flipLats=False, flipLons=False, log=None):
    '''
    Get a 2D slab of data from NC variable `var' for level `lev'.
    `var' may have a scalar timeName dimension, in which case it will be 
    sliced out. But timeName must be the first dimension in this case.
    This function will return a 
    2-D slab of data. The data will be ordered in (lon,lat), since the 
    Fortran routine expects this order and f2py will automatically copy the 
    array and transpose it (to column major) when calling Fortran routines.
    @param var The netCDF variable to get the slab from
    @param dimNames 3-tupple of (timeName or None, latname, lonName) with the
                    names given to the dimension variables in var
    @param lev index of level number to extract. Ignored for 2d vars
    @param flipLats True if latitude values need to be reversed
    @param flipLons True if longitude values are to be reversed.
    '''
    if log is None: log = _default_log()
    msg_unsupported_fmt = "Unsupported netCDF variable dimension ordering"
    (timeName, latName, lonName, levName) = dimNames
    var_dims = var.dimensions 
    latIdx = var_dims.index(latName)
    lonIdx = var_dims.index(lonName)
    assert abs(latIdx - lonIdx) == 1 # must be adjacent
    has_time_dim = False
    if timeName is not None and timeName in var_dims:
        log.debug("Field has 'time'")
        if var_dims.index(timeName) == 0:
            has_time_dim = True
        else:
            msg = "{}. {} dimension must be first".format(msg_unsupported_fmt, 
                                                    timeName) 
            raise Exception(msg)
    if levName is not None and levName in var_dims:
        log.debug("Processing level {} , field '{}'".format(lev, var.long_name))
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
        try:
            log.debug("Var '{}' being treated as 2-d".format(var.long_name))
        except AttributeError:
            log.debug("var will be treated as 2-d variable")
        if has_time_dim:
            slab2d = var[0,:]
        else:
            slab2d = var[:]
    if lonIdx > latIdx:
        msg = "Variable {v}: Lat comes before Lon. Will transpose. Supressing "\
              "further messages about this.".format(v=var.name)
        if not msg in __already_logged:
            __already_logged.append(msg)
            log.info(msg)
        slab2d = np.transpose(slab2d)
    if flipLats:
        slab2d = slab2d[::-1,:]
    if flipLons:
        slab2d = slab2d[:,::-1]
    # Convert units if necessary
    if inUnits != outUnits:
        assert not None in [inUnits, outUnits]
        msg = "Variable {v}: Conforming units from {iu} to {ou}. Supressing "\
              " further messages"\
              .format(v=var.name, iu=inUnits, ou=outUnits)
        if not msg in __already_logged:
            log.info(msg)
            # NOTE: It's appended to __already_logged after the "finished" msg
        # For soil moisture, gotta do the conversion manually
        if inUnits == "kg m-2" and outUnits in ("proportion", "proprtn",1, "1", 
                                                "proportn"):
            log.info("Converting '{0}' data from 'kg/m2' to ratio"
                     .format(var.long_name))
            # src: http://ldas.gsfc.nasa.gov/faq/#SoilMoist
            # Note we verify these thicknesses are the same as specified in 
            # the SOIL_LAYER_THICKNESSES attr in nc_to_nps_int()
            thickness = ASSUMED_SOIL_THICKNESSES[lev]
            slab2d /= (thickness * 10.) # cm -> mm
            log.debug("slab2d (min,max) = ({0},{1})".format(slab2d.min(),slab2d.max()))

        else:
            slab2d = Units.conform(slab2d, Units(inUnits), Units(outUnits))
        if not msg in __already_logged:
            log.debug("Finished conforming units")
            __already_logged.append(msg)
    return slab2d    
__already_logged = [] 
                
def g5nr_to_nps_int(filename, log=None):
    if log is None: log = _default_log()
    nc_to_nps_int(filename, source='g5nr', log=log)
    field = (inName,outName,units,description)
   
def main():
    start_date = dtime(year=2006, month=9, day=10, hour=0)
    duration = tdelta(hours=0)
    frequency = tdelta(hours=3)
    log = _default_log(log2stdout=logging.DEBUG)
    
    # general
    date = start_date
    xfcst = 0
    outFile = os.path.join('alfa', 'FILE:foobar') # TODO 
    
    #file-specific - MET
    inFile = 'g5nr_test_2_metfields_hires.20060910_0000z.nc4'
    inFile = os.path.join('/home/Javier.Delgado/scratch/nems/g5nr/data', inFile)
    #Field = (inName,outName,units,description)
    tField = ('TT','TT','K','air temperature')
    #slpField = ('SLP', 'PMSL', 'Pa', 'sea level pressure') # combined nc file has NPS fieldName
    slpField = ('PMSL', 'PMSL', 'Pa', 'sea level pressure')
    fields = [tField, slpField]
    fields = [slpField]
    fields = [tField]
    nc_to_nps_int(inFile, outFile, date, xfcst, fields, source='g5nr', 
                  geos2wrf=False, log=log)
        
    # file-specific - Soil
    inFile = '/home/Javier.Delgado/scratch/nems/g5nr/apps/full_domain/OUTPUT/SURFACEMODEL/200609/LIS_HIST_200609100300.d01.nc'
    stField = ('SoilTemp_tavg', 'ST', 'K', 'soil temperature')
    fields = [stField]
    nc_to_nps_int(inFile, outFile, date, xfcst, fields, source='lis', 
                  geos2wrf=False, log=log)
    
    # Think: idemponence, since this will likely be very time consuming
    
    # TODO  ? 
    #
    # 0. geos2wrf feature has not been tested (just modifies the hdate - but
    #     does NPS read it?
    #
    # 1. Should the xfcst be the actual g5nr forecast hour or should it be 
    # relative to when we are initializing the regional forecast? With WRF we
    # usually use the GFS forecast that starts at the same time, so I don't
    # know if there is an assumption here.
    #
    # 2. Verify date against netCDF file.
    #  For LIS: begin_time (HHMMSS) and begin_date (ymd)
    #
    # 3. Missing value replacement does not seem to be working, at least for
    #     LIS (see np.where)
    #
    # 4. Performance, processing one 3-d variable (TT, 47 levs):
    #     python nc2nps.py  272.64s user 38.81s system 91% cpu 5:39.54 total
    #     ./unpack_netcdf.exe  118.90s user 9.51s system 83% cpu 2:33.32 total
    #    -> See profile_stats.txt. The vast majority of time is spend in
    #       get_2d_slab_from_var. May just be the fact that Numpy is compiled
    #       with GCC and probably not very optimized.
    
if __name__ == '__main__':
    main()
    #import profile
    #profile.run('main()')

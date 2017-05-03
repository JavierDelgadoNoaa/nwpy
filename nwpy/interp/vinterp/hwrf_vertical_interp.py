import os
import sys
import logging
from ConfigParser import ConfigParser
from distutils.util import strtobool

import numpy as np
from cfunits import Units # !! must be imported before matplotlib !!
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
from PyNIO import Nio
#import pdb ; pdb.set_trace()


from hwrf_real_interp import interp_press2press_lin
from pycane.logging.loggers import basic_logger
from nwpy.viz.map import bling as map_bling



#def interp_press2press_lin(plevs_in, plevs_out, input_var, dim_vars,
#                            extrapolate, ignore_lowest, t_field, log=None):
def interp_press2press_lin_wrapper(dataset, plevs_out, inputVarName, dimVars,
                                   extrapolate, outLevUnits, ignore_lowest=False, 
                                   t_field=False, dimDims=None, inLevUnits=None,
                                   log=None):
    """
    Extrapolate fields from a given dataset from an input array of pressure levs 
    to an output array of pressure levels. This will use the ``interp_press2press_lin''
    routine of HWRF, which does a log-P interpolation and does smooting if any values
    need to be extrapolated. The variable being interpolated must have the 
    dimensions (lev, lat, lon), in that order. It's okay to have a single time
    dimension as well. Supporting (lon, lat, lev) should be fairly simple, 
    since it's converted to that order anyway.
    Notes on efficiency
         - Since the HWRF routines require 3-d variables, even if doing isobaric 
           interpolation, a the output plevs will be converted to 3-d
         - The HWRF routine expects data to have shape i-j-k (lon-lat-pres). Most datasets
           have data in k-j-i (pres-lat-lon) shape. Since both Nio and nc4 store data
           in row-major order, we could conveivably just pass the data as-is 
           to the Fortran routines, but f2py insists on transposing, so we 
           transpose here as well.
           TODO : Look for way around this. e.g. hacking the F_CONTIGUOUS 
                  attribute or modifying intermediate code geneated by f2py
           Since vertical interpolation involves iterating through the Z/k dimension,
           k-i-j actually works best, so the ideal solution is to change
           the HWRF routine to use k-i-j
    Known Issues:
        - The original subroutines were modified as little as possible. While running
          in parallel should be possible, that would require modifying this
          function to properly set the T,D,M indices. Currently, they are 
          just hacked to run to function correctly in serial.
    :dataset The Nio or netCDF4 dataset to get data from
    :plevs_out The output pressure levels. May have different number of vertical
               dimensions than the input but the horizontal must be the same. 
               May be a 1-d array if interpolating to isobaric levels or a 3-d 
               array. If 3-d, it should be of shape (lev, lat, lon)
    :inputVarName 3-D input variable whose data is to be interpolated. 
    :dimVars 3-tupple consisting of the name of the latitude, longitude and 
             presure-level variables in ``dataset''; in that order. These 
             dimensions should correspond with those of ``inputVarName''.
    :extrapolate True to extrapolate pressure values out of the bounds of plevs3d
                 False to just use the top/bottom layer of source data
    :outLevUnits Units of output pressure levels
    :ignore_lowest Ignore the lowest vertical dimension in the input data when doing
                 the interpolation. e.g. if it has surface data.
    : t_field - True if ``inputVarName'' a temperature Field. False if not.
                If None, guess based on ``inputVarName'' and Variable attributes

    :dimDims - Name of the dimension variables. If not passed in, assume they are the
               same as the names of the dimension variables (dimVars)
    :inLevUnits Units of input (dataset's) pressure levels. Will be ignored if the
                dataset has a "units" attribute
    """
    if log is None: log = basic_logger()
    if dimDims is None: 
        dimDims = dimVars
    # Figure out lats and lons and levs
    (lat_var_name, lon_var_name, plev_var_name) = dimVars
    (lat_dim_name, lon_dim_name, plev_dim_name) = dimDims
    lat_var = dataset.variables[lat_var_name]
    lon_var = dataset.variables[lon_var_name]
    lons = lon_var[:]
    lats = lat_var[:]
    if lats.ndim == 1:
        log.debug("Converting lat and lon to 2-d")
        lats,lons = np.meshgrid(lats, lons)  # TODO ensure order is correct
        log.debug("Finished converting to 2-d")
    plevs_in = dataset.variables[plev_var_name][:]
    input_var = dataset.variables[inputVarName]
    # Squeeze out time dimension
    data_in = input_var[:]
    shapeBefore = data_in.shape
    data_in = data_in.squeeze()
    plevs_in = plevs_in.squeeze()

    # Ensure lat and lon are in the expected order
    latIdx = input_var.dimensions.index(lat_dim_name) 
    lonIdx = input_var.dimensions.index(lon_dim_name)
    levIdx = input_var.dimensions.index(plev_dim_name)
    # The arrays themselves will have the time dimension squeezed out, so they
    # need separate indices
    arLatIdx = latIdx ; arLonIdx = lonIdx ; arLevIdx = levIdx
    if shapeBefore != data_in.shape:
        assert input_var.dimensions[0].lower() == "time"
        arLatIdx -= 1
        arLonIdx -= 1
        arLevIdx -= 1
    if plevs_in.ndim == 1:
        assert data_in.shape[arLevIdx] == plevs_in.shape[0]
    else:
        assert data_in.shape == plevs_in.shape
    
    assert lonIdx > latIdx > levIdx # should apply for var and plevs
    dims_in = (plevs_in.shape[arLevIdx], data_in.shape[arLatIdx], data_in.shape[arLonIdx])
    dims_out = (plevs_out.shape[0], data_in.shape[arLatIdx], data_in.shape[arLonIdx]) # output will be isobaric slab
    log.debug("dims_in = {i}; dims_out = {o}".format(i=dims_in, o=dims_out))
    #import pdb ; pdb.set_trace()
    
    # Create 3-d input pressure array if it is 1-d
    if plevs_in.ndim == 1:
        plevs3d_in = np.empty(dims_in)
        plevs3d_in.T[:] = plevs_in # convert to 3d
    elif plevs_in.ndim == 3:
        plevs3d_in = plevs_in
    else:
        raise Exception("Unexpected shape for plevs_in")
    #levs = dataset.variables[plev_var_name][:] # TODO : This is already defined above as plevs_in !!
    if plevs3d_in[0,0,0] < plevs3d_in[1,0,0]:
        log.info("Will flip variables to go from bottom to top of atmosphere")
        flip_vertical = True
        plevs3d_in[:] = plevs3d_in[::-1,:,:]
    log.debug("input levels = {l}".format(l=plevs3d_in))

    #if flip_vertical:
    #    log.debug("Flipping plevs3d_in. New: {}".format(plevs3d_in))
    #    plevs3d_in[:] = plevs3d_in[::-1,:,:]
        # TODO ? flip plevs3d_out if going top to bottom

    #plev_out = 775 # just interpolating one slab
    #plev_out_3d = np.empty(dims_out)
    #plev_out_3d[:] = plev_out
    #if plev_out in plevs:
    #    print 'already there, no need to interpolate'
    #else:
    #    t = dataset.variables[field_mappings["temperature"]][:]
    if flip_vertical:
        log.info("FLipping variable across the vertical")
        #import pdb ; pdb.set_trace()
        assert data_in.shape[arLevIdx] == len(plevs_in)
        assert input_var.dimensions[levIdx] == plev_dim_name
        data_in[:] = data_in[::-1,:]
    data_out = np.empty(dims_out)

    # subR
    if plevs_out.ndim == 1:
        log.info("Requested plevs_out is 1-d (i.e. isobaric")
        plevs3d_out = np.empty(dims_out)
        plevs3d_out.T[:] = plevs_out # tile horizontally
    elif plevs_out.ndim == 3:
        plevs3d_out = plevs_out
    else:
        raise Exception("Unexpected dimensionality for plevs_out")
    # Set misc variables
    its=ids=ims = 1 # TODO (v) 0-n indexing
    jts=jds=jms = 1
    kts=kds=kms = 1
    ite=ide=ime = data_in.shape[arLonIdx]
    jte=jde=jme = data_in.shape[arLatIdx]
    kte=kde=kme = plevs3d_out.shape[arLevIdx] # +1 #  add one since it iterates kde-1 
    #kme -= 1 # ; jme -= 1 ; ime -= 1 
    # add one to domain dims since it iterates d-1
    ide += 1 ; jde += 1 ; kde += 1
    # Sanity checks
    log.debug("plevs3d_out shape: {0}".format(plevs3d_out.shape))
    assert plevs3d_in.flags['F_CONTIGUOUS'] is False
    assert plevs3d_out.flags['F_CONTIGUOUS'] is False
    assert data_in.flags['F_CONTIGUOUS'] is False
    
    # Process input data as necessary
    if isinstance(extrapolate, bool):
        extrapolate = 1 if extrapolate else 0
    else:
        raise ValueError("`extrapolate' should be bool") 
    # convert to Pa if necessary
    #if plevs3d_out[-1,0,0] < 1000:
    #    log.info("plevs3d_out appears to be in mb, converting to Pa")
    #    plevs3d_out *= 100.
    #if plevs3d_in[-1,0,0] < 1000:
    #    log.info("plevs3d_in appears to be in mb, converting to Pa")
    #    plevs3d_in *= 100.

    # The routine expects i-j-k, but data is k-j-i when loaded by Nio.
    # Since Nio does not use F_CONTIGUOUS memory order, it's actually
    # already in the correct order (I think), but f2py wants 
    # it to be the same and then do it's Transpose. This is wasteful
    # so we may want to see if we can work around this
    log.debug("plevs3d_in.shape = {0}".format(plevs3d_in.shape))
    plevs3d_in = plevs3d_in.T
    log.debug("plevs3d_in.shape (after Transpose) = {0}".format(plevs3d_in.shape))
    plevs3d_out = plevs3d_out.T
    data_in = data_in.T
    num_metgrid_levs = plevs3d_in.shape[2] # !! must be done after transpose
    #log.debug("generic = {0}".format(num_metgrid_levs))

    # guess if it's temperature, if necessary
    if t_field is None:
        try:
            varAttr = getattr(input_var, "long_name")
            if 'temperature' in varAttr.lower():
                t_field = True
                log.debug("Treating variable {var} as temperature for "
                          "interpolation, based on long_name attribute '{ln}'."
                          .format(var=inputVarName, ln=varAttr))
            else:
                t_field = False
        except AttributeError:
            if inputVarName.lower() in ("t", "temperature", "tt"):
                t_field = True
                log.debug("Treating variable {var} as temperature for "
                          "interpolation.".format(var=inputVarName))
            else:
                t_field = False
    
    # Conform to common units 
    #import pdb ; pdb.set_trace()
    try:
        inLevUnits = getattr(dataset.variables[plev_var_name], "units")
    except AttributeError:
        log.debug("Dataset does not have a 'units' attribute.")
        if inLevUnits is None:
            raise Exception("Unknown 'units'. Please "
                            "specify units as argument `inLevUnits'".format())
    if outLevUnits != inLevUnits:
        log.debug("Conforming output levels (in {uIn}) to units of input levels {uOut}"
                  .format(uIn=inLevUnits, uOut=outLevUnits))
        plevs3d_out = Units.conform(plevs3d_out, Units(outLevUnits), Units(inLevUnits))
    
    data_out = interp_press2press_lin(plevs3d_in, plevs3d_out, data_in, 
                    generic=num_metgrid_levs,
                    extrapolate=extrapolate, 
                    ignore_lowest=ignore_lowest,
                    tfield=t_field,
                    ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,
                    ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,
                    its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)

    '''
    data_out = interp_press2press_lin(plevs3d_in, plevs3d_out, data_in, 
                generic=num_metgrid_levs, 
                extrapolate=extrapolate, 
                ignore_lowest=ignore_lowest, tfield=t_field,
                ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde,
                ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,
                its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
    '''
    '''
    data_out = interp_press2press_lin(plevs3d_in, plevs3d_out, data_in, 
                extrapolate, 
                ignore_lowest, t_field,
                ids, ide, jds, jde, kds, kde,
                ims, ime, jms, jme, kms, kme,
                its, ite, jts, jte, kts, kte)
    
    '''
    # hwrf interp function returns i,j,k-ordered, so transpose to k,j,i
    assert data_out.flags['F_CONTIGUOUS']
    data_out = data_out.T
    #var2d = t_out.[0,:] 
    #log.debug("data_out: {0}".format(data_out))
    return data_out
    
#ata_out = interp_press2press_lin(press_in,press_out,data_in,extrapolate,ignore_lowest,tfield,ids,ide,jds,jde,kds,kde,ims,    ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte,[generic])

# interp_press2press_lin(press_in,press_out,data_in,data_out,generic,extrapolate,ignore_lowest,tfield,ids,ide,jds,jde,kds,k    de,ims,ime,jms,jme,kms,kme,its,ite,jts,jte,kts,kte)
    
if __name__ == "__main__":

    log = _default_log(log2stdout=logging.DEBUG)
    #infile = "/scratch4/NAGAPE/aoml-osse/Javier.Delgado/nems/gfs_data/HISTORY/GFS.2012/2012102118/gfs.t18z.pgrb2f60"
    infile = "/scratch4/NAGAPE/aoml-osse/Javier.Delgado/nems/gfs_data/HISTORY/GFS.2012/2012102212/gfs.t12z.pgrb2f00"
    field_mappings = {  "temperature":"TMP_P0_L100_GLL0",
                        "mslp":"PRMSL_P0_L101_GLL0",
                        "u3d":"UGRD_P0_L100_GLL0",
                        "v3d":"VGRD_P0_L100_GLL0",
                        "ght":"HGT_P0_L100_GLL0",
                        "av":"ABSV_P0_L100_GLL0"
                     }
    plev_var = "lv_ISBL0"
    plevs = np.array([10.,    20.,    30.,    50.,    70.,   100.,   
                     150.,   200.,  250.,  300.,   350.,   400.,   
                     450.,   500.,   550.,   600.,  650.,   700.,   
                     750.,   800.,   850.,   900.,   925.,   950.,  975.,  1000.]) 
    dataset = Nio.open_file(infile, format="grib2")

    # now draw stuff
    confbasic = lambda param: conf.get("BASIC", param)
    config_files = ["test.conf"]
    #log = _default_log(log2stdout=log_level)
    conf = ConfigParser()
    #conf.readfp(open(config_file))
    conf.read(config_files)
    basemap = _create_map(conf, dataset)
    
    var2d = t_out
    print "var2d.shape: ", var2d.shape, "var2d.T.shape: ", var2d.T.shape
    plot_func = plt.contourf
    filename = "out.png"
    plotsPerFigure  = 1
    plot_settings = dict(conf.items("t_filled"))
    map_settings = dict(conf.items("map_settings"))
    _plot_helper(lons, lats, var2d, basemap, plot_func, filename, plotsPerFigure,
                 contour_levels=None, colors=None, fontsize=10, ax=None,
                 plot_settings=plot_settings, fig=None, map_settings=map_settings,
                 basic_settings=None)
    plt.savefig("out.png")

"""
This module contains functions for creating WRF Preprocessing System (WPS)/
NEMS Preprocessing System (NPS) variables from variables produced by the
GEOS-5 modelling system.
A lot of the algorithms were taken from the geos2wrf utilities found in NUWRF.
All functions take in two arguments:
    varMaps should be a dictionary that maps variables to the path
    of the dataset they can be found in

"""
from PyNIO import Nio # nc4 should work too

def create_landsea(varMaps, log=None):
    """ 
    Create the LANDSEA field. See module notes for further info.
    """
    frLand_var = Nio.open_file(varMaps["FRLAND"]).variables["FRLAND"]
    frLandIce_var = Nio.open_file(varMaps["FRLANDICE"]).variables["FRLANDICE"]
    assert frLandIce_var.dimensions == frLand_var.dimensions
    dims = frLand_var.dimensions
    '''
    dims_lc = [d.lower() for d in dims]
    try:
        timeIndex = dims_lc.index("time")
        assert frLand_var[:].shape[timeIndex] == 0 <-attribute error
    except ValueError:
        timeIndex = None
    '''
    landsea_bin = (frLand_var[:] + frLandIce_var[:]).round(decimals=0)
    # NOTE: np.round() rounds to the nearest even .5. ITC, it rounds to 0
    units = "1" # proprtn
    long_name = "Land/Sea flag (1=land, 0 or 2=sea)"

    return (landsea_bin, dims, units, long_name)

def create_soilhgt(varMaps, log=None):
    """ Create the SOILHGT field. See module notes for further info """
    phis = Nio.open_file(varMaps["PHIS"]).variables["PHIS"]
    data = phis[:] / 9.80665 # divy by gravity@msl
    units = "m" # GPM, but for cfunits, just use m
    long_name = "Soil Height"
    dims = phis.dimensions
    return (data, dims, units, long_name)

def create_rh(varMaps, log=None):
    """ 
    Create the RH (wrt liquid) field . See module notes for further info.
    Based on geos2wrf "WrfRH" module, which in turn is based on WRF.
    NOTE: This function only does the 3-d RH. To do the surface RH, you also 
    need PSFC in varMaps
    """
    
    # Constants from Teten's formula
    SVP1 = 0.6112
    SVP2 = 17.67
    SVP3 = 29.65
    SVPT0 = 273.15
    # Molecular weights of water vapor and dry air
    MW_VAP = 18.0152
    MW_AIR = 28.966

    calcMR = lambda sh : sh / (1.0 - sh)

    g5nr_qv = Nio.open_file(varMaps["QV"])
    g5nr_t = Nio.open_file(varMaps["T"])
    g5nr_pl = Nio.open_file(varMaps["PL"])
    
    mr = calcMR(g5nr_qv.variables["QV"][0,:])
    t_var = g5nr_t.variables["T"]
    pmid_var = g5nr_pl.variables["PL"] # midlayer pressure
    qv_var = g5nr_rh.variables["RH"]    
    assert getattr(pmid_var, "units") == "Pa"
    assert getattr(t_var, "units") == "K"
    assert t_var.dimensions == pmid_var.dimensions == qv_var.dimensions
    
    # Saturation vapor pressure and spechumd
    satVaporPressurePa = 100. * SVP1 * 10. * np.exp( SVP2*(t_var[0,:] - SVPT0)/(t_var[0,:]-SVP3))
    satSpecificHumidity = MW_VAP * satVaporPressurePa
    satSpecificHumidity = satSpecificHumidity /(satSpecificHumidity + MW_AIR * (pmid_var[0,:] - satVaporPressurePa))

    sat_mr = calcMR(satSpecificHumidity)
    rh_pct = mr / sat_mr * 100.
    #rh_pct2 = np.copy(rh_pct )
    rh_pct2 = rh_pct
    rh_pct2[rh_pct2 < 0.] = 0.
    rh_pct2[rh_pct2 > 100.] = 100.
    
    return (data_rh_g5nr, t_var.dimensions, "%", "relative_humidity")

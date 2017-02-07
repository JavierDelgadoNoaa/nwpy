import os
import subprocess
import logging

from objects import SpecifiedForecastDataset

"""
This module provides utilities to encapsulate information about different input
sources. 
The guess_input_source attempts to guess the input source, but note that it is
based on specific knowledge of a few different input types.

To add a new input source:
 1. (probably not necessary:) Add file type to _GENERAL_FILE_TYPES and 
    _SPECIFIC_FILE_TYPES if not already there
 2. Add source to _GENERAL_INPUT_SOURCES 
 3. Add source to _SPECIFIC_INPUT_SOURCES 
 4. Add condition in guess_input_source() to (uniquely!) identify the source
 5. (If this is a new source in _GENERAL_INPUT_SOURCES):
     i. Create a subdirectory under the "inputspec" directory with the name of the
        general input source
    ii. Create an input specification ("inputspec") file describing the contents
        the general input source. The file name should match the given 
        general input source name
 6. Create an input specification ("inputspec") file describing the contents
    of the specific input source dataset. The file name should match the
    given specific input source name.
"""

_GENERAL_FILE_TYPES = ["netcdf", "grib", "nemsio"]
_SPECIFIC_FILE_TYPES = ["nc4", "nc3", "grib1", "grib2", "nemsio"]
_GENERAL_INPUT_SOURCES = \
    ["Unknown", "wrf", "nems", "upp", "gfs", "geos5"]
_SPECIFIC_INPUT_SOURCES = \
    ["g5nr_hires", "g5nr_coarse_isobaric",  "g5nr_hires_nps_isobaric", 
     "g5nr_hires_nps", "hwrf_grib2", "upp_hwrf_grib1", "upp_hwrf_grib2", 
     "upp_nmb_grib1", "upp_nmb_grib2", "nmb_nemsio", "nmb_nemsio2netcdf"
     "gfs_generic", "gfs_prod2012", "gfs_prod2014", "hwrfprs_grib2",
     "hwrfprs_grib1", "hwrfsat_grib2", "hwrf_modelLev_netcdf"]

_CANNED_INPUTSPECS = os.path.join(os.getcwd(), "conf", "inspec")

class InputSource(object):
    """
    This class encapsulates the attributes of an input source.  It ensures
    that the attributes set by callers are valid (i.e. that they exist 
    in the values set in the global vairables _GENERAL_INPUT_SOURCES, 
    _GENERAL_INPUT_SOURCES, _SPECIFIC_INPUT_SOURCES, _SPECIFIC_FILE_TYPES)
    """
    def __init__(self, path):
        self._path = path
        self._general_file_type = None
        #self._specific_file_type = self.guess_specific_file_type # callback
        self._specific_file_type = None
        self._specific_input_source = None
        self._general_input_source = None

    @property 
    def general_file_type(self):
        return self._general_file_type
    @general_file_type.setter
    def general_file_type(self, value):
        self._general_file_type = self.__set_if_valid(value, 
                                             _GENERAL_FILE_TYPES)
    @property
    def specific_file_type(self):
        #import pdb ; pdb.set_trace()
        if self._specific_file_type is None and self._general_file_type is not None:
            return self.guess_specific_file_type()
        return self._specific_file_type
    @specific_file_type.setter
    def specific_file_type(self, value):
        self._specific_file_type = self.__set_if_valid(value,
                                              _SPECIFIC_FILE_TYPES)
    @property 
    def general_input_source(self):
        return self._general_input_source
    @general_input_source.setter
    def general_input_source(self, value):
        self._general_input_source = self.__set_if_valid(value,
                                                _GENERAL_INPUT_SOURCES)
    @property 
    def specific_input_source(self):
        return self._specific_input_source
    @specific_input_source.setter
    def specific_input_source(self, value):
        self._specific_input_source = self.__set_if_valid(value, 
                                                _SPECIFIC_INPUT_SOURCES)

    def __set_if_valid(self, value, possible_values):
        if not value in possible_values:
            raise Exception("Unknown input source '{0}'. Known values are '{1}'"
                            .format(value, possible_values))
        return value
        
    def guess_specific_file_type(self):
        assert self.general_file_type == "grib" # no logic for netcdf
        filename = os.path.basename(self._path)
        with open(os.devnull, 'w') as devnull:
            # check if wgrib2 exists
            #try:
            #    subprocess.check_call(["wgrib2", "?"])
            #except OSError:
            #    raise Exception("wgrib2 not in path?")
            #import pdb ; pdb.set_trace()
            # TODO : Put full path to wgrib2 and/or specify in some config
            cmd_args = ["wgrib2", "-s", self._path]
            #subprocess.check_call(cmd_args)
            out = subprocess.Popen(cmd_args, stderr=subprocess.PIPE, 
                                   stdout=devnull).communicate()[1]
            # wgrib does not fail if it's grib1, it just ignores the message
            if "grib1 message ignored" not in out:
                return "grib2"
            # TODO : Need full path to wgrib : /apps/wgrib/1.8.1.0b/bin/wgrib
            cmd_args = ["wgrib", "-s", self._path]
            out = subprocess.Popen(cmd_args, stderr=subprocess.PIPE,
                                   stdout=devnull).communicate()[1]
            if "grib2 message ignored" not in out:
                return "grib1"
        raise Exception("Unable to determine specific file type")

def guess_input_source(path, log=None):
    """
    Attempts to guess the input source of the given file based on its file name.
    Note that this is based on pretty specific knowledge about file naming
    conventions, so it's not a comprehensive thing.
    :return insrc: An `InputSource' object. HINT: If you need an Inputspec 
    list, use filepath2inputspec instead
    """
    if log is None: log = _default_log()
    filename = os.path.basename(path)
    insrc = InputSource(path)

    if filename.startswith("nmmb") and filename.endswith(".nc"):
        log.debug("Guessing this is an NMM-B run due to filename prefix")
        insrc.general_file_type = "netcdf"
        insrc.general_input_source = "nmmb"
        insrc.specific_input_source = "nmmb_nemsio2netcdf"
        insrc.specific_file_type = "nc3"
    elif filename.startswith("c1440"):        
        log.debug("Guessing this is a G5NR dataset due to file prefix")
        insrc.general_input_source = "geos5"
        insrc.general_file_type = "netcdf"
        insrc.specific_file_type = "nc4"
        #hires: c1440_NR.inst30mn_3d_T_Nv.20060909_0000z.nc4
        #coarse: c1440_NR.inst01hr_3d_H_Cp.20060111_1800z.nc4
        if "Nv" in filename:
            log.debug("Further guessing it is high-res due to 'Nv' in file")
            insrc.specific_input_source = "g5nr_hires"
        elif "Cp" in filename:
            log.debug("Further guessing it is low-res due to 'Cp' in file")
            insrc.specific_input_source = "g5nr_coarse_isobaric"
        elif "Cv" in filename:
            log.debug("Further guessing it is low-res due to 'Cv' in file")
            insrc.specific_input_source = "g5nr_coarse_native"
        else:
            raise Exception("Unable to determine specific g5nr input source")
    elif filename.startswith("g5nr_combined_hires"):
        log.debug("Guessing this is a `g5nr_nps' file generated with nr_input_generator")
        insrc.general_file_type = "netcdf"
        insrc.specific_file_type = "nc4"
        insrc.general_input_source = "geos5"
        insrc.specific_input_source = "g5nr_hires_nps"
    elif filename.startswith("g5nr_nps_isobaric"):
        log.debug("Guessing this is a full resolution G5NR dataset interpolated"
                  " to isobaric levels based on file name".format())
        insrc.general_file_type,insrc.specific_file_type = ("netcdf", "nc4")
        insrc.general_input_source = "geos5"
        insrc.specific_input_source = "g5nr_hires_nps_isobaric"
    elif filename.startswith("g5nr_coarse_merged"):
        log.debug("Guessing this is a coarse resolution G5NR dataset")
        insrc.general_file_type,insrc.specific_file_type = ("netcdf", "nc4")
        insrc.general_input_source = "geos5"
        insrc.specific_input_source = "g5nr_coarse_merged"
    elif filename in ("input_domain_01_nemsio.nc", "input_domain_01.nc"):
        log.debug("Guessing this is a NemsInterp file converted to netCDF."
                  " Will treat like a standard history file.".format())
        insrc.general_file_type,insrc.specific_file_type = ("netcdf", "nc3")
        insrc.general_input_source = "nmmb"
        insrc.specific_input_source = "nmb_nemsio2netcdf"
    elif filename.startswith("gfs."):
        log.debug("Guessing this is a GFS dataset based on filename")
        insrc.general_file_type,insrc.specific_file_type = ("grib", "grib2")
        insrc.general_input_source = "gfs"
        if filename.index("pgrb2f") > 0:
            log.debug("Further guessing it is a PROD-2012 GFS file")
            insrc.specific_input_source = "gfs_prod2012"
        elif filename.index("0p25") > 0:
            log.debug("Further guessing it is a PROD-2014 GFS file")
            insrc.specific_input_source = "gfs_prod2014"
        else:
            raise Exception("Unable to determine kind of GFS file")
    elif filename.startswith("nmbprs") or filename.startswith("wrfprs"):
        log.debug("Guessing this is a UPP dataset based on filename")
        insrc.general_file_type = "grib"
        insrc.general_input_source = "upp"
        # call to specific_file_type decorator will guess if grib1/2
        sft = insrc.specific_file_type
        if filename.startswith("nmb"):
            if sft == "grib1":
                insrc.specific_input_source = "upp_nmb_grib1"
            elif sft == "grib2":
                insrc.specific_input_source = "upp_nmb_grib2"
        elif filename.startswith("wrf"):
            if sft == "grib1":
                insrc.specific_input_source = "upp_wrf_grib1"
            elif sft == "grib2":
                insrc.specific_input_source = "upp_wrf_grib2"
    elif filename.startswith("wrfout"):
        log.debug("Guessing this is a WRF model output")
        insrc.general_file_type = "netcdf"
        # TODO : Need logic in guess_specific_file_type() to figure out if it's nc3/4
        insrc.general_input_source = "WRF"
        insrc.specific_input_source = "hwrf_modelLev_netcdf"
        # TODO : Need additional logic to distinguish between HWRF/NMM/ARW. e.g. 
        #        look in netCDF file
    #invest09l.2012082012.hwrfanl_i.grb2f00
    #invest09l.2012082012.hwrfprs.d1.0p25.f000.grb2
    #invest09l.2012082012.hwrfges_i.grb2f00
    #invest09l.2012082012.hwrfprs.d123.0p10.f000.grb2
    #invest09l.2012082012.hwrfprs.d2.0p10.f000.grb2
    #invest09l.2012082012.hwrfsat.d123.0p25.f000.grb2
    #invest09l.2012082012.hwrftrk.f000.grb
    elif "hwrfprs" or "hwrftrk" or "hwrfsat" or "hwrfges" in filename:
        log.debug("Guessing this file was produced by the operational HWRF")
        insrc.general_file_type = "grib"
        insrc.general_input_source = "WRF"
        if "grb" in filename:
            insrc.specific_file_type = "grib1"
        elif "grb2" in filename:
            insrc.specific_file_type = "grib2"
        log.debug("Further guessing it's {0}".format(insrc.specific_file_type))
        if "hwrfprs" in filename or "hwrfges" in filename \
           or "hwrfanl" in filename or "hwrftrk" in filename:
            insrc.specific_input_source = "hwrfprs_" + insrc.specific_file_type
        elif "hwrfsat" in filename:
            insrc.specific_input_source = "hwrfsat_" + insrc.specific_file_type
        else:
            raise Exception("Unable to determine specific HWRF source")
    else:
        raise Exception("Unable to determine the type of file used")
    
    return insrc

def filepath2inputspec(path, inputspecs_topdir=None):
    """
    Attempt to map a file, given it's `path', to one of the canned input 
    specification files included in the nwpy distribution.
    :return: -- A list of inputspecs, in order of specific-ness. e.g. 
                the field names in grib are universal (i.e. the same whether
                it was generated by UPP, GFS, etc) since they are specified
                by Nio when decoding the grib messages. Nio uses different
                naming conventions for grib1 vs. grib2.
    """
    ret = []
    if inputspecs_topdir is None:
        inputspecs_topdir = _CANNED_INPUTSPECS
    insrc = guess_input_source(path)
    general = os.path.join(inputspecs_topdir, insrc.general_input_source,
                           insrc.general_input_source + ".conf")
    specific = os.path.join(inputspecs_topdir, insrc.general_input_source, 
                            insrc.specific_input_source + ".conf")
    if insrc.general_input_source is "grib":
        grib = os.path.join(inputspecs_topdir,
                            insrc.general_file_type + ".conf")
        ret = [grib, general, specific]
    ret = [general, specific]
    # TODO? What if we want multiple inputspecs per file?
    for path in ret:
        if not os.path.exists(path):
            raise Exception("{0}: Path to inspec does not exist".format(path))
    return ret

def specdata_from_path(path, date, inputspecsTopdir=None, fcstOffset=None, 
                      domain=None, log=None):
    """
    Create a SpecifiedForecastDataset from a given file path.
    :param path: Full path to the file. filepath2inputspec will be used to 
                 guess the input source and get corresponding inputspec
    :type path: str
    :param date: Initialization date of forecast corresponding to the file
    :type date: datetime.datetime
    :param inputspecsTopdir: Directory containing input specification files. 
                             If None, use the specdata module canned ones
    :type inputspecsTopdir: str
    :param fcstOffset: Offset from the `date', if applicable
    :type fcstOffset: datetime.tdelta
    :param domain: The domain number, if applicable
    :type domain: int
    """
    inputspec = filepath2inputspec(path, inputspecsTopdir)
    dir_name = os.path.dirname(path)
    #if fcstOffset is not None:
    return SpecifiedForecastDataset(inputspec, dir_name, date,
                                    fcst_offset=fcstOffset, domain=domain, 
                                    log=log)

def _default_log(log2stdout=logging.INFO, log2file=None, name='el_default'):

    # Need this or when calling _default_log from multiple classes or else
    # if I don't pass it to constructor, I end up with multiple loggers
    global _logger
    if _logger is not None:
        return _logger

    _logger = logging.getLogger(name)
    _logger.setLevel(log2stdout)
    msg_str = '%(asctime)s::%(name)s::%(lineno)s::%(levelname)s - %(message)s'
    msg_str = '%(asctime)s::%(funcName)s::%(filename)s:%(lineno)s::%(levelname)s - %(message)s'
    date_format = "%H:%M:%S"
    formatter = logging.Formatter(msg_str, datefmt=date_format)
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
_logger = None

if __name__ == "__main__":
    # Test if guess_specific_file_type works
    #src = InputSource("/home/Javier.Delgado/scratch/nems/g5nr/data/new_beta_pl/800x800/postprd/nmbprs_d01.00")
    src = InputSource("/home/Javier.Delgado/scratch/prod2015_sample_output/GDAS1/2015/2015092706/gdas1.t06z.pgrb2.0p25.f001")
    #import pdb ; pdb.set_trace()
    src.general_file_type = "grib"
    print src.specific_file_type

    # Basic test
    src = InputSource(None)
    src.specific_file_type = "nc4"
    src.specific_file_type = "foo"
    src.foo = "bar"
    src.specific_input_source = "g5nr_coarse_isobaric"
    src.specific_input_source = "foo"

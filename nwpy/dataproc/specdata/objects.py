"""
This module defines the objects used for working with "specified" forecast data, 
including:
SpecifiedForecastDataset
ConformedField
HorizontalCoordinate
SpecificVerticalLevel

To maintain synchronicity among transforms and to keep the code as intuitive as
possible, there is no interaction between any two classes except SpecifiedForecastDataset;
e.g. if a ConformedField wants to access it's corresponding HorizontalCoordinate, it 
     must do so through it's SpecifiedForecastDataset ( ``spec_data'' ) parameter

"""

"""
NOTES:

 - open_file returns the C object NioFileObject

LIMITATIONS:
 - The NIO grib library does not read fields past 2GB in size

*************************************************************************
--l100--
XX Grib2: level_type :   Isobaric surface (Pa)
XX Grib1: level_indicator :      100
 -> In Grib1, level 100 is hPa    
    You can confirm with getattr(dimvarname, "units") - but for 2d
    this won't be possible, so use assumed value
parameter_template_discipline_category_number :        [0, 0, 1, 1]
dim name == var name

--2L106--
XX Grib2: level_type :   Depth below land surface (m)
XX Grib1: level_indicator :      112
 -> this dim has two corresponding variables (for range min and max).
    e.g. dim=lv_DBLL3 ; vars = lv_DBLL3_l0 and lv_DBLL3_l1
          print gfs.variables["lv_DBLL3_l0"][:] = [ 0.          0.1         0.40000001  1.        ]
          gfs.variables["lv_DBLL3_l1"][:] = [ 0.1         0.40000001  1.          2.        ]
          getattr(gfs.variables["lv_DBLL3_l1"], "units") = m
** For UPP Grib1, paradigm is the same, but units are in cm
parameter_template_discipline_category_number :        [0, 2, 0, 192]


--2L108--

This may be a 2-d var with a "level" attribute specifying the delta
    e.g. level = [3000, 0] -> i guess that means it's 30mb above ground
level_type: Level at specified pressure difference from ground to level (Pa)
parameter_template_discipline_category_number :        [0, 0, 1, 1]
level :        [3000, 0]
??

--l102--
May be 2d:
XX Grib2: level_type :   Specified height level above ground (m)
XX Grib1: level_indicator :      105
level: 2 (** same for grib1 and grib2, ITC) 
 -> rh @ 2m above ground
OR 3d:
   level_type: Specified height level above ground (m)
   parameter_template_discipline_category_number :        [0, 0, 2, 3] or [0, 0, 1, 1]
   Dimensions and sizes:   [lv_HTGL9 | 3] x [lat_0 | 361] x [lon_0 | 720]
      -> height dim name = var name
*** Presumably, all vars can be either 2d or 3d

** The units are implied by the parameter_number, although they
   are specified in the Nio  Variable attribute.

--l103--
level_type :   Specific altitude above mean sea level (m)
   parameter_template_discipline_category_number :        [0, 0, 2, 3]
dim name == var name (e.g. lv_AMSL1 )
----

level_type :   Isobaric surface (Pa)
 -> read dimension directly
*************************************************************************
"""

"""
TODO: 
1. Currently, we create a new SpecifiedForecastDataset for every file. Is it
worth caching some time-invariant things like field2file, coord. mappings, etc.?
"""
import os
import sys
import logging
from ConfigParser import ConfigParser,NoOptionError
from distutils.util import strtobool
import glob
from datetime import datetime as dtime
from datetime import timedelta as tdelta
import numpy as np
import copy
import collections
import inspect

from PyNIO import Nio
import ESMF # Just used for constants for now, so a proxy module 
            # can be used on systems that do not have ESMF
from cfunits import Units
# Equation parser: 'https://github.com/alphaomega-technology/Equation
#  -> added line to support log10 function
# TODO  put this comment in some dependency file
from  Equation import Expression


class SpecifiedForecastDataset(object):
    """
    This class facilitates accessing and working with data from different input
    sources with a "specified" configuration (i.e. a dataset described using 
    config file(s) known as "input specifications."
    It attempts to solve the following problems:
        1. Since different input sources have different names for different 
           model fields, it uses "standard" field names, which are mapped to native
           names via the config file (i.e. "input specification" or "inputspec")
        2. It facilitates the use of fields in different coordinate systems (i.e.
           level types, staggering, etc.) by
           encapsulating the process of retrieving the lat/lon coordinates of a 
           given field, according to mappings set in the inputspec.
        3. It also facilitates access to fields with different vertical level types
           by mapping fields to their vertical levels and specifying which fields
           in the dataset contain the level values
        4. Processes additional basic metadata from the inputspec that may not 
           be available in the data, such as fields' units
        5. Performing transforms on fields, such as shifting grid horizontally so 
           it can be properly processed with Basemap
        6. Automatically converting units
        7. This module works on top of PyNIO. Although PyNIO does a great job of 
           providing a familiar (netCDF-like) 
           interface for all datasets, files of different formats must be treated
           slightly different since they fields have different ways of describing 
           themselves. e.g. netCDF files generally require more information in the 
           inputspec to specify what vertical and horizontal variables correspond 
           to the data for a given field's dimensions.
           GriB files are self-described, but their naming convension is more 
           complex, so standard to native field mappings require more complex
           specification.
           This is something else that this object encapsulates
        8. In some datasets, separate 2-D variables are used for 3-D fields 
           at different levels of the same level type. e.g. in WRF/NMM-B, there
           are separate U10 and U30 fields for U at 10 and 30 meters above 
           the surface. NIO would treat this as a 3-D field. This module allows
           separate variables to be accessed as if they were a single 3-D 
           variable by specifying as such in the InputSpec

    For detailed information about specific 
    config parameters, see the sample config file (inputspec.conf). 
    Summary of sections:
    
    [BASIC]
    Specifies file name pattern and the input source

    [field2filename_mappings]
    Specifies what file different fields can be found in. If all fields are in a 
    single file, having this section is not necessary.
    TODO: If you do not want to specify individual field-to-file mappings, you 
    can just set the file names and the code will look for fields in each file.

    [field_mappings]
    Map native field names to standard names.

    [units]
    Needed for datasets that do not have units in the field names. Ignored otherwise.
    e.g. nC files generated with nemsio2netcdf

    [field2lev_mappings]
    Map fields to their vertical level dimension and dimension variable names
    [latlon2field_mappings]
    Maps lat,lon variables to fields that use them. Dimension mappings are not 
    needed since the dimensions are an attribute of each Variable. However, the
    variable names of the variables containing the dimension data may be different.
    e.g. HWRF has the "west_east" dimension and corresponding variables "HLON/VLON"
    TODO: If a mapping is not specified here, attempt to use a variable with the same 
    name as the dimension

    [horizontal_configuration]
    Specify dimension names for the horizontal coordinate configurations corresponding to the 
    TODO : Change section name to "horizontal_coordinate_dims"

    [vertical_configuration]
    
    [vertical_level_mappings]
    [var2dim_mappings]
    [time_dim_config]

    
    """
    
    
    def __init__(self, inputspec_paths, topdir, init_date, fcst_offset=None, 
                 domain=1, file_format=None, fieldTransforms=None, 
                 coordTransforms=None, inspecsTopdir='', log=None, 
                 logLevel=logging.INFO):
        """
        :param inputspec: Config file that specifies contents of the data files
        :type inputspec: str
        :param topdir: Root directory wherein files are found
        :type topdir: str
        :param init_date: Date the forecast was initialized
        :type init_date: datetime.datetime
        :param fcst_offset: The offset from the init_date (e.g. the forecast hour)
        :type fcst_offset: datetime.timedelta
        :param domain: The domain/grid number
        :type domain: int
        :param file_format: The file format ("netcdf", "grib", etc.) Only needed if 
                            files do not have a common extension
        :type file_format: str
        :param log: Logger object to use. 
        :type log: logging.Logger
        :param logLevel: Log level (only used if log is not passed in)
        :type logLevel: int
        """
        
        # * TODO : Need to verify init_date in dataset. Otherwise, for outputs
        #         that do not have init date in file name (e.g. UPP default and nmmb)
        #         we may be getting the wrong thing
        # ** TODO : Read thru this and ensure it makes sense for both netCdf and GriB
        self._log = log if log is not None else _default_log(logLevel)
        self.variables = []
        self.topdir = topdir
        self.init_date = init_date
        if fcst_offset is None:
            fcst_offset = tdelta(seconds=0)
        self.fcst_offset = fcst_offset
        self.domain = int(domain)
        self.file_format = file_format
        
       ##
        # Set constants
        ##
        self._date_conversion_args = {
                            "cal_date":self.init_date + self.fcst_offset, 
                            "fdate":self.init_date + self.fcst_offset,
                            "init_date":self.init_date, "dom":self.domain,
                            "adate":self.init_date,
                            "fhr": int(self.fcst_offset.total_seconds() / 3600.), 
                            "fmin": int(self.fcst_offset.total_seconds() % 3600. / 60.),
                            "fcst_offset_mins": int(self.fcst_offset.total_seconds() / 60.),
                            "dom":self.domain
                                     }
        # TODO : Support other level types - will require corresponding functionality
        # in multiple methods/functions. Search for "level_types" and "isobaric"
        # Will also require corresponding inputspec additions
        self.level_types = ["isobaric", "sfcDelta", "mslDelta", "soilDelta"]
                         #"inter_sfc_depth", "inter_iso_sfc", "mslDelta",]
        ##
        # Set state vars
        ##
        
        """
        Map native field names to the files they can be found in. Will include
        mappings for all fields in the [BASIC]file_name and for all mappings 
        specified in [field2filename_mappings].s
        """
        self._field2file = {}
        
        """
        Maintain a dictionary of standard to native field name mappings. This 
        is actually a dictionary of dictionaries. The outer dictionary maps the
        level type ("isobaric", "2d", "sfcDelta", ...) to the corresponding fields' standard-to-native mappings.
        An entry will be created for each element in self.level_types
        """
        self._standard2native = {}
        
        """
        The self._var2dim variable is a dictionary with 4 keys:
        "time", "lats", "lons", "levs"
        Each key maps to another dictionary that maps variable names to 
        dimension names for the corresponding coordinate.
        e.g. var2dim["lats"]["HLAT"] = south_north
        """
        self._var2dim = { "time":{}, "lats":{}, "lons":{}, "levs":{} }

        """
        Maintain a dictionary of SpecifiedVerticalLevel's, since multiple 
        fields can map to the same level.
        The key for each entry is a 4-tuple:
           (lev_type, stagger, dimname, dimvarname)
        """
        self._specified_levels = {}
       
        """
        Map field names to their vertical dimension and the variable containing
        the dimension data. i.e.: nativ_field_name->(dimName, dimVarName)
        """
        # TODO ? Does using the native name work for Grib?
        self.field2lev = {}

        """
        Dynamically populated mapping of field name to SpecifiedVerticalLevel 
        object
        """
        self._field2speclev = {}

        """
        Keep a dict mapping lat+lon var names to HorizontalCoordinate objects, 
        since only one HorizontalCoordinate is needed for all fields that 
        use the same coordinates. The key will be latVarName+lonVarName.
       """
        self.__horizontal_coords = {}
        
        """
        Map fields to HorizontalCoordinate objects in self.__horizontal_coords. 
        This dictionary will have the mappings that appear in the inputspec
        section [latlon2field_mappings]. Any fields that do not have a mapping 
        there should be self-described (i.e. latVarName == latDimName) and the 
        lat dimension should have "lat" in it's name and the longitude dimension 
        should have "lon" in its name.
        """
        self._field2coord = {}
        
        """
        Map fields to their ESMF stagger location
        """
        self._field2stagger = None
        
        """
        List of derived fields (i.e. fields generated from other native fields
        """
        self._derived_fields = {}
       
        """
        ESMF Grid object corresponding to this dataset
        ASSUME : Grid does not change with time
        The Grid should be obtained via get_esmf_grid_for_field() rather than
        using this object directly
        """
        self._esmf_grid = None

        """
        Dictionary containing "Multivariable" fields, which are 3-d fields made 
        up of multiple 2-d Variables.
        Each entry will be a dictionary with the following keys:
          -varNames - List of (2d) Variable names
          -levVals - List of level values
          -units - The units corresponding to the level variables
          -dimOrder - Dimension order (string of 'i','j','t')
        """
        self._multivar_fields = {}
        for levType in self.level_types: self._multivar_fields[levType] = {}

        if inspecsTopdir == '':
            self._log.warn("inspecsTopdir not passed in. Setting to current dir")
            inspecsTopdir = os.get_cwd()
        self.inputspec, self.inputspec_paths = self._set_inspec(inputspec_paths,
                                                                inspecsTopdir)

        if self.file_format is None:
            try:
                self.file_format = self.inputspec.get("BASIC", "file_format")
            except NoOptionError:
                self._log.warn("file_format not specified in inspec. May not be "
                               "able to read it")

        #
        # Populate state vars from inputspec
        #
        #import pdb ; pdb.set_trace()
        self._log.debug("Setting field2file")
        std2file_mappings = self.__set_field2file_mappings() 
        self._log.debug("Setting standard2native")
        self.__set_standard2native_mappings() 
        # Now that we have standard2native, we can set filename mappings from standard
        # name (ie std(...) entries in field2filename_mappings
        # TODO - best to split up self.__set_field2file_mappings and not have
        # so much code here
        if std2file_mappings and len(std2file_mappings) > 0:
            interp_file = lambda f: f.format(**self._date_conversion_args)
            self._log.debug("Setting std2filename")
            for (levType,fieldName,path) in std2file_mappings:
                    if not fieldName in self._standard2native[levType]:
                        raise Exception("Field with standard name '{0}' and "
                                        "level type '{1}' is not in "
                                        "standard2native mappings"
                                        .format(fieldName, levType))
                    var = self._standard2native[levType][fieldName]
                    self._field2file[var] = interp_file(path)

        self._log.debug("Setting var2dim")
        self.__set_var2dim_mappings()
        self._log.debug("Setting field2lev")
        self.__set_field2lev_mappings()
        #self.__set_derived_fields() # TODO - next version
        
        # Set transforms that will be set by default for each field
        self._field_transforms = {} 
        if fieldTransforms is not None:
            self._field_transforms = fieldTransforms
        # Set trasnforms that will be set by default for each HorizontalCoordinate 
        self._coord_transforms = {}
        if coordTransforms is not None:
            self._coord_transforms = coordTransforms 

        self._dataset_resolution = None

    def _set_inspec(self, inputspec_paths, topdir):
        """ 
        Build inputspec given a list of paths to read from. This method will
        recursively read the inspec pointed to by the parameter `parent_inspec'
        in the [BASIC] section of the given config files, up until no more
        remain. It will then create the final ConfigParser using the
        sorted list of files.
        :param inputspec_paths: List of inputspec config file paths. Path
               may be absolute (starts with "/") or relative to topdir
        :param topdir: Root directory of all inputspecs
        :return 2-tuple : ConfigParser containing inputspec and list of 
                          all files used.
                                
        """
        # All given configs and all their children, their children, etc
        combined_config_files = [] 
        #import pdb ; pdb.set_trace()
        for i,path in enumerate(inputspec_paths):
            if not path.startswith("/"):
                path = os.path.join(topdir, path)
                inputspec_paths[i] = path
            if not os.path.exists(path):
                raise IOError(path)
        for child_path in inputspec_paths:
            #child_path = os.path.join(topdir, child_path)
            conf = ConfigParser()
            config_files = [child_path]
            have_more_children = True
            #import pdb ; pdb.set_trace()
            while have_more_children:
                if not os.path.exists(child_path):
                    raise IOError("Inputspec path does not exist: {0}"
                                  .format(child_path))
                conf.read(child_path)
                #conf.read(config_files) # need less-specific ones for interp
                try:
                    child = conf.get("BASIC", "parent_inspec").strip()
                except NoOptionError:
                    have_more_children = False
                    break
                #topdir = conf.get("BASIC", "inspec_topdir")
                if child.startswith("/"):
                    child_path = child
                else:
                    child_path = os.path.join(topdir, child)
                config_files.append(child_path)
                #config_files.reverse()
                conf = ConfigParser()
            config_files.reverse()
            self._log.info("Will read the following config files for input "
                           "specification: {0}"
                           .format(config_files))
            combined_config_files.extend(config_files)
        conf = ConfigParser()
        conf.optionxform = str # don't convert option names to lowercase
        conf.read(combined_config_files)
        #import pdb ; pdb.set_trace()
        return conf,combined_config_files

    def get_filename(self, fieldName, levType=None):
        """ 
        Return name of file containing field `fieldName' of level type 'levType' 
        """
        msg_unk = "Unknown Field: '{fld}' with levType '{levType}'. If '{fld}'"\
                  "is the standard_name, ensure mapping exists in inputspec "\
                  "({inspec}) [*_field_mappings]. Also, if applicable, ensure "\
                  "there is a mapping in [field2filename_mappings]"
        
        file_name = None
        if fieldName in self._field2file:
            file_name = self._field2file[fieldName]
            native_fld_name = fieldName
            # TODO ? set standard field name, which is passed to
            #     ConformedField constructor
            #if levType is not None:
                #try: 
                    #std_fieldName = self._standard2native
            self._log.info("Since native field name given, will pass None as "
                "standard_name to ConformedField constructor")
            std_fieldName = None
        else:
            # Maybe standard_name given
            if levType is None:
                if fieldName in self._standard2native["2d"]:
                    levType = "2d"
                else:
                    raise Exception("If passing `standard_name' to `get_field()'"
                                    ", you must pass in the levType (unless "
                                    "field is 2D). If it's 2D, ensure there's a"
                                    "mapping in inputspec [2d_field_mappings]")
            try:
                native_fld_name = self._standard2native[levType][fieldName]
            except:            import pdb ; pdb.set_trace()
            std_fieldName = fieldName
            if native_fld_name in self._field2file:
                file_name = self._field2file[native_fld_name]
        if file_name is None:
	        raise Exception(msg_unk.format(fld=fieldName, levType=levType,
                            inspec=self.inputspec))
        return file_name

    def get_field(self, fieldName, units, levType=None, overrideTransforms=None):
        """
        Get a field from its dataset. First, see if it is a derived field. If 
        not, attempt to treat `fieldName' as
        the native name in the dataset. If that fails, see if there is a 
        mapping in _native2standard. 
        
        preC: The Variable has been added to self._field2file (native 
        name to file mapping). If passing in the standard name, 
        self.standard2native must be set as well.
        :param levType: One of the level type string in self.known_level_types. 
        This is MANDATORY unless passing in the native name AND the field is 2d
        :return : ConformedField
        """
        # First, check if it is a derived variable and if so process it and return
        #import pdb ; pdb.set_trace()
        derived = False
        if fieldName in self._derived_fields:
            func = self._derived_fields[fieldName]
            derived = True
            (data, lats, lons, src_units) = self.__process_derived_field(fieldName)
            dim_order = get_dim_order(nioVar, self._var2dim, self._log) # tkji
            #TODO data_callback - Callback to function that creates the derived field
            # TODO : pass in the level_type ("isobaric", "sfcDelta", etc.)
            var = ConformedField(fieldName, fieldName, units, self, dim_order,
                                 self._field_transforms, levType, data=x,
                                 srcUnits=src_units, 
                                 derived=True, log=self._log)
            return var
        # Then, check if it is a multivar field
        elif levType != None and fieldName in self._multivar_fields[levType]:
            multifield = self._multivar_fields[levType][fieldName]
            # * TODO : Pass in missing value
            # * TODO : Passing in None as the native_name, but native_name 
            #          is needed for get_coords
            c = ConformedField(fieldName, None, units, self,     
                               multifield["dimOrder"], self.field_transforms, levType, 
                               #data=np.array(multifield["levData"]), <- current logic is to just get from 2d var
                               multivar=True,
                               log=self._log)
            return c

        """  
        file_name = self.get_filename(fieldName, levTye)

        # Not derived, so first figure out which file has the variable
        file_name = None
        if fieldName in self._field2file:
            file_name = self._field2file[fieldName]
            native_fld_name = fieldName
            # TODO ? set standard field name, which is passed to
            #     ConformedField constructor
            #if levType is not None:
                #try: 
                    #std_fieldName = self._standard2native
            self._log.info("Since native field name given, will pass None as "
                "standard_name to ConformedField constructor")
            std_fieldName = None
        else:
            # Maybe standard_name given
            if levType is None:
                if fieldName in self._standard2native["2d"]:
                    levType = "2d"
                else:
                    raise Exception("If passing `standard_name' to `get_field()'"
                                    ", you must pass in the levType (unless "
                                    "field is 2D). If it's 2D, ensure there's a"
                                    "mapping in inputspec [2d_field_mappings]")
            native_fld_name = self._standard2native[levType][fieldName]
            std_fieldName = fieldName
            if native_fld_name in self._field2file:
                file_name = self._field2file[native_fld_name]
        if file_name is None:
	        raise Exception(msg_unk.format(fld=fieldName, levType=levType,
                            inspec=self.inputspec))
        """
        file_name = self.get_filename(fieldName, levType)
        # Not derived, so first figure out which file has the variable
        nio_file = get_open_file(os.path.join(self.topdir, file_name))
        
        # Process the variable
        ## this is duplicate from get_filename
        if fieldName in self._field2file:
            native_fld_name = fieldName
            std_fieldName = None
        else:
            if levType is None: 
                if fieldName in self._standard2native["2d"]:
                    levType = "2d"
            native_fld_name = self._standard2native[levType][fieldName]
            std_fieldName = fieldName
        ##
        nio_var = nio_file.variables[native_fld_name]
        # * TODO : GriB and g5nr and g5nr_nps files have the _FillValue, 
        #          nemsio2netcdf files have nothing. Don't know about native WRF
        missing_val = None
        for attrName in ("_FillValue", "missing_value", "fmissing_value"):
            if hasattr(nio_var, attrName):
                missing_val = getattr(nio_var, attrName)
                break
        #missing_val = getattr(nio_var, "_FillValue")
        dim_order = get_dim_order(nio_var, self._var2dim, self._log) # tkji
        c = ConformedField(std_fieldName, native_fld_name, units, self, dim_order,
                           self._field_transforms, levType, log=self._log,
                           missingValue=missing_val)
        return c
    
    def __set_var2dim_mappings(self):
        """
        Populate self._var2dim using the items specified in the 
        [var2dim_mappings] section of the inputspec, which should consist of 
        mappings from a dimension variable to the dimension name.
        Both horizontal and vertical coordinates are specified in this section;
        they are distinguished by the fact that the horizontal coordinates 
        specify the lat and lon in the same line. 
        e.g. 
        HLAT,HLON = south_north, west_east
        PINT = bottom_top_staggered
        
        The self._var2dim variable is a dictionary with 3 keys:
         "lats", "lons", "levs"
        Each key maps to another dictionary that maps variable names to 
        dimension names for the corresponding coordinate.
        
        These mappings are optional for Variables for which the dimension and
        the Variable corresponding to it have the same name.
        """
        
        #self._var2dim = { "lats":{}, "lons":{}, "levs":{} } # in constructor
        horzmap = dict(self.inputspec.items("var2dim_mappings"))
        latVarName = None ; lonVarName = None
        for varstr,dimstr in horzmap.iteritems():
            if varstr.startswith("__"): continue # '__' prefix == interp var
            if "," in varstr:
                # item is a horizontal coordinate mapping
                latvar  = varstr.split(",")[0].strip()
                lonvar = varstr.split(",")[1].strip()
                latdim = dimstr.split(",")[0].strip()
                londim = dimstr.split(",")[1].strip()
                self._var2dim["lats"][latvar] = latdim
                self._var2dim["lons"][lonvar] = londim
                self._log.debug("Set horz var2dim map '{latvar}->{latdim}'"
                                 " and '{lonvar}->{londim}'"
                                 .format(latvar=latvar, latdim=latdim, 
                                         lonvar=lonvar, londim=londim))
            else:
                # item is vertical coord mapping
                self._var2dim["levs"][varstr] = dimstr
                self._log.debug("Set vert var2dim mapping '{var}->{dim}'"
                                 .format(var=varstr, dim=dimstr))
                
    def set_field_transform(self, k, v, excludedVars=None):
        """
        Add a transformation that should be applied whenever get_field is 
        called.
        Known transforms are:
            multiply_by: float - Multiply all data values by this factor
            convert_units: True/False - Convert data to the given units
            monotonic_grid: True/False - Shift the grid and coordinate data so 
                        it is monotonically increasing (e.g. goes from -180->180 
                        (instead of 0->360). This is required if plotting with
                        Basemap
            squeeze_out_time: True/False - Squeeze out the time dimension
        """
        self._log.info("Adding field transform {0}->{1}".format(k,v))
        self._field_transforms[k] = v

    def set_coord_transform(self, k, v):
        """
        Add a transformation that should apply to all HorizontalCoordinate 
        objects created from now on.
        KIM that HorizontalCoordinate objects are created dynamically 
        as field data is requested, but there is only one HorizontalCoordinate 
        for all fields on the same horizontal grid
        """
        self._log.info("Adding coord transform {0}->{1}".format(k,v))
        self._coord_transforms[k] = v
        
    
    def _update_field_transforms(self, transforms):
        self._field_transforms.update(transforms)
    
    @property
    def field_transforms(self):
        return self._field_transforms
    
    def get_coords(self, fieldName, levType):
        """
        Get the HorizontalCoordinate corresponding to the given `fieldName'.
        This requires that either:
          (A) The following sections of self.inputspec to be properly configured:
                [latlon2field_mappings] - Have map of lat,lon to fieldName
                [var2dim_mappings] - Map dimension's variable name (``dimvar'') 
                                     to the dimension name. Specifically, ensure 
                                     latVarName and lonVarName correspond
                                     to the dimensions
                [field2filename_mappings] - Ensure lat,lon variables have 
                                            mappings here (unless they are in 
                                            the same file as fieldName)
          (B) The latitude and longitude dimension names correspond have 
              "lat" and "lon" in their names (respectively) AND there are
              Variables in the dataset with the same name as the dimension name 
              (this should be the case for GriB files read in with PyNIO, at
               least for latlon grids)
        
        The [latlon2field_mappings] should consist of one or more
        items specifying the name of the lat,lon variables to use (lhs) and
        the variables to use them for (rhs). An asterisk on the rhs is used to 
        signify all variables not specified by other items.
        e.g. vlat,vlon = U,V 
             hlat,hlon = *
           would assign hlat and hlon as the coordinates for all fields except U and V
        This method will update self.__horizontal_coords with the 
        HorizontalCoordinate corresponding to `fieldName'.
        It will update self._field2coord, so that no redundant work is
         performed on subsequent calls.
        
       *UPDATE* For multivar fields, just use the first level's coordinate 
       information (i.e. assume it's the same for all levels)
       -> The _field2coord will also be mapped from the native name 
          corresponding to the first field.
      :param fieldName Native field name
      :param levType Level type (only used for multivar fields)
      :return: -- A HorizontalCoordinate object
        """
        # TODO : This assumes that the coords do not change with time (or that there
        # is no time dimension in the dataset). So for moving nests it is a problem.
        #import pdb ; pdb.set_trace()
        if levType is not None and fieldName in self._multivar_fields[levType]:
            field_name = self._multivar_fields[levType][fieldName]["fieldNames"][0]
            self._log.debug("Multivar field '{0}': will use coordinates "
                            "corresponding to native Variable '{1}'"
                            .format(fieldName, field_name))
        else:
            field_name = fieldName
        
        # TODO ? iS this being populated? ( I see a lot of duplicate log messages in here)
        #import pdb ; pdb.set_trace()
        if field_name in self._field2coord:
            return self._field2coord[field_name]
            
        fileName = self._field2file[field_name]
        dataset = get_open_file(os.path.join(self.topdir, fileName))
        field_dims = dataset.variables[field_name].dimensions

        # Get latVarName and lonVarName 
        # First, try to get from inputspec
        horzmap = dict(self.inputspec.items("latlon2field_mappings"))
        latVarName = None ; lonVarName = None
        #import pdb ; pdb.set_trace()
        for dimstr,fieldstr in horzmap.iteritems():
            fields = fieldstr.strip().split()
            #import pdb ; pdb.set_trace()
            if len(fields) == 1 and fields[0] == "*":
                latVarName = dimstr.split(",")[0].strip()
                lonVarName = dimstr.split(",")[1].strip()
                self._log.info("Found wildcard indicating that non-specified "
                                "fields will map to (latvar,lonvar) = "
                                "({lat},{lon})."
                                .format(lat=latVarName, lon=lonVarName))
            elif field_name in fields:
                latVarName = dimstr.split(",")[0].strip()
                lonVarName = dimstr.split(",")[1].strip()
                break 
            
        # If no mapping for field_name, see if there are variables in the file 
        # with the same names as the lat and lon dimensions
        # TODO ? Is there overlap between here and get_dim_order()?
        if latVarName is None or lonVarName is None:
            for dimvar in field_dims:
                if dimvar in dataset.variables:
                    if "lat" in dimvar.lower():
                        self._log.info("Assuming latitude data is in Variable:"
                                         "{0} for field '{1}'"
                                         .format(dimvar, field_name))
                        latVarName = dimvar 
                    elif "lon" in dimvar.lower():
                        self._log.info("Assuming longitude data in Variable:"
                                         "{0} for field '{1}'"
                                         .format(dimvar, field_name))
                        lonVarName = dimvar 
        # If still not found, FAIL    
        if latVarName is None or lonVarName is None:
            raise Exception("Couldn't determine lat/lon variables for "
                            "field '{fld}'. Ensure it is specified under "
                            "[latlon2field_mappings] in inputspec files "
                            "'{fil}'"
                            .format(fld=field_name, fil=self.inputspec_paths))
        
        # If we already have a HorizontalCoordinate for this (latvar,lonvar), use it
        latlon_key = latVarName + lonVarName
        if latlon_key in self.__horizontal_coords:
            coord = self.__horizontal_coords[latlon_key]
            self._log.debug("Using existing HorizontalCoordinate (lat,lon)"
                            "= ({lat},{lon}) for field {f}"
                            .format(lat=latVarName, lon=lonVarName,f=field_name))
            if "dataset" in locals(): 
                close_open_file(fileName)
            return coord
            
        self._log.debug("Creating HorizontalCoordinate for (lat,lon) "
                         "= ({lat},{lon}) for field {fld}"
                        .format(fld=field_name, lat=latVarName, lon=lonVarName))    
        # Get dimension names corresponding to the dimensino variables
        if latVarName in self._var2dim["lats"] \
          and lonVarName in self._var2dim["lons"]:
            lat_dim_name = self._var2dim["lats"][latVarName]
            lon_dim_name = self._var2dim["lons"][lonVarName]
            if not "dataset" in locals():
                dataset = get_open_file(self._field2file[field_name])
                #dataset = get_open_file(os.path.join(self.topdir, fileName))
                #dataset = Nio.open_file(self._field2file[field_name])
                vardims = dataset.variables[field_name]
                assert lat_dim_name in vardims
                assert lon_dim_name in vardims
        elif latVarName in field_dims and lonVarName in field_dims:
            lat_dim_name = latVarName
            lon_dim_name = lonVarName
        else:
            raise Exception("Could not determine lat/lon dimension names."
                            "Please specify mapping for coordinate variables "
                            "{lat} and {lon} in section [var2dim_mappings] of "
                            "inputspec file ('{inspec}') for field '{fld}'"
                         .format(fld=field_name, lat=latVarName, lon=lonVarName,
                                 inspec=self.inputspec_paths))
            
        close_open_file(dataset)
        # Get files containing lat vars
        # Try to get from field_name's file if below mappings dont exist
        try:
            lat_file = self._field2file[latVarName]
        except KeyError:
            lat_file = self._field2file[field_name]
        try:
            lon_file = self._field2file[lonVarName]
        except KeyError:
            lon_file = self._field2file[field_name]
        lats_dataset = get_open_file(lat_file)
        latVar = lats_dataset.variables[latVarName]
        latDims = latVar.dimensions
        if lat_file != lon_file:
            lons_dataset = get_open_file(lonFile)
        else:
            lons_dataset = lats_dataset
        lonVar = lons_dataset.variables[lonVarName]
        lonDims = lonVar.dimensions
        if not (lat_dim_name in latDims) or (not lon_dim_name in lonDims):
            raise Exception("Specified mapping for lat,lon variables '{lat},"
                            "{lon} to dimensions {latDim},{lonDim}, but these "
                            "dimensions do not exist in the variables, whose "
                            "dimensions are: [{latDims}],[{lonDims}]"
                            .format(lat=latVarName, lon=lonVarName, 
                                    latDim=lat_dim_name, lonDim=lon_dim_name,
                                    latDims=latDims, lonDims=lonDims))
        # Get staggering and data order for lat/lon 
        # ASSUME : The lat/lon var values are at same staggerloc as the field
        stagger = self._get_staggerloc(field_name, fileName)
        # Ensure order makes sense; account for possibility of there
        # being a time dimension
        #import pdb ; pdb.set_trace()
        lat_order = get_dim_order(latVar, self._var2dim, self._log)
        lon_order = get_dim_order(lonVar, self._var2dim, self._log)
        
        # squeeze out time
        # (TODO to support multi-time datasets - gotta keep the time dimension
        #  and pass in a time to get_coords
        if "t" in lat_order:
            latData = latVar[:].squeeze()
            #import pdb; pdb.set_trace()
            assert np.all(latData == latVar[0,:]) # assume T is first idx
            lat_order = lat_order.replace("t","")
        else:
            latData = latVar[:]
        if "t" in lon_order:
            lonData = lonVar[:].squeeze()
            assert np.all(lonData == lonVar[0,:])
            lon_order = lon_order.replace("t","")
        else:
            lonData = lonVar[:]

        # Sanity checks
        if len(lat_order) == 1 or (len(lat_order) == 2 and "t" in lat_order):
            assert lat_order.replace("t","") == "j" 
            assert lon_order.replace("t","") == "i"
        else:
            assert lat_order == lon_order
        self._log.debug("Dimension order for lat: {0}; lon: {1}"
                         .format(lat_order, lon_order))
        # Create the HorizontalCoordinate and add it to self.__horizontal_coords # Get latDimName and lonDimName
        #import pdb ; pdb.set_trace()
        coord = HorizontalCoordinate(latData, lonData, lat_file, lon_file, 
                                     latVarName, lonVarName, lat_dim_name, 
                                     lon_dim_name, stagger, self,
                                     latOrder=lat_order, lonOrder=lon_order,
                                     transforms=self._coord_transforms, 
                                     log=self._log)

        self._field2coord[field_name] = coord
        # Add mapping to global dict of HorizontalCoordinate corresponding to
        # this lat and lon values
        self.__horizontal_coords[latlon_key] = coord
        if lons_dataset is not lats_dataset: 
            close_open_file(lons_dataset)
        close_open_file(lats_dataset)
        return coord
    
    def __set_derived_fields(self):
        pass
        # TODO
        #self._derived = 
        
        
    def __process_derived_field(self, fieldName):
        """
        Process a derived field (``fieldName''). The function to use to generate
        the field and the dependencies are specified in the inputspec in the 
        section [fieldName_derivation].

        """
        field_spec = dict(self.inputspec.items(fieldName + "_derivation"))
        dep_fields = []
        for fldName in field_spec["dependencies"].split():
            fld = self.get_field(fldName)
            # TODO : ensure lats,lons,levs are same for all dep fields
            dep_fields.append(fld)
        data = func_name(dep_fields)
        units = field_spec["units"]
        return (data, lats, lons, units)

    def ___get_dim_names(self, fieldName):
        """ Maps each field to it's dimensions """
        coord = self.get_coords(fieldName)
        latDimName = coord.lat_dim_name
        lonDimName = coord.lon_dim_name
        horz_dims = self.get_horz_dim_name(fieldName)
        
    __dim_name_mappings = {}
        
        
    def __set_standard2native_mappings(self):
        """ 
        Create dict from mappings defined in inputspec [field_mappings] 
        Reads [arbitrary_field_mappings], [2d_field_mappings], as well as
              [<x>_field_mappings] where <x> are all the known level types
              (i.e. self.level_types)
        See inputspec documentation for details of how the mappings are defined
        """
        #self._standard2native = dict(self.inputspec.items("field_mappings"))
        lev_types = ['arbitrary', '2d']
        lev_types.extend(self.level_types)
        for levType in lev_types:
            try:
                items = self.inputspec.items(levType + "_field_mappings")
                self._standard2native[levType] = {}
            except KeyError:
                self._log.debug("No field mappings specified for field type "
                                " '{type}'".format(type=levType))
                continue
                
            #import pdb ; pdb.set_trace()
            for std,rhs in items:
                if std.startswith("__"): continue # '__' prefix == interp var
                # If given attributes instead of the native name itself 
                # (e.g. for GriB), get native name from one of the input files
                # Look for the field with matching attributes in 
                #[BASIC]input_file and in all input files specified in
                # [field2filename_mappings]
                if rhs.startswith("attr("):
                    #import pdb ; pdb.set_trace()
                    input_files = []
                    if std in self._field2file:
                        fil = self._field2file[std]
                        fil = os.path.join(self.topdir, fil)
                        input_files.append(fil)
                    else:
                        try:
                            fil = self.inputspec.get("BASIC", "file_name")
                            fil = fil.format(**self._date_conversion_args)
                            fil = os.path.join(self.topdir, fil)
                            input_files.append(fil)
                        except NoOptionError:
                            self._log.info("No 'file_name' specified in "
                                           "inputspec::[BASIC]. Will look in"
                                           "all files in field2filename_mappings"
                                           .format())
                    if self.inputspec.has_section("field2filename_mappings"):
                        items = self.inputspec.items("field2filename_mappings")
                        infiles = [x[1] for x in items if x[0][:2] != "__"]
                        infiles = [os.path.join(self.topdir, f) for f in infiles]
                        infiles = [f.format(**self._date_conversion_args) for f in infiles]
                        input_files.extend(set(infiles))
                    if len(input_files) == 0:
                        raise Exception("No data files specified in inputspec")
                    existing_input_files = []
                    native = None
                    for fil in input_files:
                        if not os.path.exists(fil):
                            self._log.warn("Specified input file '{0}' does not exist"
                                           .format(fil))
                            continue
                        existing_input_files.append(fil)
                        try:    
                            #import pdb ; pdb.set_trace()
                            native = get_native_name_from_var_attributes(rhs, fil, fmt=self.file_format, log=self._log)
                        except UnmatchedFieldException:
                            self._log.debug("Attirubte set {0} not found in file {1}"
                                            .format(rhs, fil))
                        if native is not None:
                            break
                    if native is None:
                        self._log.warn("Could not find any field with attributes '{0}'"
                                        " from the following files: {1}. Requests "
                                        "for field '{2}' will FAIL"
                                        .format(rhs, existing_input_files, std))
                elif rhs.startswith("multivar["):
                    # Add entry to self._multivar_fields 
                    # e.g. multivar[ t200.200.hPa, t500.500.hPa, ..]
                    #      0123456789                            -1
                    #import pdb ; pdb.set_trace()
                    field_descs = rhs[9:-1].split(",")
                    if not levType in self._multivar_fields:
                        self._multivar_fields[levType] = {}
                    levUnits=[] ; levVals=[] ; varnames=[] ; dimOrders=[]
                    for field_desc in field_descs:
                        [varname,levVal,levUnit] = field_desc.split(".")
                        varnames.append(varname.strip())
                        levVals.append(float(levVal))
                        levUnits.append(levUnit.strip())
                        #if "W" in varname: import pdb ; pdb.set_trace()
                        fil = get_open_file(self._field2file[varname.strip()])
                        var = fil.variables[varname.strip()]
                        dimOrders.append(get_dim_order(var, self._var2dim, self._log))
                        close_open_file(fil)
                    assert [levUnits[0] == l for l in levUnits]
                    assert [dimOrders[0] == l for l in dimOrders]
                    native = {"varNames":varnames, "levVals":np.array(levVals), 
                              "levUnits":levUnits[0], "dimOrder":dimOrders[0]}
                    self._multivar_fields[levType][std] = native
                    # NOTE: No mapping in _standard2native since rhs aint native
                else: 
                    # rhs is not an attr string nor multifield
                    native = rhs
                self._log.debug("Creating std->native mapping: {sf}->{nf}"
                                 .format(sf=std, nf=native))
                self._standard2native[levType][std] = native
                
                          
    def __set_field2lev_mappings(self):
        """
        Set the field2lev dict, based on parameters in the [field2lev_mappings]
        section of the inputsec.
        Expected specification format for each item is one or more 
        comma-separated fields to a 2-tuple consisting of the level's
        dimension name and the variable containing the dimension data, if any.
        e.g.: TT,QQ,HGT,QC,RH = lev,PL 
        """
        if not self.inputspec.has_section("field2lev_mappings"):
            self._log.info("No [field2lev_mappings] section in inputspec. "
                            "dimension characteristics should therefore be "
                            "self-described in data file.")
        field2lev = dict(self.inputspec.items("field2lev_mappings"))
        for varNames,dimTupple in field2lev.iteritems():
            if varNames.startswith("__"): continue
            #import pdb ; pdb.set_trace()
            dimName = dimTupple.split(",")[0].strip()
            dimVarName = dimTupple.split(",")[1].strip()
            if dimVarName == "-": dimVarName = None
            for varName in varNames.split(","):
                self.field2lev[varName] = (dimName, dimVarName)

    def __set_field2file_mappings(self):
        """
        Set the self._field2file dictionary, which maps native field names to
        the files/collections they belong to. This method will read the 
        [field2filename_mappings] section of self.inputspec and for each
        specified field, create a mapping in self._field2file. The expected
        format for each line is:
        FIELD_NAME = path
        The file will be relative to self.topdir (i.e. os.path.join(self.topdir, path)
        path may itself be a relative path, in which case it must be specified as:
        FIELD_ONE = pathjoin(path, to, collection)
         -> looks in self.topdir/path/to/collection
        EXPERIMENTAL FEATURE:
          - If the native name is not known before-hand (think: GriB read with PyNIO)
            you may specify the FIELD_NAME as  standard name and level type as:
              std(level_type:field_name) = path
              e.g. std(isobaric: air_temperature) = nmbprs_TMP.grb2
          -> preC: self.standard2native must exist and the mapping must be there

        Also, read all Variables in the dataset specified in [BASIC] parameter 
        ``file_name'' and add them to self._field2file.
        The following keys are interpolated as datetime:
            cal_date - The calendar date corresponding to self.init_date + self.fcst_offset
            init_date - Initialization date of the forecast (self.init_date)
        The following keys are interpolated as strings:    
            fhr - Forecast hour
            fmin - Forecast minute
            dom - Domain/Grid number
        The left-hand side may contain multiple variables. e.g.:
            FIELD_ONE, FIELD_TWO = path_{cal_date:%Y%m%d%H%M}z.nc4
        Since this is created before the standard2native mappings, 
        the native field name must be given. 

        UPDATE:
        :return std2file_mappings: list of 3-tuples containing fieldName and level type of standard field name and level type 
                      that must be mapped after settings self._standard2native and the corresponding path
        """
        
        # TODO (for supporting GriB files in separate collections)
        #   Gotta support attr(...) on LHS, since we map native
        #   name to file name
        
        def _get_path(var, collection):
            """ inner function to resolve path to collection, which may be a
                simple string or a 'pathjoin(...)'"""
            if collection.startswith("pathjoin("):
                         #            0123456789 foo,bar,baz)
                assert collection[-1] == ")"
                #import pdb ; pdb.set_trace()
                toks = collection[9:-1].split(",")
                toks = [tok.strip() for tok in toks]
                subpath = os.path.join(*toks)
                return os.path.join(self.topdir, subpath)
            else:
                return os.path.join(self.topdir, collection)
        #import pdb ; pdb.set_trace()
        interp_file = lambda f: f.format(**self._date_conversion_args)
        
        # Add mappings for all variables in the "file_name" in [BASIC]
        # TODO : Time this for g5nr and 5000x2500 datasets
        try:
            f = self.inputspec.get("BASIC", "file_name")
            f = interp_file(os.path.join(self.topdir, f))
            #import pdb ; pdb.set_trace()
            self._log.debug("Opening file {0}".format(f))
            # TODO : create wrapper that passes in self.file_format
            nio_file = get_open_file(os.path.join(self.topdir, f), fmt=self.file_format)
            for field in nio_file.variables:
                self._field2file[field] = f
            close_open_file(nio_file)
        except NoOptionError:
            self._log.info("No 'file_name' option in section [BASIC]. Assuming all "
                     "files are specified in [field2filename_mappings]"
                     .format())
        # If there is no [field2filename_mappings] section, use the 
        # "file_name" in [BASIC] for all fields in the file
        if not self.inputspec.has_section("field2filename_mappings"):
            self._log.info("No [field2filename_mappings] section. Assuming all "
                           " fields are in the basic file '{0}'"
                           .format(self.inputspec.get("BASIC", "file_name")))
            return

        # Process "[field2filename_mappings]"
        mappings = dict(self.inputspec.items("field2filename_mappings"))
        # NOTE Don't convert to dict since there may be duplicate vars

        # Mappings for which the standard name (and levtype) is given must be
        # processed afterward since standard2native has not been called
        # yet
        std2file_mappings = []

        for vars,collection in mappings.iteritems():
            if vars.startswith("__"): continue # '__' prefix == interp var
            collection = collection.strip()
            for var in vars.split(","): # for each field (LHS) in collection (RHS)
                var = var.strip()
                path = _get_path(var, collection) # inner func
                
                # If the standard name was given, get corresponding native name (later)
                #import pdb ; pdb.set_trace()
                if var.startswith("std("):
                    levType,fieldName =  [s.strip() for s in var[4:-1].split(".")]
                    std2file_mappings.append((levType,fieldName,path))

                # Check if field was already mapped in [BASIC] file_name
                # and override if path exists for current date
                path_now = path.format(**self._date_conversion_args)
                if var in self._field2file and self._field2file[var] != path:
                    if os.path.exists(path_now):
                        self._log.warn("Overriding file mapping for field '{0}'"
                                        .format(var))
                        self._field2file[var] = interp_file(path)
                    else:
                        self._log.warn("Not overrideing field2file mapping for "
                                       "field '{0}' since the given path '{1}'"
                                       "does not exist"
                                       .format(var, path_now))
                else:
                    self._field2file[var] = interp_file(path)

                self._log.debug("Field '{fld}' expected in file '{fil}'"
                                .format(fld=var, fil=self._field2file[var]))

        return std2file_mappings

    def _get_speclev(self, fieldName):
        """
        Create (if necessary) and map (if necessary) a SpecifiedVerticalLevel
        object corresponding to the fieldName. 
        Specifically:
         0. Determine the level type
         1. Determine the dimension corresponding to `fieldName' either
            from the file itself or from the inputspec (for files that are not 
             self-described). 
         2. Retrieve the data array
         3. Create the SpecifiedVerticalLevel and add it to self._specified_levels
         4. Map the field to the level
        :return speclev: The corresponding SpecifiedVerticalLevel
        """
        # This should only be done once as the level settings don't change
        if fieldName in self._field2speclev:
            return self._field2speclev[fieldName]
        
        if fieldName in self._derived_fields:
            self.__get_lev_for_derived_field(fieldName) # TODO
            return
        try:
            field_file_name = self._field2file[fieldName]
        except KeyError:
            self._log.critical("'{0}' does not have a mapping in field2file."
                                "You must call get_field() before field2speclev()"
                                .format(fieldName))
            sys.exit(12)
        #import pdb ; pdb.set_trace()
        # Get level dim and var information
        # TODO : We can do a sanity check against the 
        #        ConformedField's lev type param
        lev_type = self.__get_lev_type(fieldName, field_file_name)
        (dimname, dimvarname) = self.__get_levdim_for_field(fieldName, field_file_name)
        stagger = self._get_staggerloc(fieldName, field_file_name)
        field_nioFile = get_open_file(field_file_name)
        if dimvarname in field_nioFile.variables:
            #dimvar_file = get_open_file(field_file_name)
            dimvar_file = field_file_name 
        elif not dimvarname in self._field2file:
            raise Exception("Dimension variable '{0}' is not in same data file "
                            "as field ('{0}'), so you must create a field2file "
                            "mapping for it in inputspec"
                            .format(dimvarname, field_file_name))
        else:
            dimvar_file = self._field2file[dimvarname]
        close_open_file(field_file_name)
        # Add to global list of SpecifiedVerticalLevel if not there.
        # Otherwise, Create appropriate SpecifiedVerticalLevel
        key = (lev_type, stagger, dimname, dimvarname)
        if key in self._specified_levels:
            self._field2speclev[fieldName] = self._specified_levels[key]
        else:
            # TODO
            # BRAINSTORMING: Do we want to pass the levFile? This could be a 
            #  problem if we later support transforms
            # OTOH, we don't want to load the level data if not necessary, but
            # we would have to do so if we do not pass inthe file - unless
            # we use a callback to SpecData, which would add complexity
            # Also consider that the dim var may be derived, so it wouldn't
            # be in the file anyway
            if lev_type == "isobaric":
                lev_object = IsobaricLevel(dimname, stagger, levVar=dimvarname,
                                           fileName=dimvar_file)
            elif lev_type == "sfcDelta":
                lev_object = SurfaceDeltaLevel(dimname, stagger, 
                                               levVar=dimvarname, 
                                               fileName=dimvar_file)
            elif lev_type == "mslDelta":
                lev_object = MSLDeltaLevel(dimname, stagger, levVar=dimvarname,
                                           fileName=dimvar_file)
            elif lev_type == "soilDelta":
                lev_object = SoilDeltaLevel(dimname, stagger, levVar=dimvarname,
                                            fileName=dimvar_file)
            #       dimName, staggerLoc, levData=None, nioFile=None, levVar=None
        
            # Update state variables
            self._specified_levels[key] = lev_object
            self._field2speclev[fieldName] = self._specified_levels[key]

        return self._field2speclev[fieldName]
    
    def __get_levdim_for_field(self, fieldName, fileName):
        """
        Get the level dimension and dimension variable corresponding to the field
        `fieldName'. 
        The dimension should always specified in the Variable (be it NioVariable or 
        something similar), but the variable containing the dimension's data
        (if it has any) is not.  
        This method will first attempt to get the variable's dimension info from
        the inputspec. Failing that, it will try to get it based on assumptions
        of how the file is self-described.
        :return (dim,dimvar): 2-tuple consisting of the dimension name and the 
            variable holding the dimension data (which may be None if it is not 
            defined)
        """
        nioFile = get_open_file(os.path.join(self.topdir, fileName))
        var = nioFile.variables[fieldName]
        # Now determine the dimension variable (if there is one)
        try:
            #import pdb ; pdb.set_trace()
            self._log.debug("Getting field->dimvar mapping from inputspec")
            (dimname,dimvarname) = self.field2lev[fieldName]
            assert dimname in var.dimensions
            if dimvarname in self._field2file:
                fil = get_open_file(self._field2file[dimvarname])
                assert dimvarname in fil.variables
                close_open_file(fil)
            else:
                assert dimvarname in nioFile.variables
                self._log.debug("Adding mapping in feld2file for dim var '{0}'"
                                 .format(dim_var_name))
                self._field2file[dim_var_name] = fileName
        except KeyError:
            self._log.debug("Field '{fld}' has no specified vertical coord in "
                             "inputspec [field2lev_mappings]. Will attempt "
                             "to get mapping from Variable's attributes"
                             .format(fld=fieldName))
            #self._log.debug("Getting field->dimvar mapping from dataset")
            (dimname,dimvarname) = get_lev_dim_info_from_nio_file(fieldName, nioFile)
       
        return (dimname,dimvarname)
    
    def __get_lev_type(self, fieldName, nioFile):
        """
        Given a field name and file, determine it's level type.
        :param fieldName: The native field name. 
        :return type: String representing the type of level:
            This will match one of the level types in known_level_types:
            isobaric, sfcDelta, mslDelta, soilDelta
        """
        try:
            #return get_lev_type_from_attr(fieldName, nioFile, self.level_types)
            return get_lev_type_from_attr(fieldName, nioFile)
        except UnknownLevelTypeException:
            # not self-described, so the field should be in one of the 
            # field_mappings for the known level types. Note that it must
            # be the native field name since the standard names are not necessarily
            # unique across all level types
            # TODO ? have a native2std mapping?
            for levType in self.level_types:
                #if fieldName in self.standard2native[levType]:
                for std,native in  self._standard2native[levType].iteritems():
                    if native == fieldName:
                        self._log.debug("Treating field '{0}' is as type '{1}'"
                                        .format(fieldName, levType))
                        return levType
            raise Exception("Field '{0}' not in any of the known level_types"
                            .format(fieldName))

    def get_lev_idx(self, varName, levType):
        """
        Get index of 3d data array corresponding to a given level value
        """
        field = self.fields[levType][varName]
        nio_var = field.nio_var
        dims = nio_var.dimensions
        lev_var = self.__get_inputspec_lev(varName, levType)
        levIdxInDims = dims.index(lev_var)
        #assert dims[levIdxInDims] == 
        assert field.dimensions["k"] == lev_var # dimensions is dict with keys i,j,k
        # NOTE : ISBL can be different in GriB files (e.g. RH is ISBL0 and others are ISBL1)
    
    #def __getattribute__(self, attrib):
        #"""
        #Override the "variables" attribute functionality to try
        #to obtain from standar_variable mappings
        #"""
        #if attrib == "variables":
            #try:
                #return getattr(self._nio_file, attrib)
            #except KeyError:
                #try:
                    #var = _native2standard[XX]
                #except KeyError:
                    #var = _derived_fields[XX]
                #except:
                    #raise Exception("Unknown variable: XX")
    
    #def __setattr__(self, k, v):
        #return setattr(self._nio_file, k, v)
       
    def __get_level_data(self, nioFile, fieldName, levType):
        """
        Get level data for the given field, either from the 
        inputspec or from the input file attributes/metadata
        """
        field_name = fieldName
        if fieldName in self.standard2native[levType]:
            field_name = self.standard2native[fieldName][levType]

        try:
            self._data = self.__get_level_data_from_var(nioFile, field_name)
        except AttributeError:
            self._log.debug("Attempting to get level data from inputspec")
            self._data = self.__get_level_data_from_inputspec()

    def __get_level_data_from_var(self, nioFile, fieldName):
        """
        Get level data from a data file. This currently only works with PyNio 
        Open File's of GriB files, since they're self described. Theoretically 
        should work with any object that has "level_type" attributes, a level
        dimension and corresponding variable in the same nioFile with the same
        name (and the name should start with "lv"), and some other things I can't 
        think of right now.
        Current limitation is that the dimension order is assumed to be 'kji'
        :param fieldName: The native name of the field in `nioFile'
        :type fieldName: str
        """
        if isinstance(nioFile, str):
            nioFile = get_open_file(nioFile)
            must_release = True
        var = nioFile.variables[fieldName]
        vardims = var.dimensions
        if must_release:
            close_open_file(nioFile)
        # Get lev_type, which may be needed to determine units 
        if hasattr(var, "level_type"):
            # Grib 2
            lev_type = getattr(var, "level_type")
        elif hasattr(var, "level_indicator"):
            # Grib 1
            lev_type = int(getattr(var, "level_indicator"))
        else:
            raise  AttributeError("No 'level_type' or 'level_indicator'. "
                                  "Try getting level data from inputspec")
        
        
        # If only one level value in dataset, the field will be 2-d and 
        # it will have a 'level' attribute
        # TODO : may have a time dimension (or some other dimension, which 
        #  must be accounted for. For now, just assert
        if hasattr(var, "level"):
            assert len(vardims) == 2
            # data will be an ndarray of 1 (e.g. surface) or 2 (e.g. 
            # layer between 2 layers).
            data = getattr(var, "level") # ndarray
            units = get_units_from_grib_defaults(lev_type)
            self._log.debug("Assuming units of '{0}'".format(units))
            # ensure ji order  - can't use get_dim_order() since that 
            # needs var2dim mappings, which are typically only specified for
            # non-GriB
            dimvardims = dim_var.dimensions
            assert "lat" in dimvardims[0] and "lon" in dimvardims[1]
            order = "ji"
        else:
            assert len(vardims) == 3
            dimIdx = [d[:2] for d in vardims].index("lv") 
            # Get dimension variable. Variables with data on the grid points has 
            # var name == dim name. Level types that are between layers have two 
            # corresponding dimension variables with suffixes l0 and l1
            levTypeStr = get_lev_type_from_attr(fieldName, nioFile)
            if levTypeStr == "soilDelta": # TODO : there are others
                self._log.warn("Only returning one of the dimension variables "
                                "but there are two (min/max)")
                dim_var = nioFile.variables[vardims[dimIdx] + "_l0"]
                #dim_var2 = nioFile.variables[vardims[dimIdx] + "_l1"]
            else:
                dim_var = nioFile.variables[vardims[dimIdx]]
            units = getattr(dim_var, "units")
            data = dim_var[:]
            # ensure kji order. See also note above
            dimvardims = dim_var.dimensions
            assert "lv" in dimvardims[0] and "lat" in dimvardims[1] \
                    and "lon" in dimvardims[2]
            order = "kji"
        return (data, units, order)
  
    def _get_staggerloc(self, fieldName, fieldFile):
        """
        Populate self._field2stagger from inputspec section 
        [field2stagger_mappings] (if not already done)
        and/or from the level_indicator/level_type attributes of 
        Variable ``fieldName'' in ``fieldFile''
        :return stag: One of the constants defined in
                      ESMF.api.constants.StaggerLoc.
                      If specified in the inputspec, just the name of the 
                      constant is given.
                      If getting from Variable attributes, assume CORNER_VFACE for 
                      isobaric, sfcDelta, and mslDelta fields.
        """
        # Populate self._field2stagger if necessary
        if self._field2stagger is None:
            self._field2stagger = {}
            f2s = dict(self.inputspec.items("field2stagger_mappings"))
            #import pdb ; pdb.set_trace()
            for fields,staggerLoc in f2s.iteritems():
                if fields.startswith("__"): continue
                for field in fields.split(","):
                    staggerVal = getattr(ESMF.api.constants.StaggerLoc, 
                                         staggerLoc)
                    self._field2stagger[field.strip()] = staggerVal
            
        if fieldName in self._field2stagger:
            return self._field2stagger[fieldName]
        
        # Not in mappings (yet), try to determine from Variable
        # attributes (e.g. GriB level indicator)
        # NOTE : we assume grib here
        self._log.debug("field2stagger not specified, attempting to get from "
                        "dataset attributes")
        lev_type = get_lev_type_from_attr(fieldName, fieldFile)
        if lev_type in ("isobaric", "sfcDelta", "mslDelta"):
            self._log.info("Assuming stagger = CORNER_VFACE for {0} "
                            "field {1}".format(lev_type, fieldName))
            staggerVal = ESMF.api.constants.StaggerLoc.CORNER_VFACE
        elif lev_type in ("soilDelta"):
            self._log.info("Assuming stagger = CORNER_VCENTER for {0} "
                            "field {1}".format(lev_type, fieldName))
            return ESMF.api.constants.StaggerLoc.CORNER_VCENTER # TODO : verify
        elif lev_type == "2d":
            return ESMF.api.constants.StaggerLoc.CORNER
        else:
            raise Exception("Unknown level type. Try specifying in inputspec")
        self._field2stagger[fieldName] = staggerVal
        return staggerVal
   
    def get_esmf_grid_for_field(self, cslice):
        """
        Get an ESMF_Grid object for the given ConformedDataSlice.
        This method will create/update self._esmf_grid on demand based on 
        the stagger location of the given field.

        TROUBLE: staggerloc in inspec for 3d fields is the 3d stagger, but we
        create 2d Grid here - gotta flatten 
        """
        #import pdb ; pdb.set_trace()
        coords = cslice.coords
        #staggerloc = coords.staggerloc
        staggerloc = cslice.staggerloc_2d
        if self._esmf_grid is None:
            (lats,lons) = coords.latlons
            nx = cslice.data.shape[cslice.dim_order.index("i")]
            ny = cslice.data.shape[cslice.dim_order.index("j")]
            dims = np.array([nx-1,ny-1]) # * TODO : Set this correctly - had to substract one since esmf adds one for corner stagger
            self._esmf_grid = ESMF.Grid(dims, staggerloc=[staggerloc], 
                                        coord_sys=ESMF.CoordSys.CART)
            #max_index = np.array([nx, ny, nz])
            # grid = ESMF.Grid(max_index, staggerloc=[ESMF.StaggerLoc.CENTER_VCENTER], coord_sys=ESMF.CoordSys.SPH_DEG)

            # Populate coordinates. Data order for ESMF is always x,y,z and 2d
            #import pdb ; pdb.set_trace()
            gridX = self._esmf_grid.get_coords(0, staggerloc=staggerloc)
            if coords.lon_dim_order == "ij":
                gridX[:] = lons
            elif coords.lon_dim_order == "ji":
                gridX[:] = lons.T
            else:
                raise Exception("Unexpected lon_dim_order: {0}".format(coords.lon_dim_order))
            gridY = self._esmf_grid.get_coords(1, staggerloc=staggerloc)
            if coords.lat_dim_order == "ij":
                gridY[:] = lats
            elif coords.lat_dim_order == "ji":
                gridY[:] = lats.T
            else:
                raise Exception("Unexpected lat_dim_order: {0}".format(coords.lat_dim_order))
            #gridXCenter[...] = x_par.reshape(x_par.size, 1, 1)

        if self._esmf_grid.coords[staggerloc][0] is None:
            self._esmf_grid.add_coords(staggerloc)
            # * TODO : populate coords as above - gotta put in separate loop or inner function

        return self._esmf_grid 

    @property
    def dataset_resolution(self):
        """
        Get resolution, in degrees, for the dataset. This method will base
        its value on global attributes. You should only use it if unable
        to get from the ConformedField's get_resolution() method if wanting
        to get the resolution of a specific field. If there are multiple 
        input files, use only the first one
        :return : 2-tuple (dx, dy)
        """
        if self._dataset_resolution is not None:
            return self._dataset_resolution
        # Find the first existing file in self.field2file
        for fil in self._field2file.values():
            args = self._date_conversion_args.update({"fhr":0})
            fil = fil.format(args)
            if os.path.exists(fil):
                dataset = get_open_file(fil)
                break
        # Get resolution from global attribute, assuming it's the same for dx and dy
        for possible_attr_name in ["DPHD", "resolution"]:
            try:
                res = float(getattr(dataset, possible_attr_name))
                self._log.debug("Returning global resolution from attr '{0}'"
                                .format(possible_attr_name))
                close_open_file(fil)
                return (res,res)
            except AttributeError:
                pass
        # Get resolution from global attribute, assuming it is different for X and Y
        dx = -1
        dy = -1
        for possible_attr_name in ["dx", "dy", "DX", "DY", "Dx", "Dy", "Di", "Dj"]:
           try:
                res = getattr(dataset, possible_attr_name)
                self._log.debug("Returning global resolution from attr '{0}'"
                                 .format(possible_attr_name))
                found = True
           except AttributeError:
                found = False
           if found:
               if possible_attr_name[1] in ["i","I", "x", "X"]:
                dx = float(getattr(dataset, possible_attr_name))
               elif possible_attr_name[1] in ["j", "J", "Y", "y"]:
                dy = float(getattr(dataset, possible_attr_name))
               else:
                raise Exception("Couldn't determine if attribute is for "
                                "lat or lon resolution")
               if dx > -1 and dy > -1:
                    return (dx, dy)
            
        # Try to get resolution by calculating spacing between lat and lon
        # TODO : Rather than just the 4 items in the list below, try to to
        # get from file.dimensions (KIM : lat usually has "lat in name.
        # same for lon)
        for possible_var_name in ["lon", "longitude", "lons", "HLON"]:
            if possible_var_name in dataset.variables:
                self._log.info("Calculating resolution from lat variable "
                               "'{0}' data".format(possible_var_name))
                lon_var = dataset.variables[possible_var_name]
                # figure out dimension name(s) to see if its 1d or 2d
                # TODO : finish

        # if at this point we haven't returned anything, we were unable to 
        # determine the resolution
        raise Exception("Unable to determine dataset resolution")
    
    # SpecData::
    def get_horz_slice(self, specfield, levIdx, levValue, levUnits, 
                       timeIdx, fieldTransforms, sliceRetrievalSetting=None,
                        outDimOrder=None):
        """
        Helper method for getting a horizontal slice of a ConformedField.
        This is needed since the SpecifiedDataset manages the transformations
        that must be performed on field data as a result of transformations
        that were performed on the coordinate data.
        :param fieldName: Native field name 
        """

        # TODO (for derived fields): Since derived fields do not have an 
        # associated Variable, extended selection won't work...

        #if "t" in self.dim_order:
        #    raise Exception("Not implemented with current (simple) slicing scheme") # TODO
        msg_unk_lev = "Unknown level variable. Make sure is being passed to"\
                      " the constructor of `speclev'. Behavior is "\
                      " not implemented for derived level values"
        
        if sliceRetrievalSetting is None: 
            sliceRetrievalSetting = self.inputspec.get("TODO")
        
        # TODO : Check if field is derived. If so, we'll need to retrieve its
        #  data from the specfield, unless we go with approach menioned in
        # above TODO message
        if specfield.multivar:
            levType = specfield.level_type
            std = specfield.standard_name
            multivar = self._multivar_fields[levType][std]
            inUnits = multivar["levUnits"]
            levVals = multivar["levVals"]
            levVals = Units.conform(levVals, Units(inUnits), Units(levUnits))
            if levValue is None:
                raise Exception("Only levValue suported for multivar")
            try:
                idx = levVals.tolist().index(levValue)
            except IndexError:
                raise Exception("Need identical level value for multivar")
            native_field_name = multivar["varNames"][idx]
            self._log.info("Getting slice for multivar from native var '{0}'"
                            .format(native_field_name))
            
        else:
            native_field_name = specfield.native_name
        
        filName = self._field2file[native_field_name]
        nioFile = get_open_file(filName)
        var = nioFile.variables[native_field_name]
        
        # If there is a time index, set it for the PyNIO extended
        # selection string
        if timeIdx is None:
            fmt = ""
        else:
            #Get time dimension name
            vardims = var.dimensions
            # TODO ? allow for other time dimension names
            idx = [v.lower() for v in vardims].index("time")
            time_dim_name = vardims[idx]
            fmt = "{timeDim}|i{idx} ".format(timeDim=time_dim_name, idx=timeIdx)
        # TODO? I don't know if behavior is defined for PyNIO extended selection
        # if Pressure variable is 3-d
        
        # Get coordinates
        coords = self.get_coords(native_field_name, specfield.level_type)
                
        # Get 2-d slice of data
        if specfield.multivar:
            slice_dat = var[fmt]
            self._log.info("Slice retrieved")
        elif specfield.level_type != None and levIdx is not None:
            # level index given
            speclev = self._get_speclev(native_field_name)
            self._log.debug("level index given, returning data slice @ idx={0}"
                            .format(levIdx))
            # format string for PyNIO extended selection
            fmt += "{dim}|i{idx}" #  i => get by index
            fmt = fmt.format(dim=speclev.dim_name, idx=levIdx)
            #import pdb ; pdb.set_trace()
            self._log.info("Getting data slice using the following fmt: '{0}'"
                           .format(fmt))
            slice_dat = var[fmt]
            self._log.info("Slice retrieved")
        elif specfield.level_type != None and levValue is not None:
            # level value given.
            speclev = self._get_speclev(native_field_name)
            self._log.debug("Will get slice of data at level '{0}'"
                             .format(levValue))
            lev_dim_name = speclev.dim_name
            lev_var_name = speclev.var_name
            if lev_var_name is None:
                try:
                    lev_var_name = nioFile.variables[lev_dim_name]
                except KeyError:
                    raise Exception(msg_unk_lev)
                
            # Extended selection is not working for g5nr_nps_isobaric ...
            #if lev_dim_name == lev_var_name:
            #    fmt += "{dim}|{levVal}" # get by value
            #else:
            #    fmt += "{dim}|{var}|{levVal}"
            # TODO : Get units from inspec if necessary
            inUnits = getattr(nioFile.variables[lev_var_name], "units")
            available_levs = nioFile.variables[lev_var_name][:]
            # Conform level's units
            available_levs = Units.conform(available_levs,
                                           Units(inUnits), Units(levUnits))
            try:
                #First try to get exact value
                idx = available_levs.tolist().index(levValue)
                self._log.info("Found exact vertical level = {0} {1}"
                               .format(levValue, levUnits))
                fmt += "{dim}|i{idx}".format(dim=lev_dim_name, idx=idx)
                slice_dat = var[fmt]
            except IndexError:
                if sliceRetrievalSetting == "exact":
                    # Not working for g5nr_nps_isobaric
                    #fmt = fmt.format(dim=lev_dim_name, var=lev_var_name, 
                    #                 levVal=levValue) 
                    #import pdb ; pdb.set_trace()
                    raise Exception("Wanted vertical level = {0} {1}"
                                    "Available levels = {2} {1}"
                                    .format(levValue, levUnits, available_levs))
                elif sliceRetrievalSetting == "nearest":
                    nearest_val = speclev.get_nearest(levValue, levUnits)
                    #fmt = fmt.format(dim=lev_dim_name, levVal=nearest_val) 
                    idx = available_levs.tolist().index(nearest_val)
                    fmt += "{dim}|i{idx}".format(dim=lev_dim_name, idx=idx)
                    slice_dat = var[fmt]
                elif sliceRetrievalSetting == "interpolate":
                    # TODO
                    raise Exception("Not implemented - next version")
                else:
                    raise Exception("Unknown slice_retrieval_setting")
        elif specfield.level_type is not None:
            raise Exception("Must pass in levVals or levIdx for 3d Variables")
        else:
            # no vertical levels
            speclev = None # TODO ? Have a speclev for 2d or use specific type like "MSL"
            self._log.info("Getting data slice using the following fmt: '{0}'"
                           .format(fmt))
            slice_dat = var[fmt]


        # Convert to specfield's units
        # TODO ? Move to perform_transforms_on_slice()
        # -> why are we only doing it for 3d, given-lev units anyway?
        try:
            in_units = getattr(var, "units")
        except AttributeError:
            try:
                in_units = getattr(var, "Units")
            except AttributeError:
                # TODO : Move this to a separate method
                units_map = dict(self.inputspec.items("units"))
                if native_field_name in units_map:
                    in_units = units_map[native_field_name]
                elif specfield.standard_name in units_map:
                    self._log.warn("Using inputspec specified units for "
                                   "standard field name {0}. This is not "
                                   "recommended since there could be multiple "
                                   "native fields(of different lev types) with "
                                   "the same standard name"
                                   .format(specfield.standard_name))
                    in_units = units_map[specfield.standard_name]
                else:
                    raise Exception("Unable to determine units of field '{0}'. "
                                    "If the Variable does not have a 'units' or"
                                    "'Units' attribute, you must specify it "
                                    "manually in the inputspec [units] section"
                                    .format(native_field_name))
        out_units = specfield.out_units
        
        # change units to names cfunits knows
        if in_units == "gpm": in_units = "m"
        if out_units == "gpm": out_units = "m"

        if in_units != out_units:
            slice_dat = Units.conform(slice_dat, Units(in_units), Units(out_units))
        
        # Synchronize data
        #import pdb  ; pdb.set_trace()
        # TODO ? Support for level data transforms - will be needed if support for
        # 3d lev values (e.g. for PINT)
        #if specfield.level_type != None:
        #    speclev.perform_transforms(self._lev_transforms)
        # Note that this will constantly overwrite the coord data if different 
        # coord_transforms are pushed from Field/lev, but currently we are not 
        # doing this - only the Coordinates transforms are done - and thos are
        # only done once 
        coords.perform_transforms(self._coord_transforms)
        # Note that if the follow call updates coords' transforms, the data will 
        # need to be synchronized again, so this may need to be recursive somehow
        slice_dat, fieldTransforms = specfield.perform_transforms_on_slice(fieldTransforms, slice_dat)
        
        # Transpose slice if outDimOrder != native dim order
        if outDimOrder is not None:
            if specfield.multivar:
                # use coords for corresponding 2-d field, which was defined above
                latIdx = var.dimensions.index(coords.lat_dim_name)
                lonIdx = var.dimensions.index(coords.lon_dim_name)
            else:
                # TODO ? Necessary ? Can't we just use coords defined above?
                latIdx = var.dimensions.index(specfield.coords.lat_dim_name)
                lonIdx = var.dimensions.index(specfield.coords.lon_dim_name)
            if latIdx > lonIdx:
                native_order = "ij"
            else:
                native_order = "ji"
            if native_order != outDimOrder:
                slice_dat = slice_dat.T

        if specfield.multivar:
            speclev = None
        # Return the slice
        return (slice_dat, coords, speclev)

class ConformedDataSlice(object):
    """
    Class to encapsulate data that has already been conformed 
    """
    def __init__(self, nativeFieldName, fieldUnits, levType, levIdx, levValue, 
                 levUnits, timeIdx, transforms, dimOrder, retrievalSetting, 
                 conformedField):
        self.native_name = nativeFieldName
        self.field_units = fieldUnits
        self.lev_type = levType 
        self.lev_idx = levIdx
        self.lev_value = levValue
        self.lev_units = levUnits 
        self.time_idx = timeIdx
        self.transforms = transforms
        self.dim_order = dimOrder
        self.retrieval_setting = retrievalSetting
        self._esmf_field = None
        self._conformed_field = conformedField
        stagger3d = conformedField.get_coords(levValue, levUnits).staggerloc
        self.staggerloc_2d = get_2d_staggerloc_from_3d_staggerloc(stagger3d)
        self._data = None
        self._coords = None
        self._speclev = None

    def __eq__(self, other): 
        """
        To save time, only compare the attributes, not the data. If attributes
        are the same the data should be the same
        """
        d = copy.deepcopy(self.__dict__)
        d.pop("_data")
        d.pop("_coords")
        d.pop("_speclev")
        return d == other.__dict__
    
    @property
    def data(self):
        if self._data is None:
            raise Exception("Haven't called get_horz_slice() on SpecField object")
        return self._data
    @data.setter
    def data(self, data):
        assert data.ndim == 2
        self._data = data
    
    @property
    def coords(self):
        if self._coords is None:
            raise Exception("Haven't called get_horz_slice() on SpecField object")
        return self._coords
    @coords.setter
    def coords(self, coords):
        self._coords = coords
    @property 
    def speclev(self):
        if self._speclev is None:
            raise Exception("Haven't called get_horz_slice() on SpecField object")
        return self._speclev
    @speclev.setter
    def speclev(self, value):
        self._speclev = value

    def get_esmf_field(self):
        if self._esmf_field: 
            return self._esmf_field 
        fld = self._conformed_field
        #staggerloc = fld.coords.staggerloc
        esmf_grid = fld._specdata.get_esmf_grid_for_field(self)
        esmf_field = ESMF.Field(esmf_grid, staggerloc=self.staggerloc_2d)
        if self.dim_order == "ij":
            esmf_field.data[:] = self.data 
        elif self.dim_order == "ji":
            esmf_field.data[:] = self.data.T
        else:
            raise Exception("Invalid dim order for slice")
        #esmf_field.dataI#
        self._esmf_field = esmf_field
        return esmf_field


class ConformedField(object):
    """
    Encapsulates data access for a Variable whose data may be conformed in some way.
    e.g. converted to different units, shifted, multiplied by a factor, etc.
    :param standardName: The standard name assigned to the variable (from the inputspec; e.g. "air_temperature")
    :param nativeName: The native name of the variable in the dataset (or the value automatically assigned by PyNIO). For derived fields, this should be None
    :param specData: SpecifiedDataset associated with this ConformedField. It will be used for mapping coordinates to the field when retrieving data.
    :param units: 
    :param multivar: Is this a 3-d field made up of individual 2-d Variables in the dataset?
    :param levType: Level type indicator string: None or "2d" for 2d, "isobaric", 
                   "sfcDelta", "mslDelta", "soilDelta"
    TODO  :: Some methods assume there is a native_name. But this is not the case for derived_fields
    """
    def __init__(self, standardName, nativeName, outUnits, specData, 
                 dimOrder, varTransforms, levType, data=None, derived=False, multivar=False,
                 log=None, missingValue = None): # TODO ? Make missingVal mandatory
        self._log = log if log is not None else _default_log()
        self.out_units = outUnits
        self.standard_name  = standardName
        self.native_name = nativeName
        self._specdata = specData
        self._transforms = varTransforms
        self.derived = derived
        self.multivar = multivar
        self.dims_order = dimOrder
        if levType == "2d": levType = None
        self.level_type = levType
        if data is not None:
            self._data = data
        self.missing_value = missingValue

        # have lats,lons already been conformed?
        self.__coords_conformed = False
        # Keep track of what transforms have been performed to which rows
        # of data
        # TODO : Modify self._log message to have prefix of self.name

        ##
        # Other state variables
        ##
        """
        dict that maps conformed data slices obtained from get_horz_slice() 
        from the settings that were used to conform the data. This way, 
        if data is requested multiple times, the data is not re-conformed.
        Also, the files do not need to remain open since new data is created.
        """
        self._conformed_data_slices = {}

        """
        2-tuple containing (longitudinal grid resolution, latitudinal grid resolution),
        in degrees. Not for public access. Use self.resolution
        """
        self._resolution = None

    def _get_native_name(self, levValue=None, levUnits=None):
        # Get native variable name
        if self.multivar:
            #import pdb ; pdb.set_trace()
            multivar_tup = self._specdata._multivar_fields[self.level_type][self.standard_name]
            idx = multivar_tup["levVals"].tolist().index(levValue)
            assert levUnits == multivar_tup["levUnits"]
            native_name = multivar_tup["varNames"][idx]
        else:
            native_name = self.native_name
        return native_name

    @property
    def coords(self):
        """
        Returns the HorizontalCoordinate associated with this Field
        """
        # NOTE : native_name is None for multivar; caller should know
        # how to handle it
        self._log.warn("coords is deprecated use get_coords()")
        if self.multivar: raise Exception("Must use get_coords")
        return self._specdata.get_coords(self.native_name, self.level_type)
    def get_coords(self, levValue=None, levUnits=None):
        """
        Returns the HorizontalCoordinate associated with this Field.
        Compatible with multivar fields
        """
        #import pdb ; pdb.set_trace()
        native_name = self._get_native_name(levValue, levUnits)
        return self._specdata.get_coords(native_name, self.level_type)

    def _conform_array(self):
        """
        Conform an array by performing any of the known transforms specifiedin 
        self._var_transforms
        TODO : Do we want this here or in an external function in a "transforms" module ??
        """
        data = self.__data
        if "multiply_by" in self._var_transforms:
            factor = self._var_transforms["multiply_by"] 
            if factor is not None:
                data *= factor
        else:
            self._log.debug("Key 'multiply_by' not specified")
        if "squeeze_out_time" in self._var_transforms:
            squeeze = self._var_transforms["squeeze_out_time"]
            if squeeze is True:
                if "t" in self.dims_order:
                    shape_orig = data.shape()
                    data = data.squeeze()
                    if shape_orig == data.shape:
                        raise Exception("Unable to squeeze time dimensions")
        else:
            self._log.debug("No 'squeeze_out_time' specified")
        if "units" in self._var_transforms:
            out_units = self._var_transforms["units"]
            in_units = None # TODO ?
            if out_units != in_units:
                data = Units.conform(data, Units(in_units), Units(out_units))
        else:
            self._log.debug("Key 'units' not specified")
        
        if "horz_roll" in self._var_transforms:
            shift = False
            shift = bool(self._var_transforms["monotonic_grid"])
        else:
            self._log.debug("Key 'horz_roll' not specified. Will not shift")
                
        return data

    #@property
    #def data(self):
        #"""
        #Perform any of the known transforms specifiedin self._var_transforms
        #on self.__data and return the transformed data
        #"""

        ## TODO : Don't want to do all data since we often only need a 
        ##        particular level (or maybe just 3 levels)

        #if np.all(self.__data_conformed): # data_conformed is an array specifying if each row has been conformed
            #return self._data

        #data = self.__data
        #data = self._conform_array(data)
        #self.__conformed[:] = True
        
        #self.__conformed = True
        #return self.__data

    def _get_data_in_order(self, expectedOrder):
        """
        Return self.data, possibly transposed depending on the value of
        `expectedOrder', which should be a string (e.g. "ij", "kji", etc.)
        """
        if len(expectedOrder) == 2 and expectedOrder in ("ij", "ji"):
            if self.data_order == order:
                return self.data
            elif self.data_order[::-1] == order:
                return self.data.T
            else:
                raise Exception("Unexpected data order")

#   ConformedField::
    def get_horz_slice(self, levValue=None, levUnits=None,  
                       levIdx=None, outDimOrder=None, 
                       #lookupSetting=None,
                       sliceRetrievalSetting="exact",
                       fcstOffset=None):
        """
        Retrieve a horizontal slab of data from self.data. The arguments set 
        depend on the type of level wanted.
        :param units: Units given for the level to be retrieved. e.g. to get 
                      the slab corresponding to 800 hPa, units="hPa" 
                      and levValue=800.
                      HINT: Do not use "mb", specify "millibars" or "hPa"
        :type units: str
        :param fcstOffset: The time offset from analysis time wanted. 
                          This is only needed if the dataset has multiple times
                          AND you want some time besides the analysis time.
        :type fcstDate: datetime.timdelta
        :param isobaricLev: Isobaric Level to get slab for
        :param outDimOrder: Dimension order (either "ij" or "ji"). If not passed
                            in, use native order
        :return -- The slab of data and it's coordinates. Ensures the transformations are consistent with the coordinates
        """
        
        native_name = self._get_native_name(levValue, levUnits)
        
        # See if any global transforms have been added (e.g. by the 
        # HorizontalCoordinate, to do horz_shift) for this variable
        self._transforms.update(self._specdata._field_transforms)

        # * TODO : Instead of getting the time slice here, just get the
        # time index and then use Nio extended selection or whatever
        # also: ensure time slice is in the ConformedDataSlice
        # time index will need to be padded to _specdata.get_horz_slice 
        if "t" in self.dims_order:
            time_idx = self.get_time_slice_idx(fcstOffset=fcstOffset)
        else:
            #assert data.ndim in (2,3)
            time_idx = None
        
        if outDimOrder is None:
            outDimOrder = self.dims_order.replace("t","").replace("k","")

        # If slice has already been retrieved, return it
        # TODO : derived_field_name for unique hash of derived field?
        #import pdb ; pdb.set_trace()
        slice_key = ConformedDataSlice(native_name, self.out_units,
                                       self.level_type,
                                       levIdx, levValue, levUnits,
                                       time_idx, self._transforms,
                                       outDimOrder,
                                       sliceRetrievalSetting, self)
        if slice_key in self._conformed_data_slices:
            self._log.debug("Using previously-requested slice with settings:"
                             " '{0}'"
                             .format(slice_key))
            return self._conformed_data_slices[slice_key]            
        
        
        self._log.debug("Getting data slice with the following specs: '{0}'"
                        .format(self._transforms))

        if "k" in self.dims_order and levIdx is None and levValue is None:
            raise ValueError("Must specify levIdx OR levValue")
        
        
        (slice_dat, coords, speclev) = \
            self._specdata.get_horz_slice(self, levIdx, levValue, levUnits, 
                                          time_idx, self._transforms, 
                                          sliceRetrievalSetting, outDimOrder)
        """
        if levIdx is not None:
            slice_dat = self.specdata.get_horz_slice(self, levType, levIdx,
                             levValue, levUnits, sliceRetrievalSetting)
        #if levType in (None, "-"):
                #self._log.debug("I'm 2-d, returning all of me")
                #return self.__get_data_in_order(order)
        elif levType == "index":
            slice_dat = ... # probably want to use the lev class to encapsulate 
                    # different possible positions of the lev dim
    
        elif levType in ("isobaric", "sfcDelta"):
            #native_name = self._get_native_name(levValue, levUnits)
            lev = self._spec_data.field2lev[native_name] # TODO : field2lev in specdata
            var = self._spec_data.variables[native_name]
            slab =lev.get_field_data_slab(var, levValue, lookupSetting,
                                        slabRetrievalSetting)
        else:
            raise Exception("Unknown level type '{0}'".format(levType))
        """
        
        # Add to 'database' in case the data is requested again.
        self._conformed_data_slices[slice_key] = (slice_dat, coords, speclev)
        
        # Populate data in ConformedDataSlice, so it can return the ESMF_Field
        slice_key.data = slice_dat
        slice_key.coords = coords
        slice_key.speclev = speclev

        return slice_key
        
    @property 
    def resolution(self):
        if self.multivar: raise Exception("Must use get_resolution")
        self._log.warn("resolution is deprecated. use get_resolution()")
        return self.get_resolution()
    def get_resolution(self, levValue=None, levUnits=None):
        """
        Determine resolution for the Field.
        ASSUME : values are degrees
        TODO : lots of testing
        """
        if self._resolution is not None:
            return self._resolution
            
        native_name = self._get_native_name(levValue, levUnits)

        srcFile = self._specdata._field2file[native_name]
        
        dataset = get_open_file(srcFile)
        var = dataset.variables[native_name]
        # Try to get from Variable's attributes 
        dx = -1.
        dy = -1.
        for possible_attr_name in ["Di", "dx", "Dx", "DX"]:
            try:
                lonvar = dataset.variables[self.coords.lon_var_name]
                dx = float(getattr(lonvar, possible_attr_name))
                assert getattr(lonvar, "units") == "degrees_east"
                self._log.info("Determining resolution from Variable attribute"
                               " '{0}'".format(possible_attr_name))
            except AttributeError:
                pass # try next
        for possible_attr_name in ["Dj", "dy", "Dy", "DY"]:
            try:
                latvar = dataset.variables[self.coords.lat_var_name]
                dy =  float(getattr(latvar, possible_attr_name))
                assert getattr(latvar, "units") == "degrees_north"
                self._log.info("Determining resolution from Variable attribute"
                                 " '{0}'".format(possible_attr_name))
            except AttributeError:
                pass # try next
        if dx > 0 and dy > 0:
            self._resolution = (dx,dy)
            return (dx, dy)

        # Try to get from lat&lon values
        #import pdb ; pdb.set_trace()
        (lats,lons) = self.coords.latlons
        if lats.ndim == 1 and lons.ndim == 1:
            cenlat = len(lats) / 2
            cenlon = len(lons) / 2
            dy = abs(lats[cenlat+1] - lats[cenlat])
            dx = abs(lons[cenlon+1] - lons[cenlon])
            assert abs(dy - dx) < 0.001
            self._log.info("Calculated resolution: (dx,dy) = ({0},{1})"
                           .format(dx,dy))
            self._resolution = (dx,dy)
            return (dx,dy)
        elif lats.ndim == 2:
            # HorizontalCoordinate class will put lat dim before lon dim if making 2d
            assert self.coords.lat_dim_order.index("j") < self.coords.lat_dim_order.index("i")
            assert self.coords.lon_dim_order.index("j") < self.coords.lon_dim_order.index("i")
            latIdx = self.coords.lat_dim_order.index("j")
            lonIdx = self.coords.lat_dim_order.index("i")
            cenlat = lons.shape[latIdx] / 2
            cenlon = lons.shape[lonIdx] / 2
            if latIdx > lonIdx:
                dy = abs(lats[cenlon,cenlat] - lats[cenlon,cenlat+1])
            else:
                dy = abs(lats[cenlat,cenlon] - lats[cenlat+1,cenlon])
            # repeat for lons
            latIdx = self.coords.lon_dim_order.index("j")
            lonIdx = self.coords.lon_dim_order.index("i")
            cenlat = lons.shape[latIdx] / 2
            cenlon = lons.shape[lonIdx] / 2
            if latIdx > lonIdx:
                dx = abs(lons[cenlon,cenlat] - lons[cenlon+1,cenlat])
            else:
                dx = abs(lons[cenlat,cenlon] - lons[cenlat,cenlon+1])
            self._log.info("Calculated resolution: (dx,dy) = ({0},{1})"
                           .format(dx,dy))
            self._resolution = (dx,dy)
            assert dx > 0.001 and dy > 0.001
            return (dx,dy)
        else:
            raise Exception("Unexpected dimensions")

        # If at this point we have not returned the resolution, we have failed
        raise Exception("Unable to determine resolution")
        
    def perform_transforms_on_slice(self, transforms, slice_dat):
        # TODO ? Should this be public
        alreadyPerformed = {}
        for xformName,xformValue in transforms.iteritems():
            # NOTE: For the horz_roll transform, gotta make sure you're 
            # rolling along the "i" axis. Since `slice_dat' is 2d, 
            # ignore k and t dims
            # ** TODO
            #import pdb ; pdb.set_trace()
            #shift_axis = self.dims_order.index("i")
            if xformName == "horz_roll":
                assert slice_dat.ndim == 2
                if self.dims_order.index("i") < self.dims_order.index("j"):
                    shift_axis = 0
                else:
                    shift_axis = 1
                self._log.debug("Assuming data has not been shifted")
                #import pdb ; pdb.set_trace()
                slice_dat,performed = shift_data(transforms, alreadyPerformed, slice_dat, 
                                                 shiftAxis=shift_axis, log=self._log)
                alreadyPerformed.update(performed)
            elif xformName == 'convert_units':
                # This should have already been performed (TODO : (v) )
                self._log.debug("Ignoring 'convert_units' transform")
            elif xformName == "equation":
                #import pdb ; pdb.set_trace()
                fn = Expression(xformValue, ["fieldSlice"])
                slice_dat = fn(slice_dat)
                alreadyPerformed["equation"] = transforms["equation"]
            elif xformName == "multiply_by":
                self._perform_multiply_xform()
            else:
                raise Exception("Uknown transform: {0}".format(xformName))

        transforms.update(alreadyPerformed)
        return slice_dat,transforms 
    
    def get_raw_data_arrays(self):
        """ 
        Return a 3- or 4- tupple consisting of the data, lats, lons,
        and (if 3d) level values. Each element in the tuple is an 
        an array
        """
        # TODO
        raise Exception("Not implemented")
    
    def get_raw_data(self):
        """
        Like get_raw_data_arrays() but return a Nio Variable for each
        of the 3 or 4 fields (field, lat, lon, lev). For derived fields,
        the corresponding element in the tuple will be None (since there 
        is no Variable associated with it
        """
        # TODO
        raise Exception("Not implemented")
    
    def __get_data_in_order(self, slab, order):
        """
        Convert 2d slab of data to desired order
        TODO
        """

    def get_time_slice_idx(self, fcstOffset=None):
        """
        This doesn't really do anything. It is a placeholder for future functionality
        that supports retrieving specific time offsets from data containing
        a (non-flat) time dimension
        :return: index along time dimension corresponding to fcstDate
        """
        # TODO : This is not really implemented. Assuming the data
        #        dimension is flat
        tIdx = self.dims_order.index("t")
        if fcstOffset is not None:
            # TODO : for this to work, need to passin some information about the
            # time dimension
            raise Exception("Not implemented")
            #i = get_index_for_time(tIdx, fcstOffset, self._log)
            if tIdx == 0:
                #data = data[i,:]
                return tIdx
            # Just returning index here so we don't have to read in
            # the entire dataset
            #elif tIdx == len(self.dims_order) - 1:
            #    if data.ndim == 3:
            #        data = data[:,:,i]
            #    elif data.ndim == 4:
            #        data = data[:,:,:,i]
            #    else:
            #        raise Exception("Unexpected dimensionality")
            else:
                raise Exception("'time' dimension must be first or last")
        else: # fcstDate not passed in
           return 0
           # if tIdx == 0:
           #     return tIdx
                #data = data[0,:]
           # else:
           #     raise Exception("Unexpected dimensionality")


class HorizontalCoordinate(object):
    """
    This class encapsulates functionality specific to the horizontal coordinate
    system pertaining to one or more variables in a SpecifiedForecastDataset.
    It's main functionality is providing an interfrace to the lat,lon data while
    supporting the ability to perform transforms to the data and letting the 
    SpecifiedForecastDataset object know of corresponding transforms needed 
    for the Fields that correspond to this coordinate system.
    """
    def __init__(self, lats, lons, latFileName, lonFileName, latVarName,
                 lonVarName, latDimName, lonDimName, stagger, specdata, 
                 transforms={}, latOrder=None, lonOrder=None, log=None):
        """
        Instantiate a HorizontalCoordinate with the given lat/lon values.
        The lat/lon source file and Variable names are given in case
        attributes need to be looked up, but they should __not__ be used
        to access data, since the data may be conformed using transforms.
        
        :param latFileName: The file containing the latitude data. 
        :type latFileName: str 
        :param lonFileName: Same as latNioFile, for longitude
        :type lonFileName: str or Nio.OpenFile
        :param specdata: Dataset that I am bound to. Fields and coordinates are 
                        synchronized through the SpecifiedDataset
        :type specdata: SpecifiedForecastDataset
        """
        self._log = log if log is not None else _default_log()
        #import pdb ; pdb.set_trace()
        self.specdata = specdata
        self.__lat_file_name = latFileName
        self.__lon_file_name = lonFileName
        self.lat_var_name = latVarName
        self.lon_var_name = lonVarName
        self.lat_dim_name = latDimName
        self.lon_dim_name = lonDimName
        self.staggerloc = stagger
        self.lat_dim_order = latOrder
        self.lon_dim_order = lonOrder
        assert not "k" in latOrder 
        assert not "k" in lonOrder
        
        self._transforms = {}
        self._performed_transforms = {}
        self._lats = lats[:]
        self._lons = lons[:]
        self.__is_latlon_grid = False
        if self._lons.ndim == 1:
            self.__is_latlon_grid = True
         
        # Keep track of when transforms must be performed
        self._transforms_dirty = True
        # Map known transforms to the corresponding methods to call
        #import pdb ; pdb.set_trace()
        self._known_transforms = {
            "pm_relative_longitude": self._perform_pm_rel_lon_xform, 
            "make_2d": self._perform_2d_latlon_xform,
        }
        for k,v in transforms.iteritems():
            self.set_transform(k,v)
        
    def set_transform(self, name, value):
        """
        Update self._transforms with the given name,value pair
        """
        if not name in self._known_transforms:
            raise Exception("Unknown transform")
        lon_idx = self.lon_dim_order.index("i")
        lat_idx = self.lat_dim_order.index("j")
        try:
            currValue = self._transforms[name]
        except KeyError:
            prevValue = None
        self._transforms[name] = value 
        if value != prevValue:
            self._log.info("Updating transform {name} from {orig} to {new}"
                            .format(name=name, orig=prevValue, new=value))
            
        self._transforms_dirty = True
        
    def _perform_transforms(self):
        """Perform any transforms to the lat and lon values"""
        if not self._transforms_dirty:
            return
        self.perform_transforms(self._transforms)
        self._transforms_dirty = False
        # TODO ? - eliminate this - either just leave it up to
        # caller to manually perform them (as necessary) or use this
        # global thing that is automatically called when data is
        # requested

    def perform_transforms(self, transforms):
        for name,value in transforms.iteritems():
            new_value = False # previously-performed transform with new value?
            if name in self._performed_transforms:
                if value == self._transforms[name]:
                    self._log.debug("Transform [{0}={1}] already done"
                                    .format(name,value))
                    continue
                else:
                    self._log.info("Will perform transform '{0}' with new "
                                    " value '{1}'".format(name,value))
                    new_value = True
            if name == "pm_relative_longitude":
                self._perform_pm_rel_lon_xform()
            elif name == "make_2d":
                self._perform_2d_latlon_xform()
            else:
                raise Exception("Unknown transform '{0}'".format(name))
        
        
    def _perform_2d_latlon_xform(self):
        """ Convert lat and lon to 2-d using meshgrid """
        assert self._lons.ndim == self._lats.ndim
        if self._lons.ndim == 1:
            self._log.debug("Making lat and lon 2-d")
            #import pdb ; pdb.set_trace()
            (self._lons, self._lats) = np.meshgrid(self._lons,self._lats)
            if "t" in self.lat_dim_order or "t" in self.lon_dim_order:
                assert self.lat_dim_order.index("t") == 0
                assert self.lon_dim_order.index("t") == 0
                self._log.info("Changing dimension order for lat/lon to 'tij' due to make2d transform")
                self.lat_dim_order = "tji"
                self.lon_dim_order = "tji"
            else:
                self._log.info("Changing dimension order for lat/lon to 'ij' due to make2d transform")
                self.lat_dim_order = "ji"
                self.lon_dim_order = "ji"
                
    def _perform_pm_rel_lon_xform(self):
        """
        Perform the "Prime Meridian relative longitude" transform:
        (1) Ensure that the values of self._lons is a real number relative 
            to the prime meridian (e.g. -180W to 180W instead of 0 to 360)
        (2) Ensure longitude values go from westernmost point to easternmost 
            point of the domain by shifting the data. If this is done, 
            update self.specdata.field_transforms() with the mapping 
        {"horz_shift":<shift amount>}
        So that fields will know they need to be shifted accordingly to 
        synchrhoize with the longitude coordinate array.
        """
        #import pdb ; pdb.set_trace()
        # Change to -180->180 if necessary
        if np.max(self._lons) > 180.0 or np.min(self._lons) < -180.0:
            self._log.info("Changing lon to be relative to prime meridian")
            self._lons += 180.
            self._lons %= 360.
            self._lons -= 180.
        
        # Shift grid to be west->east
        diffs = np.diff(self._lons)
        #if diffs.ndim == 2:
        #    diffs = diffs[0]
        if np.all(diffs > 0):
            self._log.debug("Shifting is not necessary")
        ##if self._lons.min()  == self._lons[0]:
        ##     self._log.debug("Shifting is not necessary")
        else:
            if "t" in self.lon_dim_order:
                # The best solution to this is just squeeze out the time dim,
                # since the HorizontalCoordinate doesn't need it
                raise Exception("Not implemented: Cannot have time dimension "
                                "in self._lons to perform pm_rel xform")
            # Not sure what exactly must be done with non latlon grids. At the 
            # very least, the lat grid must also be shifted, which is not 
            # currently done.
            assert self.__is_latlon_grid 
            
            # If lon is 2d and 'i' comes before 'j', temporarily tranpose
            transposed = False
            if "j" in self.lon_dim_order and self.lon_dim_order.index("j")!= 0:
                self._lons = self._lons.T
                transposed = True
                
            # want to shift along the "i" axis, which is the inner axis
            # if we did the 2d meshgrid 
            #import pdb ; pdb.set_trace()
            idx = 1 if self._lons.ndim == 2 else 0
            # shift point is wherever diff is below 0 
            # e.g. 100, 150, -150, -100) -> [50,-300,50]
            #shift = np.where(diffs < 0.)[idx] # diff -> same-dim-array as lons
            shift = np.where(diffs[idx] < 0.)[0] # diff -> same-dim-array as lons
            ##shift = np.where(self._lons == self._lons.min())
            if self._lons.ndim == 2:
                # there should only be one horizontal point where diff<0
                assert np.all(shift == shift[0]) 
            shift += 1 # since it's the diff'd array
            self._transforms["horz_roll"] = shift
            
            self._lons,performed = shift_data(self._transforms, self._performed_transforms,
                                              self._lons, shiftAxis=idx, log=self._log)
            self.specdata.field_transforms.update(performed)
            self._performed_transforms.update(performed)
            
            if transposed:
                self._lons = self._lons.T
                transposed = False
                
    @property 
    def latlons(self):
        """
        Return data with whatever the currently set transforms are. Use the 
        set_transform() method to add transforms to be performed. If any of 
        the transforms require corresponding tranforms to fields, call 
        the update_field_transforms() of self.specdata accordingly.
        """
        self._perform_transforms()
        #self.specdata.update_field_transforms(self._field_transforms)
        # -> updating directly in method calls
        return (self._lats, self._lons)
    
class VerticalLevel(object):
    pass 

class ArbitraryVerticalLevel(VerticalLevel):
    def __init__(self):
        pass

class SpecifiedVerticalLevel(VerticalLevel):
    """ 
    Vertical level with a specified type (e.g. isobaric,
    hybrid, surface-delta, etc.
    Class fields:
        data - Array containing the level data
        dim_name - Name of the dimension. If there is a 
        corresponding Nio file, this should be the same
        as the dimension name in the file.
    """
    def __init__(self, dimName, staggerLoc, levData=None, fileName=None,
                 levVar=None, inDataUnits=None, inDataOrder=None):
        """
        :param dimName: Name of the dimension 
        :type dimName: str
        :param staggerLoc: Stagger location (should be one of the ESMF types)
        :type staggerLoc: str
        :param levVals: Array with level values. If None, use nioFile and levVar
        :type levVals: numpy.ndarray
        :param fileName: Name of file containing the level data. This should 
                         __NOT__ be used to access the data directly by outside
                         programs, since the data may be transformed. Use the 
                         `data' member instead. The only reason this is passed in
                         is so that data is not read unnecessarily before it is
                         needed.
                         This only needs to be passed in if `levData' is not.
        :type fileName: str 
        :param levVar: Name of the variable with level values. It should be
                       specified unless it is derived, since other methods 
                       of SpecifiedForecastDataset use it. It is only used here
                       if fileName is passed in (and levData is not)
        :type levVar: str
        """
        self.stagger_loc = staggerLoc
        self.dim_name = dimName
        self.var_name = levVar
        if levData is not None:
            self.__data = levData
            self.__data_units = inDataUnits
            self.__data_order = inDataOrder
        elif fileName is not None:
            assert levVar is not None
            self.__data_file = fileName
            self._lev_var = levVar
            # self.__data and self.__data_units will be populated 
            # on first call to self._data
        else:
            raise Exception("Must specify either levData OR fileName and levVar")
        self._performed_transforms = {}
        self.__data_order = inDataOrder
        self.transforms = None
        
    @property
    def _data(self):
        if self.__data is None:
            var = nioFile.variables[self.__data_file]
            self.__data = var[:]
            self.__data_units = getattr(var, "units")
            if self.__data_order is None:
                self.get_data_order() # TODO
        if self.transforms is not None:
            raise Exception("Not implemented")
        return self.__data
    
    def get_nearest(self, wanted, wantedUnits):
        """
        Find the closest value in self.data to the given value
        :param wanted: The value we want to get closest or equal to.
        Note: You probably don't want this for hybrid levels.
        :type wanted: Number
        :return -- a 2-tuple of (nearest value, its index)
        """
        data = self._data
        idx = (np.abs(data - wanted)).argmin()
        return (data[idx], idx)

    def get_index(self, wanted, wantedUnits):
        """
        Find index of value in self.data
        """
        data = self._data
        # Convert lev data to given units if necessary
        if self.__data_units.replace(" ","") != wantedUnits.replace(" ",""):
            data = Units.conform(data, Units(self.__data_units), 
                                 Units(wantedUnits))
        idx = np.where(data == wanted)
        assert len(idx) == 1
        return idx[0]
    
    def get_data(units, transforms=None):
        if self._data is None:
            fil = get_open_file(self.__data_file)
            var = fil.variables[self._lev_var]
            self._data = var[:]
            self._data_dims = get_dim_order(var, dimnames, self._log)
        if transforms is not None:
            self.transforms.update(transforms)
        if self._data.ndim > 2:
            raise Exception("Not implemented - gotta shift data")
            #performed = shift_data(self.transforms, self._performed_transforms,
                                   #self._data)
        return self._data
   
    def set_data(self):
        """
        Populates self._data. This is needed if there is no variable with the
        data. For example: In WRF, the 2-, 10-, 30-meter winds are all separate
        2-d variables. But we represent them as a 3-d "sfcDelta" field, so we
        use this routine to create the data
        """
        # TODO

class IsobaricLevel(SpecifiedVerticalLevel):
    
    def interpolate(target_lev):
        raise Exception("Not implemented - use my press2presss")
   
    def get_2d_slice(wantedPres, wantedUnits, srcData, srcDataDimOrder,
                     sliceRetrievalSetting ):
        """
        Get a 2-D slice of 3-D array `srcData'.
        :param wantedPres: The desired pressure
        :param wantedUnits: The units of the given desired pressure
        :param srcData: The 3-d array of data being sliced.
        :param srcDataDimOrder: 3-character string specifying the dimension 
               order ('k' for level, 'j' for lat, 'i' for lon)
        """
        
        if self.__data_units != wantedUnits:
            wantedPres = Units.conform(wantedUnits, Units(wantedUnits),
                                       Units(self.__data_units))
            
        if sliceRetrievalSetting == 'exact':
            try:
                idx = self.data.index(wantedPres)
            except:
                self._log.critical("Selected `sliceRetrievalSetting' of 'exact'"
                                    " but level {wLev} is not in level data"
                                    .format(wLev=wantedPres))
                sys.exit(44)

        elif sliceRetrievalSetting == 'nearest':
            _,idx = self.get_nearest(wantedPres)
        elif sliceRetrievalSetting == 'interpolate':
            raise Exception("Not implemented") # TODO (see generic_plot::vinterp_var)
        else:
            raise Exception("Unknown slab_retrieval_setting: {0}"
                            .format(slabRetrievalSetting))

        if sliceRetrievalSetting in ('exact', 'nearest'):
            assert len(srcDataDimOrder) == 3
            kIdx = srcDataDimOrder.index("k")
            if kIdx == 0:
                return srcData[idx]
            elif kIdx == 2:
                return srcData[:,:,idx]
            else:
                raise Exception("Unexpected dimensionality")

class HybridLevel(SpecifiedVerticalLevel):
    
    def __init__(self):
        """
        Initialize a HybridLevel. Since pressure is 3d for hybrid levels, keep track of whether the data has been "rolled" horizontally via 
        self.transforms["horz_roll"]
        """
        self.transforms = {}
        
        self.transforms["horz_roll"] = 0
        super(HybridLevel, self).__init__()
        
    def interpolate(target_lev):
        raise Exception("Not implemented - use my press2press")
    
class SurfaceDeltaLevel(SpecifiedVerticalLevel):
    def get_field_data_slab(var, levValue, lookupSetting,
                            slabRetrievalSetting):
        """
        :param var: The Nio variable containing the field's data
        :type var: Nio.Variable
        :param levValue: The vertical value we need the slab for
        :type levValue: float
        :param lookupSetting: The level lookup setting ... inputspec or field_data
        :type lookupSetting: str
        :param slabRetrievalSetting: The ... exact, nearest, interpolate
        :type slabRetrievalSetting: str

        """
        # If the Variable is 2-d, there should be an attribute that
        # specifies the level value
        if var.ndim == 2: # TODO : can't rely on this since there may be a time dim - gotta use specData 
            try:
                lev = getattr(var, "level")
                if float(lev) == levValue:
                    return var[:]
                else:
                    raise Exception("Mismatch level: wanted={0}, data={1}".format(levValue, lev))
            except AttributeError:
                self._log.debug("No 'level' attribute in {0}".format(var))
            if hasattr(self, "level"):
                lev = getattr(self, "level")
            else:
                raise Exception("IDK")
            if float(lev) == levValue:
                return var[:]
            else:
                raise Exception("Mismatch level: wanted={0}, data={1}".format(levValue, lev))
        elif var.ndim == 3: # TODO : see above
            if slabRetrievalSetting == "exact":
                idx = self._lev_var.index(levValue)
                # use PyNio's extended selection indexing so we don't need to know order
                #fmt = "{levVar}|i2 {lonVar}|i0:3 {latVar}|i0:3" # getting full slab so no need to pass in lat/lon
                #fmt = "{levVar}|i{idx} "
                #      .format(levVar=self._lev_dim_name, idx=idx)  #TODO lev_dim_name
                    #return var[fmt]
                #Actually, this won't work since we need to use the conformed data
                #specdata needs to have a "native_dim_order" and a could also have a 
                #get_slab(idx) that gets a vertical index

                # Also TODO: Need a data() routine/property for conformedField that
                #  gets the 3d data, which means we need a separate routine for conforming data (2d/3d/4d)

    def interpolate(target_delta):
        raise Exception("Not implemented")
    
#
# COMMON FUNCTIONS
#
def shift_data(needed_transforms, performed_transforms, data, 
               shiftAxis=0, log=None):
    """
    shift `data' horizontally if needed_transforms["horz_roll"] is 
    present and does not equal performed_transforms["horz_roll"]
    This function updates performed_transforms and data
    :return 2-tuple The shifted data and dict of performed transforms
    """
    if log is None: log=_default_log()
    shift_needed = 0
    shift_performed = 0
    try:
        shift_needed = needed_transforms["horz_roll"]
        shift_performed = performed_transforms["horz_roll"]
        shift_needed -= shift_performed
    except KeyError:
        pass
    if shift_needed != 0:
        data = np.roll(data, shift_needed, axis=shiftAxis)
        log.info("Shifting data values @ axis={0} by {1}"
                 .format(shiftAxis, shift_needed))
    return data,{"horz_roll": shift_needed + shift_performed}

def get_dim_order(nioVar, known_dims, log):
    """
    Get a string representing the dimension order of a Nio variable
    :type nioVar: Nio variable
    :param known_dims: A dictionary with 4 keys:
                        "time", "lats", "lons", "levs"
              Each key maps to another dictionary that maps dimension variable 
              names to their corresponding dimension names. Only the dimension
              name is used here (not the key)
    :return: -- A string consisting of the letters t (time), k (vertical),
                j (latitude), i (longitude), or a subset thereof
    """
    #latDimName = None
    #lonDimName = None
    #levDimName = None
    #timeDimName = None
    #dim_names = [timeDimName, levDimName, latDimName,lonDimName]
    timeIdx,levIdx,latIdx,lonIdx = range(4)
    dim_names = [None for i in range(4)]

    #import pdb ; pdb.set_trace()
    #{'levs': {'PRESSURE': 'PINT'}, 'lons': {'HLON': 'west_east', 'VLON': 'west_east'}, 'lats': {'VLAT': 'south_north', 'HLAT': 'south_north'}, 'time': {}}
    var_dims = nioVar.dimensions
    # Determine dimension names.
    # If dimension starts with "lat", "lats", "lev", or "lv_", 
    # assume accordingly
    for dim in var_dims:
        if "lat" in dim: 
            dim_names[latIdx] = dim
            # TODO : Do this verification for others as well
            try:
                expected = "degrees_north"
                units = getattr(nioVar.file.variables[dim], "units") 
                if units != expected:
                    log.warn("Unexpected units for assumed lat var: {0}={1}"
                             .format(dim, units))
            except AttributeError:
                log.warn("Unable to verify var {0} is in units of {1}"
                         .format(dim, expected))
        elif "lon" in  dim: dim_names[lonIdx] = dim # should be degrees_east
        elif "lev" in dim or "lv_" in dim: dim_names[levIdx] = dim
        elif dim.lower().startswith("time"): dim_names[timeIdx] = dim
    # Otherwise use known_dims
    #for i,currDimName,dimKey in zip(range(len(dim_names)), dim_names, ["time", "lats", "lons", "levs"]):
    for i,currDimName,dimKey in zip(range(len(dim_names)), dim_names, ["time", "levs", "lats", "lons"]):
        if currDimName is not None: continue
        for varname,dimname in known_dims[dimKey].iteritems():
            if dimname in var_dims: 
                currDimName  = dimname
                dim_names[i] = dimname 
        #if currDimName is None:
        #    log.debug("Unable to determine the '{0}' dimension from "
        #                    "variable's dimensions: '{1}' based on known "
        #                    "dimensions: '{2}'"
        #                    .format(dimKey, var_dims, known_dims)) 
    # Ensure all of Variable's dimensions' have been determined
    if len(var_dims) != len([d for d in dim_names if d is not None]):
        raise Exception("Unable to determine one or more dimension types. "
                        "Variale dims: {0} ; guessed dimension names: {1}"
                        .format(var_dims, dim_names))

    # Create the string
    retstr = "0123"
    idc = (timeIdx, levIdx, latIdx, lonIdx) = (None, None, None, None)
    dim_strings = ["t", "k", "j", "i"]
    #dim_strings = ["t", "j", "i", "k"]
    #import pdb; pdb.set_trace()
    for idx,dimName, dimStr in zip(idc, dim_names, dim_strings):
        if dimName is not None:
            try:
                where = var_dims.index(dimName)
                retstr = retstr.replace(str(where), dimStr)
            except ValueError:
                log.debug("No {dimType} dimension in variable {var}"
                          .format(dimType=dimName, var=nioVar.varname)) 
    
    # strip out dims that are not present
    for n in (0,1,2,3):
        retstr = retstr.replace(str(n), '')
    
    #import pdb ; pdb.set_trace()
    if len(retstr) != len(var_dims):
        raise Exception("Unable to determine dimension string for variable '{0}'"
                        "with dimensions '{1}' and known_dimensions='{2}'"
                        .format(nioVar.varname, var_dims, known_dims)) 
    return retstr

def ___get_dim_order(nioVar, dimMappings, log):
    """
    Get a string representing the dimension order of a Nio variable
    :type nioVar: Nio variable
    :param dimMappings: Dictionary containing vertical configuration settings
                    for the variable, including: 
                        timeDimName - Name of the time dimension
                        levDimName - Name of the level dimension
                        latDimName - Name of the latitude dimension
                        lonDimName - Name of the longitude dimension
                           
    :return: -- A string consisting of the letters t (time), k (vertical),
                j (latitude), i (longitude)
    """
    dims = nioVar.dimensions
    retstr = "0123"
    idc = (timeIdx, levIdx, latIdx, lonIdx) = (None, None, None, None)
    dim_names = ["timeDimName", "levDimName", "latDimName", "lonDimName"]
    dim_strings = ["t", "k", "j", "i"]
    for idx,dimName, dimStr in zip(idc, dim_names, dim_strings):
        try:
            where = dims.index(dimMappings[dimName])
            retstr = retstr.replace(str(where), dimStr)
        except ValueError:
            log.debug("No {dimType} dimension in variable {var}"
                      .format(dimType=dimName, var=nioVar))
    # strip out dims that are not present
    for n in (0,1,2,3):
        retstr = retstr.replace(str(n), '')
    return retstr

def get_index_for_time(dataset, timeDim, timeWanted):
    """
    Get the index in the time dimension of `dataset' corresponding to
    `timeWanted'
    """
    if timeWanted == 0:
        return 0
    else:
        raise Exception("Not implemented")
 
def get_units_from_grib_defaults(lev_type):
    """
    Determine the units of a grib field's vertical level based on other
    metadata. This is useful for fields that do not have a corresponding
    level dimension because there is only one such GriB message in the 
    file, so Nio represents it as a 2-d field with an extra 'level' attribute.
    :param lev_type: Either a string representing the "level_type" (Grib2) 
                     OR a an int representing a "level_indicator" (GriB 1)
    TODO: Read from GriB tables instead of assuming. This is particularly 
          important from GriB 1 since we just use the level_indicator. For GriB 2,
          the lev_type stringcontains the units (e.g. 'Isobaric surface (Pa)')
    """
    pres_diff_str =  "Level at specified pressure difference from ground"\
                     "to level (Pa)"

    if lev_type == 100:
        return "hPa"
    elif lev_type == "Isobaric surface (Pa)" or lev_type == pres_diff_str:
        return "Pa"
    elif lev_type == 112:
        return "cm"
    elif lev_type == "Depth below land surface (m)" \
      or lev_type == "Specified height level above ground (m)" \
      or lev_type == "Specific altitude above mean sea level (m)" \
      or lev_type == 105:
        return "m"
    else:
        raise Exception("Unable to determine units")

def get_lev_dim_info_from_nio_file(varName, nioFile):
    """
    Get dimension information for a given variable. This should work for GriB
    files and for other data files whose dimension variable is contained in the
    same file and whose name matches a certain paradigm.
    :return (dim_name, var_name): 2-tuple consisting of the dimension name and
                the variable name(s). The latter will be a string if only one 
                corresponding variable or a 2-element arrray if two corresponding
                variables (e.g. for value-between-2-layers level types)
    """
    if isinstance(nioFile, str):
        nio_file = get_open_file(nioFile)
        must_release = True
    else:
        nio_file = nioFile
        must_release = False
    lev_type = get_lev_type_from_attr(varName, nio_file)
    var = nio_file.variables[varName]
    vardims = var.dimensions
    # Get dim name
    try:
        dim_idx = [dimName[:3] for dimName in vardims].index("lv_")
    except ValueError:
        log.critical("Unable to determine vertical dimension name. Field "
                    "dimensions: '{0}'. If this field is being represented "
                    " as 2-d because there is only one vertical dimension,"
                    " this is a bug as the code should not have come here"
                    .format(vardims))
        sys.exit(46)
    dim_name = vardims[dim_idx]
    
    # Get var name
    if lev_type in ("isobaric", "sfcDelta", "mslDelta"):
        var_name = dim_name
    elif lev_type == "soilDelta":
        # lev-between-2-soil-depths has two dimension variables - one 
        # for the starting point of the layer and on for the ending point
        var_name = [var_name+"_l1", var_name+"_l2"]
    else:
        raise Exception("Uknown level type")
    if must_release:
        close_open_file(nioFile)
    return (dim_name, var_name)
        
def get_native_name_from_var_attributes(attrs, nioFileName, log=None, fmt=None):
    """
    Given a string of attributes, find the variable in the Nio file that matches
    the attributes. The attributes list must consist of the following format:
        attr(foo=bar&&bar=baz || me=myself  )
     || is for "or" a. Each substring preceeding one of these is an "attribute set"
     && is for "and"
     The algorithm will go from left to right and match the first 
     attribute set with exactly one matching field in the file.
     TODO : When it does not find a matching field (or if there are multiple identical
            fields in the same grib file, it ultimately crashes with an error saying
            that the last file in the inspec coud not be found
    """
    #log = _default_log() if log is None else log
    #import pdb ; pdb.set_trace()
    if log is None:
        log = _default_log()
    nio_file = get_open_file(nioFileName, fmt=fmt)
    assert attrs.startswith("attr(") and attrs.endswith(")")
    attrs = attrs[5:-1]
    # Try to match any ...
    possible_attr_sets = attrs.split("||")
    for attrSet in possible_attr_sets:
        attrList = attrSet.split("&&")
        matches = copy.copy(nio_file.variables) # dict; reading only, so shallow copy
        # Try to match all (by recursively removing non-matching vars)
        for attrRequirements in attrList:
            non_matching_keys = []
            k,v = [ x.strip() for x in attrRequirements.split("=")]
            #matches = [{nioVarName:nioVarObj} for nioVarName,nioVarObj \
            #           in matches.iteritems() if getattr(nioVarObj, k) == v ]
            for nioVarName,nioVarObj in matches.iteritems():
                #if "MSL" in nioVarName: import pdb ; pdb.set_trace()
                # NOTE : All attributes are converted to string
                if not hasattr(nioVarObj, k) or not str(getattr(nioVarObj, k)) == v:
                    #matches.pop(nioVarName)
                    #if nioVarName == "PRMSL_GDS0_MSL": import pdb ; pdb.set_trace()
                    #log.debug("Not a match: {0}:{1} for attrs {2}={3}"
                    #          .format(nioVarName, nioVarObj.__dict__, k, v))
                    non_matching_keys.append(nioVarName)
                #else:
                    #log.debug("Found match: {0}:{1}".format(nioVarName, nioVarObj.__dict__))
            for nmk in non_matching_keys:
                #log.debug("Popping {0}".format(nmk))
                matches.pop(nmk)
                
        #import pdb ; pdb.set_trace()
        if len(matches) == 1:
            return matches.keys()[0]
        elif len(matches) == 0:
            log.debug("No matches found for '{0}'. Trying next attribute set"
                      .format(attrSet))
        else:
            log.debug("Found more than one matching Variable for attribute set "
                      "'{0}'. Matches: {1}.....Continuing search"
                      .format(attrSet, matches))
    raise UnmatchedFieldException("No variables found in file '{fil}' that "
                                  "match attribute set '{attrs}'"
                    .format(fil=nioFileName, attrs=possible_attr_sets))
    close_open_file(nioFileName)

def get_lev_type_from_attr(fieldName, nioFile):
    """
    Given a field name and file, determine it's level type by opening it
    with PyNIO and reading it's level_type or level_indicator attribute
    :return type: String representing the type of level:
        isobaric, sfcDelta, mslDelta, soilDelta
    """
    if isinstance(nioFile, str):
        nio_file = get_open_file(nioFile)
        must_release = True
    else:
        nio_file = nioFile
        must_release = False
    var = nio_file.variables[fieldName]
    if hasattr(var, "level_type"):
        # Grib 2
        lev_type = getattr(var, "level_type")
    elif hasattr(var, "level_indicator"):
        # Grib 1
        lev_type = int(getattr(var, "level_indicator"))
    else:
        raise UnknownLevelTypeException()
    if must_release:
        close_open_file(nio_file)
    # determine string based on self-described file's attributes
    if lev_type == 100 or ("isobaric surface" in str(lev_type).lower()):
        return "isobaric"
    elif (lev_type == 112) or ("Depth below land surface" in str(lev_type)):
        return "soilDelta"
    elif (lev_type == 103) or ("height level above ground" in str(lev_type)):
        return "sfcDelta"
    elif (lev_type == 102) or ("altitude above mean sea level" in str(lev_type)):
        return "mslDelta"
    elif "Mean sea level (Pa)" in str(lev_type) or "Ground or water surface" in str(lev_type) \
         or lev_type in (1, 102):
        # TODO ? Currently set MSLP as 2d, should it be mslDelta?
        return "2d"
    # TODO : others http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-5.shtml
    #elif (lev_type 
    else:
        #import pdb ; pdb.set_trace()
        raise UnknownLevelTypeException(lev_type)

def get_2d_staggerloc_from_3d_staggerloc(stagger3d):
    """
    Determine stagger location of a horizontal slice of 3D data 
    at a given stagger location
    """
    from ESMF.api.constants import StaggerLoc as sl
    #mappings = { "CENTER" : ["CENTER_VCENTER", "CENTER_VFACE"],
    #             "EDGE1" : ["EDGE1_VCENTER", "EDGE1_VFACE"],
    #             "EDGE2" : ["EDGE2_VFACE", "EDGE2_VCENTER"],
    #             "CORNER" : ["CORNER_VFACE", "CORNER_VCENTER"] }
    mappings = {sl.CENTER_VCENTER:sl.CENTER, sl.CENTER_VFACE:sl.CENTER,
                sl.EDGE1_VCENTER:sl.EDGE1, sl.EDGE1_VFACE:sl.EDGE1,
                sl.EDGE2_VCENTER:sl.EDGE2, sl.EDGE2_VFACE:sl.EDGE2,
                sl.CORNER_VFACE:sl.CORNER, sl.CORNER_VCENTER:sl.CORNER }
    return mappings[stagger3d]

class UnknownLevelTypeException(Exception):
    pass
class UnmatchedFieldException(Exception):
    pass

_logger = None
def _default_log(log2stdout=logging.INFO, log2file=None, name=None):
    # Need this or when calling _default_log from multiple classes or else
    # if I don't pass it to constructor, I end up with multiple loggers
    #global _logger
    #if _logger is not None:
    #    return _logger
    if name is None:
        name = inspect.stack()[1][3]
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

def get_open_file(fileName, fmt=None):
    """
    Get a Nio OpenFile object for the given path and increment the
    global counter to indicate another task is using it
    """
    global _g_open_files 
    #import pdb ; pdb.set_trace()
    if fileName in _g_open_files.keys():
        _g_open_file_count[fileName] += 1
        return _g_open_files[fileName]
    else:
        _g_open_file_count[fileName] = 1
        if fmt is None:
            ext = os.path.splitext(fileName)[1]
            if not ext in ["nc", "nc4", "grib", "grib2", "grb", "grb2", "nc3"]:
                if "grb" in fileName or "grib" in fileName:
                    fil =  Nio.open_file(fileName, format="grib")
                else:
                    fil = Nio.open_file(fileName)
            else:
                fil =  Nio.open_file(fileName)
        else:
            fil =  Nio.open_file(fileName, format=fmt)
        _g_open_files[fileName] = fil
        return _g_open_files[fileName]

def close_open_file(key):
    """
    Close an open file if there are no tasks using it
    @param key Either the file name (as string) or the Nio.NioFile object
    """
    global _g_open_files 
    if isinstance(key, str):
        fileName = key
    #elif isinstance(key, Nio.NioFile):
    else: # Since NioFile is a proxy object, can't use isinstance
        el_id = id(key)
        filenames = _g_open_files.keys()
        # get index of file name by matching ID of OpenFile
        #idx = [id(v) for k,v in filenames.iteritems()].index(el_id)
        idx = [id(v) for v in _g_open_files.itervalues()].index(el_id)
        fileName = filenames[idx]
    #else:
    #    raise Exception("Invalid key")
    _g_open_file_count[fileName] -= 1
    if _g_open_file_count[fileName] == 0:
        _g_open_files[fileName].close()
        _g_open_files.pop(fileName)

"""
Global dictionaries that keep opened files and the number of times they've been
requested
"""
_g_open_files = collections.OrderedDict()
_g_open_file_count = collections.OrderedDict()

if __name__ == "__main__":
    inspec = ["/home/Javier.Delgado/apps/pycane_dist/trunk/scripts/map/generic_plotter/inputspec_g5nr_hires.conf"]
    topdir = "/home/Javier.Delgado/scratch/nems/g5nr/data/gamma/combined_nc"
    # test after: topdir = "/home/Javier.Delgado/scratch/nems/g5nr/data/raw_collections"
    date = dtime(year=2006, month=9, day=6) 
    log = _default_log(log2stdout=logging.DEBUG)
    grib = SpecifiedForecastDataset(inspec, topdir, date, log=log)
    el_var = grib.get_field("air_temperature", "K", "isobaric")
    #[slab,coords,speclev] = 
    fslice = el_var.get_horz_slice(levIdx=2, outDimOrder=None, sliceRetrievalSetting="exact", fcstOffset=None)
                           # levValue=None, levUnits=None,
    slab = fslice.data
    coords = fslice.coords

    # we use get_horz_slice for 2d fields too, since it also conforms the data and gets the time slice
    mslp = grib.get_field("mslp", "hPa")
    #(slab,coords,_) = mslp.get_horz_slice()
    fslice = mslp.get_horz_slice()
    slab = fslice.data
    print slab.min(), slab.max()

    # Needed Tests:
    # retrieving basic data
    # conform/transform data
    # Multiple grib files
    # Multiple netcdf files
    # interpolating
    
    #http://code.activestate.com/recipes/496741-object-proxying/
    #http://stackoverflow.com/questions/17478284/python-is-there-a-way-to-implement-getitem-for-multidimension-array

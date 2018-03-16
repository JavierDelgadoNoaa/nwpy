"""
This module implements a factory object that loads/pre-loads large datasets
from the file system.
It is used in SpecData to mask long load times experienced with the
large datasets (e.g. > 15GB) it was designed to handle.
"""

import os
import sys
import threading
import logging
import time
from datetime import datetime as dtime
from datetime import timedelta as tdelta
from collections import defaultdict

from PyNIO import Nio

##
# CONSTANTS
##
# How long to wait when threshold is reached (seconds)
POLL_INTERVAL = 2

##
# MODULE FUNCTIONS
##
def _default_log(log2stdout=logging.INFO, log2file=None, name='el_default'):
    global _logger
    if _logger is None:
        _logger = logging.getLogger(name)
        _logger.setLevel(log2stdout)
        msg_str = '%(asctime)s::%(name)s::%(lineno)s::%(levelname)s - %(message)s'
        #msg_str = '%(asctime)s::%(funcName)s::%(filename)s:%(lineno)s::%(levelname)s - %(message)s'
        # TODO : Use this one if running in debug mode
        #msg_str = '%(levelname)s::%(asctime)s::%(funcName)s::%(filename)s:%(lineno)s - %(message)s'
        msg_str = '%(levelname)s :: %(filename)s:%(lineno)s - %(message)s'
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

def timeit(method):
    """
    Wrapper to be used to decorate methods with timing information
    """
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        elapsed_sec = (te - ts) * 1000.
        if "timingLogger" in kw:
           logger = kw["timingLogger"]
           logger.info("Timing for {0}: {1}".format(method.__name__, elapsed_sec))
        return result
    return timed

##
# CLASSES
##
class DatafileFactory(threading.Thread):
    """
    This class implements a factory for providing datafiles.
    It was designed to mask the potentially long load times of SpecifiedForecastDataset
    file loads by preloading sequences of files (assuming they are loaded 
    chronologically by the calling program).
    The factory runs as a thread and it creates threads to load the data files
    in parallel.
    Example:
        factory = DatafileFactory(sdate, sDelta, durationDelta, intervalDelta, 
                                  7, 3, {"domain":1})
        # Specify filegroups to load
        factory.add_filegroup("U10_dom{domain}_f{fhr:3d}.grb2")
        factory.add_filegroup("T10_dom{domain}_f{fhr:3d}.grb2")
        factory.close_doors() # let it know it does not need to load any additional
                              #  datasets - so when all files from all forecast 
                              # offsets for U10 and V10 are loaded, it can terminate
        # Start thread that loads all files from U10 and T10 datasets
        factory.start()
        for pattern in patterns:
            for foffset in foffsets: # list of timdeltas
                fname = fil_pattern.date_interp(foffset)
                # Wait for file to load and then use it
                dataset = factory.wait_for_file(fname, logTiming=True)
                do_stuff_with_dataset
    """
    def __init__(self, initDate, fcstOffset, duration, interval, maxPreloadTotal,
                 maxPreloadPerDataset, group=None, target=None, name=None, 
                 verbose=None, xtraInterpArgs={}, log=None):
        """
        Initialize member variables and Set up all the date-related parameters.
        The date-related parameters will determine the sequence of files
        to be loaded for each filegroup added.
        :param initDate: Initialization date
        :type initDate: datetime.datetime
        :param fcstOffset: The initial forecast offset to start analyzing from 
        (difference from analysis time)
        :type fcstOffset: datetime.timedelta
        :param duration: Duration fo the the forecast
        :type duration: datetime.timedelta
        :param maxPreloadTotal: Max number of data files to load at once
        :type maxPreloadTotal: int
        :param maxPreloadPerDataset: Max number of data files to load for
         each dataset. This is useful since multiple file may correspond
         to a given time
        :param: xtraInterpArgs : Dictionary containing additional parameters
         that should be interpolated for each file name
        """
        super(DatafileFactory,self).__init__(group=group, target=target, 
                                      name=name, verbose=verbose)
        ## Initialize member variables
        # Initialize parameters passed to constructor
        self.init_date = initDate
        self._xtra_fname_interp_args = xtraInterpArgs 
        self._duration = duration
        self._first_fcst_offset = fcstOffset
        self._interval = interval
        self.max_preload_total = maxPreloadTotal
        self.max_preload_per_dataset = maxPreloadPerDataset
        if log is None:
            log = _default_log()
        self._log = log
        # Bool to let the factory know if it is done
        self._closed_doors = False
        # List of dataset file patterns to process
        self._datafile_patterns = []
        # Dictionary that maps a file pattern (i.e. pre-date-interpolated file name)
        # to the list of files it consists of. 
        self.filegroups = {}
        # List of files that have been loaded/opened. Each item will be a dict
        # with keys "path" (path to file) and "dataset" (Nio Dataset object)
        self.loaded_files = []
        # List of forecast offsets to process
        self.fcst_offsets = self.get_fcst_offset_list()
        # Dict to keep track of #files being loaded for each filegroup
        self._num_loading = {}
    
    @property 
    def _num_loading_total(self):
        return sum(self._num_loading.values())

    def get_fcst_offset_list(self):
        """ Get list of forecast offsets that will be processed.
        The range returned will be from initial date plust initial offset 
        specified in constructor through the duration, with the interval
        specified in the constructor."""
        # XXXIt's in reverse order since
        # XXX datasets to be loaded will be popped from a list in chronologic order
        #offset_secs = range(duration.total_seconds(), -1, interval)
        duration_sec = int(self._duration.total_seconds())
        interval_sec = int(self._interval.total_seconds())
        first_fcst_offset_sec = int(self._first_fcst_offset.total_seconds())
        offset_secs = range(first_fcst_offset_sec, duration_sec + 1, interval_sec)
        offsets = [ tdelta(seconds=s) for s in offset_secs ]
        #date_list = [ (initDate + fcstOffset + tdelta(seconds=t)) \
        #              for t in self.fcst_offsets) ]
        return offsets

    def add_filegroup(self, pattern, fileFormat=None):
        """
        Add a datafile pattern to process for all forecst offsets. 
        (i.e. add a generic dataset file name to preload)
        This will update self.datasets with the pattern and point it to 
        its corresponding list of files.
        :param pattern: File pattern. i.e. file path before being interpolated.
                        Includes the entire path.
        :param fileFormat: The file format (e.g. "grib", "netcdf"). Only needed
                           if it's not included in file name. 
                           See DatafileLoader::_open_file
        """
        if self._closed_doors:
            raise Exception("close_doors() called; Unable to accept more patterns")
        #self._dataset_patterns.append(pattern)
        self.filegroups[pattern] = []
        self._num_loading[pattern] = 0
        for offset in self.fcst_offsets:
            filename = self.date_interp(pattern, offset)
            loader = DatafileLoader(filename, pattern, 
                                    self.notify_of_loaded_file, log=self._log)
            self.filegroups[pattern].append(loader)

    def date_interp(self, s, fcstOffset):
       """
       Interpolate file name based on current date/forecast offset and any 
       parameters set via self._xtra_fname_interp_args
       """
       curr_date = self.init_date + fcstOffset 
       date_args = {
                    "cal_date":curr_date,
                    "fdate":curr_date,
                    "init_date":self.init_date, 
                    "adate":self.init_date,
                    "fhr": int(fcstOffset.total_seconds() / 3600.), 
                    "fmin": int(fcstOffset.total_seconds() % 3600. / 60.),
                    "fcst_offset_mins": int(fcstOffset.total_seconds() / 60.),
                    }
       date_args.update(self._xtra_fname_interp_args)
       return s.format(**date_args)

    def run(self):
        """
        Run thread that spawns DatafileLoader threads to load files
        Threads will be created as long as the constraints of
        ``max_preload_per_dataset'' and ``max_preload_total'' are satisfied.
        This thread will run until the caller lets it know there are no more
        filegroups via call to close_dors(). Once close_doors() is 
        called, no new datasets may be added, but it will continue loading
        all files pertaining to existing datasets. It will then exit when 
        all files are either loading or loaded.
        Files are loaded in round-robin order. i.e. hour 0 from filegroup 1, 
        hour 0 from filegroups 2, hour 1 from filegroup 1, etc.
        """
        filegroup_iters = {} # iterator for each filegroup
        all_files_loading = defaultdict(lambda: False)
        while(1):
            #import pdb ; pdb.set_trace()
            for filePattern,fileList in self.filegroups.iteritems():
                if all_files_loading[filePattern]: 
                    continue
                if self._num_loading[filePattern] < self.max_preload_per_dataset:
                    # Create iterator if this is the first round
                    if not filePattern in filegroup_iters:
                        filegroup_iters[filePattern] = iter(fileList)
                    # Start thread to load the file
                    try:
                        curr_loader = filegroup_iters[filePattern].next()
                        with threading.Lock():
                            # critical sect since notification callback 
                            # can occur any time
                            self._num_loading[filePattern] += 1
                        curr_loader.start()
                    except StopIteration:
                        #completed[filePattern] = True
                        all_files_loading[filePattern] = True
                        self._log.info("All files of pattern {0} are loading "
                                       "or loaded".format(filePattern))
                    # If exceeding the total max, wait
                    while self._num_loading_total > self.max_preload_total:
                        self._log.info("Reached max number of concurent loads. "
                                       "Sleeping for {0}.".format(POLL_INTERVAL))
                        time.sleep(POLL_INTERVAL)
            # Loop again unless factory is closed _and_ all files loading/loaded
            if self._closed_doors and (not False in all_files_loading.values()):
                break
        # Wait for threads to finish loading files
        self._log.info("Finished looping. Waiting for threads to finish")
        for filePattern,fileList in self.filegroups.iteritems():
        	for loader in fileList:
        		loader.join()
        self._log.info("All files loaded. My work here is done")
        
    def notify_of_loaded_file(self, loader):
        """
        This is called by the DatafileLoader when it finishes loading a file.
        It updates the internal counters for number of loading files
        :param loader: The file loader object that has finished. 
        :type loader: DatafileLoader
        """
        with threading.Lock():
            self.loaded_files.append({"path":loader.path, "dataset":loader.dataset})
            #self._num_loading_total -= 1 # using property 
            self._num_loading[loader.file_pattern] -= 1
    
    @timeit # TODO - giving NoneType retun error
    def wait_for_file(self, filepath, pollInterval=POLL_INTERVAL, timingLogger=None):
        """
        Wait for a given file to load and return.
        This is a BLOCKING call - it will wait until the file is loaded.
        """
        counter = 0
        while(1):
            try:
                idx = [f["path"] for f in self.loaded_files].index(fname)
                return self.loaded_files[idx]["dataset"]
            except ValueError:
                self._log.info("Still waiting for file {0} to load"
                               .format(filepath))
                counter += 1
                if counter == 5:
                    self._log.warn("BLOCKING while waiting for file {0} to lood"
                                   .format(filepath))
                    counter = 0
            time.sleep(pollInterval)
    
    def close_doors(self):
        """
        Let the factory know that no more dataset patterns will be
        accepted and it can close shop when all files from already-added
        patterns are loaded.
        This is needed since the factory doesn't know how many patterns
        it will need to process ahead of time.
        """
        #import pdb ; pdb.set_trace()
        self._closed_doors = True

class DatafileLoader(threading.Thread):
    def __init__(self, path, filePattern, callback, status="queued", format=None, 
                 log=None):
        """
        :param path: Path to dataset file to load
        :param filePattern: The file pattern (i.e. uninterpolated file). 
         THis is needed since it acts as a key in the DatafileFactory.
        :param callback: function to call back when file is loaded.
         The function shoudl take two arguments: The dataset path
         and the dataset pattern
        :param status: The initial status of the datafile.
        """
        super(DatafileLoader, self).__init__()
        self.status = status
        self.path = path
        self._dataset = None
        self.callback = callback
        self.file_pattern = filePattern
        self._file_format=None
        self._log = _default_log() if log is None else log

    def run(self):
        self.status = "loading"
        self._dataset = self._open_file(self.path)
        self._log.info("File {0} has been loaded".format(self.path))
        self.status = "loaded"
        self.callback(self)

    def _open_file(self, path):
        """ Open the file. Since PyNIO cannot always guess the file type, 
            there is additional logic here to guess it or to use 
            self._file_format"""
        kw = {}
        if self._file_format is None:
            ext = os.path.splitext(path)[1]
            if not ext in ["nc", "nc4", "grib", "grib2", "grb", "grb2", "nc3"]:
                if "grb" in path or "grib" in path:
                    kw["format"] = "grib"
        else:
            kw["format"] = self._file_format
        try:
            return Nio.open_file(path, **kw)
        except:
            self._log.error("File {0} could not be opened".format(path))
            sys.exit(5)
        
    @property
    def dataset(self):
        if self._dataset:
            return self._dataset
        else:
            raise Exception("Dataset not loaded. Current status is '%s'"
                            .format(self.status))


## 
# TEST
##
if __name__ == "__main__":
        sdate = dtime(year=2006, month=9, day=4, hour=0, minute=0)
        sDelta = tdelta(seconds=0)
        durationDelta = tdelta(days=1)
        intervalDelta = tdelta(hours=6)
        intervalDelta = tdelta(hours=12)
        maxTotal = 7
        maxPerFilegroup = 3
        factory = DatafileFactory(sdate, sDelta, durationDelta, intervalDelta, 
                                  maxTotal, maxPerFilegroup,
                                  xtraInterpArgs={"domain":1})
        logger = factory._log # bad practice, I know
        # Specify filegroups to load
        #patterns = ["U10_dom{domain}_f{fhr:03d}.grb2", "T10_dom{domain}_f{fhr:03d}.grb2"]
        patterns = ["nmbprs_UGRD_2006090400.f{fhr:03d}.{fmin:02d}.grb2",
                    "nmbprs_VGRD_2006090400.f{fhr:03d}.{fmin:02d}.grb2"]
        topdir = "/home/Javier.Delgado/scratch/nems/g5nr/data/gamma/plan_c/2006090400/nam_physics/postprd"
        patterns = [os.path.join(topdir, p) for p in patterns]
        factory.add_filegroup(patterns[0])
        factory.add_filegroup(patterns[1])
        factory.close_doors() # let it know it does not need to load any additional
                              #  datasets - so when all files from all forecast 
                              # offsets for U10 and V10 are loaded, it can terminate
        # Start thread that loads all files from U10 and T10 datasets
        factory.start()
        foffsets = factory.get_fcst_offset_list()
        for fil_pattern in patterns:
            for foffset in foffsets: # list of timdeltas
                fname = factory.date_interp(fil_pattern, foffset)
                #while(factory.filegroups[pattern].status != "loaded"): time.sleep(20)
                #while(not fname in factory.loaded_files): time.sleep(20)
                #dataset = factory.loaded_files[...]
                # Wait for file to load and then use it
                logger.info("Processing file: {0}".format(fname))
                dataset = factory.wait_for_file(fname, timingLogger=logger)
                #print dataset
                try:
                    v = getattr(dataset.variables["UGRD_P0_L100_GLL0"], "forecast_time")
                except:
                    v = getattr(dataset.variables["VGRD_P0_L100_GLL0"], "forecast_time")
                logger.info("Variables: {0}. ftime: {1}".format(v, dataset.variables.keys()))


'''
Encapsulates the execution of a WRF run.

NOTES
 - When attempting to use on newer (or older) versions of HWRF, pay particular
   attention to areas of code under string "Verify for new HWRF versions"
 - There are multiple ways to determine the end date. The following scheme is used:
     - if the duration is passed into the constructor, use that 
       (actually, since it is by default, this is currently the only value used,
        but logic for the following exists in case we want to change this in the
        future)
        ->This has first priority since the user may be intrested in a certain
          integration time 
     - The timestamp associated with the 'SUCCESS COMPLETE WRF' message gets 
       second priority since that is the official end_date
     - If any of the run_* parameters in the namelist are set, use those
     - Last resort: rely on the end_* parameters in the namelist
CAVEATS
 - The exch_time (i.e. halo_times[]) uses only the d01 halo times, since the
   counter (in HWRF) that produces the RSL files seems to kind of combine the 
   two (or something - I need to figure out what exactly it does.
   NOTE: The rsl file has an entry per forecast hour for exch_time (and other
   parameter), but does not indicate which domain its for
'''
import os
import sys
import logging as log
from datetime import datetime, timedelta
import gzip
import numpy as np
import time
from modelrun import ModelRun

# Exit with non-zero exit code if the "SUCCESS COMPLETE WRF" string is absent
FAIL_IF_UNCONFIRMED_FORECAST = True

class WrfRun(ModelRun):
    ''' 
    Object to encapsulate configuration settings for a WRF run as well as various
    timing statistics.

    The following performance-related variables can be used.

        io_server_selections - List of I/O servers selected, in the order 
                               they were selected
        hist_output_times - Time it took to write out each of the main
                            history files (i.e. wrfout_d*)
        gsi_auxhist_output_times - Time it took to write out the OSSE-specific
                                   aux_hist file with the prefix "gsi_auxhist"
        halo_times - Time it takes to perform halos. This is a dictionary mapping
                     the domain to the number of halos performed. This data is 
                     output once per hour. It corresponds to the "exch_tim" in
                     both HWRF and NMMB
        integrate_times - Time it takes to perform the entire call to integrate(), 
                          for each timestep. This is a dictionary mapping domains 
                          to their times. Each entry corresponds to the integration
                          time for a domain and all its subdomains for a single 
                          timestep. This corresponds to the "timing for main" 
                          data in WRF and the "Finished timestep ..." for NMMB
        filter_times - timing for writing the filter output 
    '''
    def __init__(self, path, duration=timedelta(hours=126)):
        '''
        Instantiate a WrfRun. Does the standard instantiation of the 
        parent class as well as instantiation of parameters specific to WRF.
        See docstring for parent (ModelRun) for important information.
        '''
        super(WrfRun, self).__init__(path, duration)

        # Set variables sepecific to WRF runs 
        
        # Track what I/O servers are used 
        self.io_server_selections = []
        # write time for auxhist files with prefix gsi_auxhist
        self.gsi_auxhist_output_times = {}
         # time to write the filter?
        self.filter_times = {}

        #self._set_defaults()
        self._set_params_from_namelist()

        # Instantiate the dictionaries using the domain number as keys and 
        # an array as values.
        for domNum in range(1, self.num_domains+1):
            self.hist_output_times[domNum] = []
            self.gsi_auxhist_output_times[domNum] = []
            #halo_times.append([])
            self.integrate_times[domNum] = []
            self.filter_times[domNum] = []
            self.halo_times[domNum] = []
        self._set_params_from_rsl_files()


    def _set_start_and_end_dates(self, start_params, end_params, runtime_params):
        '''
        Set the start and end dates as datetime objects, given input parameters
        for start and end date.
        However, if self._desired_duration is less than the calculated duration,
        override the end date accordingly.
        
        INPUT
            * Contains various parameters read from the namelist
            - start_params : Dictionary with keys for 'start_year', 'start_month', 
               'start_day', 'start_hour', (optional) 'start_minute', 
               (optional) 'start_second' 
            - end_params : Dictionary with keys for 'end_year', 'end_month', 
               'end_day', 'end_hour', 'end_minute', (optional) 'end_second' 
            - runtime_params : Dictionary with keys for  
               'run_days', 'run_hours', 'run_minutes', 'run_seconds' .
              All keys in runtime_params are optional. The runtimes will override the
              end_* parameters when determining actual forecast duration, ala
              WRFv3.5 namelist in which the run_* parameters override the end_*.
              This end time, in turn, may be overriden when reading the log file 
              later on, since the log file is the most accurate source for determining
              the actual end time.
        POSTCONDITION
            - The self.start_date and self.end_date parameters will be set
        '''

        # Set default parameters 
        for key in ('minute', 'second'):
            if not start_params.has_key('start_' + key):
                start_params['start_' + key] = 0
            if not end_params.has_key('end_' + key):
                end_params['end_' + key] = 0
        self.start_date = datetime(start_params['start_year'], 
                                   start_params['start_month'],
                                   start_params['start_day'],
                                   start_params['start_hour'],
                                   start_params['start_minute'],
                                   start_params['start_second'])
        # Determine end_date and runtime. If a runtime is specified
        # in namelist, WRF will use it. Otherwise, it will use
        # the specified end_year/month/day/...
        specified_end_date = datetime(end_params['end_year'], 
                                   end_params['end_month'],
                                   end_params['end_day'],
                                   end_params['end_hour'],
                                   end_params['end_minute'],
                                   end_params['end_second'])
        #import pdb ; pdb.set_trace()
        runtime = timedelta(0)
        if runtime_params.has_key('run_seconds'):
            runtime += runtime_params['run_seconds']
        if runtime_params.has_key('run_minutes'):
            runtime += runtime_params['run_minutes'] * 60
        if runtime_params.has_key('run_hours'):
            runtime += runtime_params['run_hours'] * 3600
        if runtime_params.has_key('run_days'):
            runtime += runtime_params['run_days'] * 3600 * 24
        if runtime.total_seconds() == 0:
            runtime = specified_end_date - self.start_date
            #self.end_date = specified_end_date
        overriden_for_user = False
        if runtime > self._desired_duration:
            log.info("Overriding runtime in namelist with the runtime passed " \
                     "to the constructor")
            runtime = self._desired_duration
            overriden_for_user = True
        self.end_date = self.start_date + runtime
        # sanity check/warn
        if (not overriden_for_user) and self.end_date != specified_end_date:
            log.warn("The date specified by the end_* parameters in the namelist"\
                     " does not correspond to the run_* parameters. Since WRF "\
                     " uses the latter, it will be used. Hence, the end_date is %s."\
                     " However, this end_date may be overriden when reading the logfile."
                     %self.end_date)

    def _set_params_from_namelist(self, nml_file="namelist.input"):
        '''
        Set fields that can be read from the namelist
        '''
        
        # Some values' length will depend on max_dom, so they'll be stored in
        # this dictionary temporarily
        pv_list = {}
        # temporarily store start/end/run time params as local dict with keys for year, month,
        # day, and whatever other start_* parameter is defined in the namelist
        start_time = {}
        end_time = {}
        run_time = {}

        with open(os.path.join(self.path,nml_file)) as nl:
            nml_lines = nl.readlines()
        for line in nml_lines:
            toks = [t.strip() for t in line.split("=")]
            if len(toks) < 2: continue
            param = toks[0]
            values = [t.strip() for t in toks[1].strip().split(",")]
            if values[-1] == '': values = values[:-1]
            # Some parameters' values is a list of the number of domains.
            # For those, we add their name:value to pv_list and set them
            # later (when we know max_dom). For parameters
            # whose values are of length one, set the value here.
            if param == 'max_dom':
                self.num_domains = int(values[0])
                if len(values) != 1: raise Exception("Invalid value for max_dom")
            # If parent_ids are unique, it is not multistorm
            elif param == "parent_id":
                parent_ids = [int(x) for x in values]
                if parent_ids == list(set(parent_ids)):
                    self.is_multistorm = False
                else:
                    self.is_multistorm = True
            elif param == 'time_step':
                assert len(values) == 1
                self._dt_sec = int(values[0])
            elif param == 'time_step_fract_num':
                assert len(values) == 1
                self._dt_frac_num = int(values[0])
            elif param == 'time_step_fract_den':
                assert len(values) == 1
                self._dt_frac_den = int(values[0])
            # set the _start_(year|month|day|etc.) parameters
            elif param.startswith('start_'):
                if len(values) > 1 and len(set(values)) != 1:
                    raise Exception("All of the namelist values for %s"\
                                   " should be the same." %param)
                start_time[param] = int(values[0])
            # set the _end_year|month|day|etc. params
            elif param.startswith('end_'):
                if len(values) > 1 and len(set(values)) != 1:
                    raise Exception("All of the namelist values for %s"\
                                   " should be the same." %param)
                end_time[param] = int(values[0])
            elif param.startswith("run_"):
                assert len(values) == 1
                if len(values) > 1 and len(set(values)) != 1:
                    raise Exception("All of the namelist values for %s"\
                                   " should be the same." %param)
                run_time[param] = int(values[0])
            # parameters whose values is a list of size max_dom
            elif param in ('nproc_x', 'nproc_y', 'history_interval',
                            'interval_seconds', "nest_pes_x", "nest_pes_y"):
                pv_list[param] = values

            #nproc_x = [-1 for x in range(self.num_domains)]
            #nproc_y = [-1 for y in range(self.num_domains)]

        for (param,values) in pv_list.iteritems():
            if param == 'nproc_x':
                if self.is_multistorm:
                    log.debug("Ignoring nproc_x for multistorm")
                    if int(values[0]) > -1: 
                        log.warn("Repeat: Ignoring nproc_x for multistorm")
                    continue
                if len(values) > 1:
                    raise Exception("nproc_x array not supported")
                if int(values[0]) > -1:
                    log.debug("nproc_x was passed in namelist")
                    self.nproc_x = [int(values[0])]
            elif param == 'nproc_y':
                if self.is_multistorm:
                    log.debug("Ignoring nproc_y for multistorm")
                    if int(values[0]) > -1: 
                        log.warn("Repeat: Ignoring nproc_x for multistorm")
                    continue
                if len(values) > 1:
                    raise Exception("nproc_y array not supported")
                if int(values[0]) > -1:
                    log.debug("nproc_y was passed in namelist")
                    self.nproc_y = [int(values[0])]
            elif param == 'history_interval':
                # history_interval will be list with index corresponding
                # to domain number will be the value for that domain
                self.history_interval = [None] + [int(v) for v in values]
            elif param == "nest_pes_x" and self.is_multistorm:
                self.nproc_x = values[0:self.num_domains]
                self.nproc_x = [int(x) for x in self.nproc_x]
                log.debug("Setting nproc_x from nest_pes_x = {0}"
                          .format(self.nproc_x))
            elif param == "nest_pes_y" and self.is_multistorm:
                self.nproc_y  = values[0:self.num_domains]
                self.nproc_y = [int(y) for y in self.nproc_y]
                log.debug("Setting nproc_y from nest_pes_y = {0}"
                          .format(self.nproc_y))

        #import pdb; pdb.set_trace()
        self._set_start_and_end_dates(start_time, end_time, run_time)
        # Now that we know num_domains, we can set self.parent_ids
        self.parent_ids = parent_ids[0:self.num_domains]

    def _set_params_from_rsl_files(self):

        # Local vars
        forecast_succeeded = False
        desired_duration_lt_actual = False

        # open RSL file
        if os.path.exists(os.path.join(self.path,'rsl.out.0000')):
            rsl_file = open(os.path.join(self.path,'rsl.out.0000'))
        elif os.path.exists(os.path.join(self.path, 'rsl.out.0000.gz')):
            rsl_file = gzip.open(os.path.join(self.path, 'rsl.out.0000.gz'))
        else:
            raise Exception("Could not find rsl.out.0000 or rsl.out.0000.gz in %s" %self.path)

        #
        # variables for tracking cummulative times
        #
        # Counter that tracks outputs that are done hourly for each domain
        # (i.e. exch_tim, sum_tim, radiation_tim, etc.)
        # - tracks how many times _each_ domain's output has been printed
        cumm_hourly_output_ctr = 0
        # variable to keep track of cummulative exchange times
        prev_cumm_exch_time = np.zeros(self.num_domains+1, dtype=np.float)
        # variable to keep track of the current domain being reported in the
        # hourly timing output, since they don't specify the domain in each line
        curr_timer_dom = 0

        if self.num_domains > 3:
            log.warn("More than 3 domains...this has not been tested extensively."
                      " Known caveats: (1) Timer information (e.g. exch_tim)"
                      " is only up to domain 3"
                     .format())
        # Begin long loop through RSL file
        for line in rsl_file.readlines():
            
            toks = line.strip().split()
            #56 I/O server 2 is ready for operations.
            if toks[0:2] == ['I/O', 'server']:
                self.io_server_selections.append(int(toks[2]))

            # 'Timing for writing' outputs
            elif toks[0:3] == ['Timing', 'for', 'Writing']:

                #
                # Determine domain number
                #
                if toks[3] == 'filter':
                    domNum = int(toks[7][:-1]) # remove trailing :
                else:
                    domNum = int(toks[6][:-1]) # remove trailing :

                #
                # process output times for different output file types
                #
                # Timing for Writing wrfout_d01_2005-08-01_12:00:00 for domain 1: 1.95625
                if toks[3].startswith("wrfout"):
                    self.hist_output_times[domNum].append(float(toks[7]))
                # Timing for Writing gsi_wrfhwrf_d01_2005-08-01_12: for domain 1: 0.88487
                elif toks[3].startswith("gsi_wrfhwrf"):
                    self.gsi_auxhist_output_times[domNum].append(float(toks[7]))
                # Timing for Writing filter output for domain 1: 1.29152 seconds
                elif toks[3] == 'filter':
                    self.filter_times[domNum].append(float(toks[8]))

            #
            # Process the timers built into HWRF. These show the _cummulative_ execution times
            # for various functions.
            # These are hardcoded to be printed every hour in solve_nmm.F
            #
            elif toks[0] == 'grid%ntsd=' and toks[2] == 'solve_tim=':
                # This is the first timer printed in HWRF 3.5, so increment
                # the domain counter (the output will be printed in a separate
                # segment for each domain
                #[ Verify for new HWRF versions]
                curr_timer_dom = curr_timer_dom + 1
            elif toks[0] == 'hifreq_tim=':
                # This is the last timer printed in HWRF 3.5, so reset the domain
                # counter.
                #[ Verify for new HWRF versions]
                #if currDom == self.num_domains: 
                # For multistorm, hifreq timers are only printed for d01,d02,d03
                if curr_timer_dom in (self.num_domains, 3):
                    curr_timer_dom = 0
            elif toks[0] == 'exch_tim=':
                # exch_tim=     0.083923 pct=  7.013%
                # This is the commulative time of all halos (i.e. HALO_NMM_F + HALO_NMM_I + ...).
                # The percentage shown is the percentage shown in the log file is that
                # of solve_nmm (i.e. solve_interface), so it excludes things like forcing, feedback,
                # allocating domains, med_after_solve_io, med_nest_move, etc. (see module_integrate)
                #if cumm_hourly_output_ctr % self.num_domains == 0:
                curr_cumm_exch_time = float(toks[1])
                dt = curr_cumm_exch_time - prev_cumm_exch_time[curr_timer_dom]
                prev_cumm_exch_time[curr_timer_dom] = curr_cumm_exch_time
                self.halo_times[curr_timer_dom].append(dt)
            ####
            # End parsing of HWRF built-in timers
            #######
            

            # Timing for main: time 2005-08-01_12:00:05 on domain   1:   30.89756
            # Note: timing for main includes all operations performed in module_integrate
            #       for a given time step, for all grids, including forcing and feedback
            #  -> So the time of d01 includes the times of all nests
            #  -> The vast majority of execution time should be included here, since
            #     only wrf_init and wrf_finalize are done outside of here
            elif toks[0:3] == ['Timing', 'for', 'main:']:
                domNum = int(toks[7][:-1])
                currDuration = wrf_timestamp_to_datetime(toks[4]) - self.start_date
                if currDuration > self._desired_duration:
                    # if we're over the desired duration passed in at 
                    # instantiation, there's no need to go further
                    desired_duration_lt_actual = True
                    break
                self.integrate_times[domNum].append(float(toks[8]))
            # newer style of main timing log:
            #  main: 2016-09-11_09:19:56 dom 3 <dP/dt>=5.022 hPa/3hr timing .03612 s
            elif toks[0] == "main:" and toks[2] == "dom" and toks[6] == "timing":
                domNum = int(toks[3])
                currDuration = wrf_timestamp_to_datetime(toks[1]) - self.start_date
                if currDuration > self._desired_duration:
                    desired_duration_lt_actual = True
                    break
                self.integrate_times[domNum].append(float(toks[7]))

            # Completo
            elif toks[3:6] == ['SUCCESS', 'COMPLETE', 'WRF']:
                #d01 2005-08-09_00:00:00 wrf: SUCCESS COMPLETE WRF
                # mark forecast as succeeded and update end_time
                forecast_succeeded = True
                newEndDate = wrf_timestamp_to_datetime(toks[1])
                if newEndDate != self.end_date:
                    log.warn("The end_date determined in the log file does not"\
                             " match the one in the namelist. The log file is"\
                             " a more reliable source, hence the end_date is %s"\
                             %(newEndDate) )
                self.end_date = newEndDate
            ######
            # end loop thru rsl.out.0000 readlines()
            ####

            
            # if nproc_x and nproc_y were not not set in namelist, get from RSL file
            if not 'nproc_x' in self.__dict__:
                if toks[0:3] == ['Ntasks', 'in', 'X']:
                    log.debug("Setting nproc_x and nproc_y based on RSL output")
                    if toks[3][-1] == ',': toks[3] = toks[3][0:-1]
                    if toks[7][-1] == ',': toks[3] = toks[3][0:-1]
                    self.nproc_x = int(toks[3])
                    self.nproc_y = int(toks[7])

        # Set sibling doms' timers (i.e. d04=d02...)
        for domNbr in range(4, self.num_domains+1):
            # TODO : Generalize for more than 2 pairs of siblings
            if domNbr in (4,6,8,10,12,14,16):
                haloTime = self.halo_times[2]
            elif domNbr in (5,7,9,11,13,15,17):
                haloTime = self.halo_times[3]
            self.halo_times[domNbr]  = haloTime

        #
        # Some debugging outputs that also serve as sanity checks
        #
        stats = lambda x : ( np.mean(x), np.min(x), np.max(x) )
        #print integrate_times[1]
        for domNum in range(1, self.num_domains+1):
            for param in ('integrate_times', 'hist_output_times',
                          'gsi_auxhist_output_times', 'filter_times',
                          'halo_times'):
                a = getattr(self, param)
                if len(a[domNum]) == 0:
                    log.debug("No entries in log file for parameter %s "%param)
                    continue
                log.debug("(mean,min,max) time for domain %i for %s = (%f,%f,%f)"
                           #%( domNum, param, stats(a[domNum])))
                           %( domNum, param,
                           stats(a[domNum])[0], stats(a[domNum])[1], stats(a[domNum])[2]))
        # now for params that do not have an element per domain
        #for param in ['halo_times']:
        #    dom = 1
        #    a = getattr(self,param)
        #    log.debug("(mean,min,max) time for %s = (%f,%f,%f)"
        #              %(param, stats(a[dom])[0], stats(a[dom])[1], stats(a[dom])[2]))

        #
        # sanity checks
        #
        # TODO gotta formulate this for multidom runs
        #print self.forecast_duration
        #print self.integrate_times[1]
        #print len(self.integrate_times[1])
        #print len(self.halo_times[1])
        #import pdb ; pdb.set_trace()
        num_timesteps = self.forecast_duration.total_seconds() / self.time_step
        if len(self.integrate_times[1]) != num_timesteps:
            raise Exception("Num integrate logs {0} != num_timesteps {1}"
                            .format(len(self.integrate_times[1]), num_timesteps))
        # for params that vary by dom number
        for domNum in range(1,self.num_domains+1):
            #import pdb ; pdb.set_trace()
            expected_hist_outs = self.forecast_duration.total_seconds() \
                                 / (self.history_interval[domNum] * 60) + 1
            # TODO - This needs to be modified to account for auxhist files
            #  with "wrfout" prefix in filename and/or with auxhistN_end 
            # params. e.g. with BS HWRF there is hourly wrfout_d01 output
            # for the first 9 hours and then 3 hourly
            if  len(self.hist_output_times[domNum]) != expected_hist_outs:
                #raise Exception("Num hist outs {0} != expected {1}"
                #                .format(len(self.hist_output_times[domNum]), expected_hist_outs))
                log.warn("Num hist outs {0} != expected {1}"
                         .format(len(self.hist_output_times[domNum]), expected_hist_outs))
            # exch_tim et al are hardcoded to be printed every hour,
            # they print at the 0th hour, but not at the final forecast hour
            # so if self._desired_duration is lower than the actual forecast duration,
            # the value will be  one more than if it's the same as the actual duration
            # TODO: Make this better using the desired_duration_lt_actual var
            expected_num_halos = self.forecast_duration.total_seconds() / 3600
            if not ( len(self.halo_times[domNum]) == expected_num_halos) \
                or (len(self.halo_times[domNum]) == expected_num_halos + 1):
                raise Exception("Domain {0}: Num halos {1} != expected {2}"
                         .format(domNum, len(self.halo_times[domNum]), 
                                 expected_num_halos))
        # Ensure forecast succeeded
        if not desired_duration_lt_actual and not forecast_succeeded \
           and FAIL_IF_UNCONFIRMED_FORECAST:
            log.error("Did not find 'SUCCESS COMPLETE WRF' string in output."\
                      " Cowardly bailing. Set 'FAIL_IF_UNCONFIRMED_FORECAST'"\
                      " to False to bypass this sanity check.")
            sys.exit(1)
       # now that we read the output, ensure the duration  passed in during
       # instantiation is not larger than the actual forecast duration
        if self._desired_duration > self.forecast_duration:
            raise Exception("Desired duration {0} exceeds actual forecast duration {1}"
                            .format(self._desired_duration, self.forecast_duration))

#
# Utility functions
#
def wrf_timestamp_to_datetime(timestamp):
    '''
    Given a timestamp in the usutal WRF format (YYYY-MM-DD_HH:MM:SS), 
    return the corresponding datetime object 
    '''
    ts = time.strptime(timestamp,"%Y-%m-%d_%H:%M:%S")
    return datetime(ts.tm_year,
                    ts.tm_mon,
                    ts.tm_mday,
                    ts.tm_hour,
                    ts.tm_min,
                    ts.tm_sec)

if __name__ == '__main__':
    log.basicConfig(level=log.DEBUG)
    wrfrun = WrfRun('./sample_wrf_dir', duration=32400)
    for domIdx in range(1, wrfrun.num_domains + 1):
        print 'Timing for integrate() for dom ', domIdx, ':', sum(wrfrun.integrate_times[domIdx])
    
    # The following should be very close to integrate_times[1]
    #print 'Total execution time:', open('./sample_wrf_dir/execution_time.txt').readline()
    
    #print wrfrun.forecast_duration
    #print wrfrun.io_server_selections

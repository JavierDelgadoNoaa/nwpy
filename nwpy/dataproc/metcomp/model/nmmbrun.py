#!/usr/bin/env python
'''
This module defines the NmmbRun class. Which encapsulates useful information
about an NMM-b run.
'''
import sys
import os
import re
import logging as log
import numpy as np
import gzip
from datetime import datetime, timedelta
#from . import ModelRun
from modelrun import ModelRun



class NmmbRun(ModelRun):
    
    def __init__(self, path, duration=432000, logfile_name='nmmb_fcst.log'):
        '''
        Instantiate an NmmbRun. Does the standard instantiation of the 
        parent class as well as instantiation of parameters specific to NMM-B.
        See docstring for parent (ModelRun) for important information.
        '''
        super(NmmbRun, self).__init__(path, duration)
        self._set_params_from_configure_file(domain=1)
        # Set variables specific to NMMB runs
        self.logfile_name = logfile_name
        # communications performed during physics
        self.physics_comms = {}
        # communications performed during dynamics 
        self.dynamics_comms = {}        
        
        for domNum in range(1, self.num_domains+1):
            self.integrate_times[domNum] = []
            self.halo_times[domNum] = []
            self.physics_comms[domNum] = []
            self.dynamics_comms[domNum] = []
            
        self._set_params_from_output()
            
    def _set_params_from_configure_file(self, domain):
        '''
        Parses the configura file for the given domain and
        sets member variables accordingly.
        If domain == 1, call this method recursively according to 
        the number of domains specified in configure_file_01
        '''
        curr_cfg = 'configure_file_' + str(domain).zfill(2)
        with open(os.path.join(self.path,curr_cfg)) as cfg_file:
            pv_list = {}
            for line in cfg_file.readlines():
                if line.startswith("#"): continue
                toks = line.split(":")
                if len(toks) < 2: continue
                param = toks[0].strip()
                #values = [ v.strip() for v in toks[1].split(",")]
                values = toks[1]#.strip("\n")
                # strip trailing comments
                commentIdx = values.find('#')
                #if commentIdx > -1 : print values ; import pdb ; pdb.set_trace()
                values = values[0:commentIdx]
                #if line.startswith("start"): import pdb ; pdb.set_trace()
                pv_list[param] = values
                #print values
        if int(domain) == 1 and not pv_list.has_key('num_domains_total'):
            raise Exception("num_domains_total not found in configure_file_01")
        self.num_domains = int(pv_list['num_domains_total'])
            
        if domain == 1:
           # For d01, instantiate variables that have an element for each domain
           # and call this routine for the other domains
            self.nproc_x = {}
            self.nproc_y = {}
            self.gridsize_we = {}
            self.gridsize_ns = {}
            for d in range(2, self.num_domains):
                self._set_params_from_config_file(d)
        
        # Now set the values
        to_bool = lambda x : str(x).lower() in ('true', '1', 1)
        start_params = {} # dict to termporarily store all the start_date params
        for (param,value) in pv_list.iteritems():
            if param == 'inpes':
                #self.nproc_x[domain] = int(value)
                self.nproc_x = int(value)
            elif param == 'jnpes':
                log.warn("TODO : Make nproc_x and y a dict mapping domain to num_compute_tasks for it")
                #self.nproc_y[domain] = int(value)
                self.nproc_y = int(value)
            elif param == 'dt_int':
                self._dt_sec = int(value)
            elif param == 'dt_num': 
                self._dt_frac_num = int(value)
            elif param == 'dt_den':
                self._dt_frac_den = int(value)
            elif param == 'im':
                self.gridsize_we[domain] = int(value)
            elif param == 'jm':
                self.gridsize_ns[domain] = int(value)
            elif param == 'quilting':
                self.quilting = to_bool(value)
            elif param == 'write_groups':
                self.write_groups = value
            elif param == 'write_tasks_per_group':
                self.write_tasks_per_group = value
            elif param.startswith('start_'):
                print param, value
                start_params[param] = int(value)
            elif param == 'tstart':
                self.start_fhr = value
            elif param == 'nhours_fcst':
                # the actual forecast duration is set as an attribute in parent
                _forecast_duration = float(value) * 3600.
            elif param == 'minutes_history':
                self.history_interval = float(value) * 60.
            elif param == 'nrads':
                self.nrads = int(value) #  #dynamics timesteps between calls to swrad
            elif param == 'nradl':
                self.nradl = int(value) #  #dynamics timesteps between calls to lwrad
            elif param == 'nphs':
                self.nphs = int(value) #   ... calls to landsurface and turbulence
            elif param == 'nprecip':
                self.nprecip = int(value) # ... calls to convection and microphysics
            elif param == 'my_domain_moves':
                self.is_static = to_bool(value)
            elif param == 'generation':
                self.generation = int(value)
            elif param == 'nest_mode':
                self.nest_mode = value

        # set the start_date
        self.start_date = datetime(start_params['start_year'],
                                   start_params['start_month'],
                                   start_params['start_day'],
                                   start_params['start_hour'],
                                   start_params['start_minute'],
                                   start_params['start_second'])

        self.end_date = self.start_date + timedelta(seconds=_forecast_duration)


    def _set_params_from_output(self):
        '''
        Sets the following attributes from the output file:
            - integration_times - list holding integration times for each time step
            - (communication_times) - if communication times are output, this list
                      will be created and it will hold communication times for
                      each time step.
        '''
        if self.logfile_name.endswith(".gz"):
           log_file = gzip.open(os.path.join(self.path,self.logfile_name), 'r')
        else:
            log_file = open(os.path.join(self.path, self.logfile_name), 'r')
        # Read entire file now
        try:
            log_file_lines = log_file.readlines()
        except:
            log.error("Unable to read log file")
            sys.exit(2)
        log_file.close()
        
        #
        # variables for tracking cummulative times
        #
        # Counter that tracks outputs that are done hourly for each domain
        # (i.e. exch_tim, sum_tim, radiation_tim, etc.)
        # - tracks how many times _each_ domain's output has been printed
        cumm_hourly_output_ctr = 0
        # variables to keep track of cummulative exchange times
        prev_cumm_phys_exch_times  = np.zeros(self.num_domains+1, dtype=np.float)
        prev_cumm_dyn_exch_times = np.zeros(self.num_domains+1, dtype=np.float)
        prev_cumm_exch_time = np.zeros(self.num_domains+1, dtype=np.float)
        # variable to keep track of the current domain - for outputs that do
        # not specify the domain in the output itself
        iCurrDom = 0

        for line in log_file_lines:
            # BEGIN: Iterate through lines of log file
            # We usually use -prepend-rank (or Intel equivalent), so check for 
            # cases where it was not used and append dummy element there
            m = re.search('^\[[0-9]*\]', line)
            if m is None:
                line = '[dummy] %s' %line
            toks = line.strip().split()
            if len(toks) == 1: continue
            #print toks
            #print line
            #import pdb ; pdb.set_trace()
            if toks[1:3] == ['Finished', 'Timestep']:
                #[0]  Finished Timestep 11327 for domain   1 ending at
                #    9.440 hours: elapsed integration time   0.05037
                if float(toks[9]) * 3600  > self._desired_duration:
                    # if we're over the desired duration passed in at 
                    # instantiation, there's no need to go further
                    break
                self.integrate_times[int(toks[6])].append(float(toks[14]))
            elif toks[1:4] == ['Clocktimes', 'for', 'domain']:
                #[1]  Clocktimes for domain #01
                iCurrDom = int(toks[4][1:])
            elif toks[1] == 'exch_dyn=':
                # [1]    exch_dyn=               377.46     pct=   11.26
                curr_cumm_t = float(toks[2])
                dt = curr_cumm_t - prev_cumm_dyn_exch_times[iCurrDom]
                prev_cumm_dyn_exch_times[iCurrDom] = curr_cumm_t 
                self.dynamics_comms[iCurrDom].append(dt)
            elif toks[1] == 'exch_phy=':
                #[1]    exch_phy=               84.798     pct=    2.53
                curr_cumm_t = float(toks[2])
                dt = curr_cumm_t - prev_cumm_dyn_exch_times[iCurrDom]
                prev_cumm_phys_exch_times[iCurrDom] = curr_cumm_t 
                self.physics_comms[iCurrDom].append(dt)
            elif toks[1] == 'exch_tim=':
                #[1]    exch_tim=               462.25     pct=   13.79
                curr_cumm_t = float(toks[2])
                dt = curr_cumm_t - prev_cumm_exch_time[iCurrDom]
                prev_cumm_exch_time[iCurrDom] = curr_cumm_t
                self.halo_times[iCurrDom].append(dt)
            elif toks[3:6] == ['Write', 'time', 'is']:
                # for history:  Write Time is    3.19939638790601       
                #  at Fcst           12 :           0 :
                #  -> only gets output if 'print_output' or 'print_all' is true 
                #     in configure_file_NN
                # TODO : See how this works for multiple domains
                self.history_output_times.append(float(toks[6]))
        # now that we read the output, ensure the duration  passed in during
        # instantiation is not larger than the actual forecast duration
        if self._desired_duration > self.forecast_duration:
            raise Exception("Desired duration %i exceeds actual forecast duration %i"
                            %(self._desired_duration, self.forecast_duration))
        
if __name__ == '__main__':
    log.basicConfig(level=log.DEBUG)
    try:
        err_run = NmmbRun('./sample_nmmb_dir', duration=1e6)
    except:
        print ';) Successfully caught error'
    nmmb_run = NmmbRun('./sample_nmmb_dir', duration=345600)
    for domIdx in range(1, nmmb_run.num_domains + 1):
        print 'Timing for integrate() for dom ', domIdx, ':', sum(nmmb_run.integrate_times[domIdx])
        print 'Timing for halos (exch_tim) for dom ', domIdx, ':', sum(nmmb_run.halo_times[domIdx])
        print nmmb_run.get_x_tile_size(domIdx)
        print nmmb_run.y_tile_size[domIdx]
    

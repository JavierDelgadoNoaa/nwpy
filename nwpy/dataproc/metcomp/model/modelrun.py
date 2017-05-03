'''
This module contains the definition of the ModelRun, which provides methods
to extract various information from model run directories. The ModelRun
is a an abstract class containing common/required functions. 
'''

import sys
import os
import logging as log
import numpy as np
import math

class ModelRun(object):
    def __init__(self, path, duration=432000):
        '''
        Instantiate an ModelRun. This will set the parameters that are
        required for all derived types, which are the following:
          -path - path containing model output data
          -duration - how much of the forecast data to account for
          -hist_output_times - time for writing the main "history" files
          -halo_times - time spent exchanging data for halo points. 
          -integrate_times - time spent integrating. 
        These last 3 fields are dictionaries wherein the keys correspond to 
        grids/nests and each key has a corresponding array of times.
        All times are in seconds. If the forecast duration is greater than
        the passed in `duration' (default 432k seconds = 120 hours), only 
        up the passed in `duration' will be read. If the passed in `duration'
        exceeds that actual forecast duration, an Exception will be thrown
        after reading the log file.
        '''
        self.path = path
        if not os.path.exists(path):
            raise Exception("Invalid Path given: %s" %path)
        self._desired_duration = duration

        # Instantiate variables that must be present
        # for all derived types
        
        # write time for main history files 
        self.hist_output_times = {}
        # halo times
        self.halo_times = {}
        # 'timing for main' in WRF
        self.integrate_times = {}

    def get_x_tile_size(self, dom_num):
        ''' Max number of "tiles"/"subgrids" along the west-east axis 
            for domain number `dom_num' per compute task '''
        dom_num = int(dom_num)
        return int(math.ceil(
            float(self.gridsize_we[dom_num]) /  self.nproc_x[dom_num])) 
    
    def get_y_tile_size(self, dom_num):
        ''' Max number of "tiles"/"subgrids" per compute task along the
            north-south axis for domain number `dom_num' '''
        dom_num = int(dom_num)
        return int(math.ceil(
            float(self.gridsize_ns[dom_num]) /  self.nproc_y[dom_num])) 

    @property
    def y_tile_size(self):
        ''' Max number of "tiles"/"subgrids" per compute task along the
            north-south axis for domain number `dom_num' '''
        d = {}
        a = [int(math.ceil(
            float(self.gridsize_ns[dom_num]) /  self.nproc_y[dom_num])) 
            for dom_num in range(1,self.num_domains + 1)] 
        for k,v in enumerate(a):
           d[k+1] = v
        return d
    @property
    def forecast_duration(self):
        ''' Forecast duration in seconds '''
        timeDelta = self.end_date - self.start_date
        return timeDelta.days * 3600 * 24 + timeDelta.seconds

    @property
    def time_step(self):
        ''' Integration time step in seconds '''
        return float(self._dt_sec) + float(self._dt_frac_num) / float(self._dt_frac_den)

    @property
    def num_compute_tasks(self):
        ''' Number of processes used for computation - nproc_x * nproc_y '''
        return self.nproc_x * self.nproc_y

    def add_scalability_datapoint_to_axes(self, ax, dom=None, forecast_duration=None, 
                                           x_axis='total_procs', **kwargs):
        '''
        Plot a point in the given axes corresponding to 
        the number of processors used and the execution time.
        The following optional parameters can be used:
            -dom - The domain to get times for. None means get sum of times for
                    all domains.
            -forecast_duration (default: self.forecast_duration): For what forecast 
             duration? (Note this will be an approximation based on stdout)
            -x_axis (default: total_procs; options: nproc_x, nproc_y) - What
             parameter to use for the x axis
            -kwargs - arguments to MPL's ax.plot()
        '''
        if forecast_duration != None:
            pass
            #y = calculate based on start time in header and time step counter
        else:
            y = self.execution_time
        
        if x_axis == 'total_procs':
            x = self.nproc_x[dom] * self.nproc_y[dom]
        elif x_axis == 'nproc_x':
            x = self.nproc_x[dom]
        elif x_axis == 'nproc_y':
            x  = self.nproc_y[dom]
        else:
            raise Exception("Invalid value passed for x_axis")
            
        ax.plot(x,y, **kwargs)
                 

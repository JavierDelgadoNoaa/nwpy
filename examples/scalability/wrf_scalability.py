#!/usr/bin/env python
"""
Plot one or more metrics from a set of WRF run paths.
Most of the magic is done using the  nwpy.viz.performance.model_plot.ModelPerformancePlotter.
This script takes a pattern of paths and plots data from each path matching the pattern
"""

import os
import sys
from glob import glob
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import plot, gca, figure, show, savefig, legend, tight_layout, xlim
import logging as log
import numpy as np
from datetime import datetime, timedelta

from nwpy.dataproc.metcomp.model.wrfrun import WrfRun
from nwpy.viz.performance.model_plot import ModelPerformancePlotter

if __name__ ==  '__main__':
    
    ##
    # SETTINGS
    ##
    
    #log.basicConfig(level=log.DEBUG)
    log.basicConfig(level=log.INFO)
        
    domain_of_interest = 1
    #fcst_duration = 120 * 3600
    fcst_duration = timedelta(hours=60)
    
    plotter = ModelPerformancePlotter(x_axis_metric = "num_compute_tasks", 
                                      x_axis_label = "#cores",
                                      y_axis_metrics = ["halo_times", "integrate_times"], 
                                      y_axis_labels = ["halo times", "integration times"],
                                      # Reduction functions to apply the y_axis metric(s).
                                      # e.g. if you want the total amount of time spent
                                      # integrating, use np.sum along with integrate_times
                                      y_axis_reducers = [np.sum, np.sum],
                                      # Which index of the metric to use. Generally, this means
                                      # "which domain to plot the metric's value for"
                                      y_axis_keys = [1,1])
    '''
    plotter.set_plot_options_for_groups(
        { 'wrf-noNest' : 
            { 'marker' : 'x' , 'color' : 'blue' }  
        , 'wrf-nest' :
            { 'marker' : '+' , 'color' : 'red' }
        })
    '''
    #path_pattern = "/home/Javier.Delgado/scratch_lfs2/pytmp_2016-rev3_history-jza/hwrf-basinscale_multistorms-rev3_2016_history-jza/2016091018/00L/run_sjet/cfg_two"
    path_pattern = "/home/Javier.Delgado/scratch_lfs2/pytmp_2016-rev3_history-jza/hwrf-basinscale_multistorms-rev3_2016_history-jza/2016091018/00L/run_sjet/cfg*"
  

    ##
    # LOGIC
    ##
    cfg = 3 # 1-2 do not have rsl files
    for wrfPath in glob(path_pattern):
        #import pdb ;pdb.set_trace()
        wrfRun = WrfRun(wrfPath, duration=fcst_duration)
        #groupName = wrfRun.write_tasks_per_group # not implemented - group paths in glob according to the number of i/o tasks/group
        #plotter.add_datapoint(wrfRun, groupName)
        label = "cfg" + str(ctr) # just one dataset
        ctr += 1
        plotter.add_datapoint(wrfRun, label)
    plotter.create_line_plot()
    legend(loc='best')
    #show()
    savefig("wrf_basin_scal.png")

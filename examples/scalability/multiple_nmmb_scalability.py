#!/usr/bin/env python

import os
import sys
from glob import glob
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import plot, gca, figure, show, savefig, legend, tight_layout, xlim
import logging as log
import numpy as np
from wrfrun import WrfRun
from nmmbrun import NmmbRun
from model_plot import ModelPerformancePlotter

if __name__ ==  '__main__':
    
    log.basicConfig(level=log.DEBUG)
        
    domain_of_interest = 1
    #fcst_duration = 120 * 3600
    fcst_duration = 345600
    
    plotter = ModelPerformancePlotter(x_axis_metric = "num_compute_tasks", 
                                      x_axis_label = "#cores",
                                      y_axis_metrics = ["halo_times"], 
                                      y_axis_labels = ["halo times"],
                                      y_axis_reducers = [np.sum],
                                      y_axis_keys = [1])
    '''
    plotter.set_plot_options_for_groups(
        { 'wrf-noNest' : 
            { 'marker' : 'x' , 'color' : 'blue' }  
        , 'wrf-nest' :
            { 'marker' : '+' , 'color' : 'red' }
        })
    '''
    for nmmbPath in glob("/home/Javier.Delgado/nmmb_3km_scalability/output*"):
        nmmbRun = NmmbRun(nmmbPath, duration=fcst_duration, logfile_name='nmmb_fcst.log.gz')
        groupName = nmmbRun.write_tasks_per_group
        plotter.add_datapoint(nmmbRun, groupName)
    plotter.create_line_plot()
    legend(loc='best')
    #show()
    savefig("nmmb_3km_scal.png")

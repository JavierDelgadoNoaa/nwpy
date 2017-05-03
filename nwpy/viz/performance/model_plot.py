#!/usr/bin/env python

import os
import sys
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import plot, gca, figure, show, savefig, legend, tight_layout, xlim, bar
import numpy as np
import logging as log 
from pycane.met_components.modelrun import ModelRun

'''
This module contains utilities for creating plots based on ModelRun data. 
'''

class ModelPerformancePlotter(object):
    '''
    Object to facilitate the generation of plots that show performance metrics
    of a ModelRun. 
    Basic Usage:
     1. Instantiate a ModelPerformancePlotter, defining the x-axis metric and 
        y-axis metric(s) that you would like to plot. The chosen values will
        then be queried from the ModelRun objects passed to the add_datapoint
        method.
     2. Add `Datapoint's to it using the add_datapoint object. Distinguish
        different datasets using the 'label' attribute to this method. Each
        label will be assigned a different color/marker. If multiple y_axis_metric
        values are used, each will get a different color/marker.
     3. Create the plot using method corresponding to the plot you would like
        (or create your own)
     ** NOTE: The X and Y arrays containing the data to be plotted are only 
        generated once - on the first call to a plotting method. So adding 
        data points after calling a plotting method will have no effect on 
        future plots ***
    '''
    
    def __init__(self, x_axis_metric, x_axis_label, y_axis_metrics, y_axis_labels,
                 y_axis_reducers, y_axis_keys):
        '''
        PARAMS
         * x_axis_metric - Metric to use for the x-axis. This should correspond
                           to a field of the ModelRun objects being used.
                           Should be a scalar. 
                           (e.g. 'num_compute_tasks', 'y_tile_size')
         * x_axis_label - Label to give the x-axis (e.g. 'tile size', '#cores')
         * y_axis_metrics - List of metrics to be used for the y-axis. They should
                            correspond to member vars of the ModelRun objects being 
                            used. They can be lists or scalars. In the case of the 
                            former, the corresponding y_axis_reducer will be used
                            to reduce the list to a scalar for plots that expect
                            scalars. 
                            (e.g. [integration_times], [halo_times, compute_times] )
        * y_axis_labels - List of labels corresponding to the `y_axis_metrics'
        * y_axis_reducers - List of reduction methods. Each element corresponds to the
                            method that should be applied to the same index in `y_axis_metrics'. This is needed because a lot of the data
                            is stored in the ModelRun object as a list. This can be any
                            method that works on a list or None if no reduction is
                            necessary.
                            (e.g. sum, np.mean, None)
        * y_axis_keys - Many of the fields in the ModelRun are dictionaries that map a
                       domain/nest to its data. If this is the case, pass the key 
                       in this parameter. In other cases (i.e. when the field is a 
                       scalar or list with all the data, use None for the 
                       corresponding index.
        '''
        assert len(y_axis_labels) == len(y_axis_metrics) == len(y_axis_reducers) \
                                  == len(y_axis_keys)
        self.x_axis_metric = x_axis_metric
        self.x_axis_label = x_axis_label
        self.y_axis_metrics = y_axis_metrics
        self.y_axis_labels = y_axis_labels
        self.y_axis_reducers = y_axis_reducers
        self.y_axis_keys = y_axis_keys
        self.datapoints = []
        self._xvals = None
        self._yvals = None
        # These will be used for groups/labels/metrics that do not have them 
        # specified vi group_plot_settings or metric_plot_settings
        self._hatch_list = [ '/' , '+' , ',' , '-' , '\\' , 'x' , 'o' , 'O' , 
                             '.' , '*' ]
        self._color_list = [ 'red', 'green', 'blue', 'orange', 'cyan', 
                             'purple', "#006600", "#9999ff", "#ff0099",
                             "666666", 'cccc99' ] 
        self._marker_list = [ 'o', 'v', '^', '<', '>', '8', 's', 'p', '*',
                              'x', 'D', '+', 'H' ]
        self._color_generator = self._init_color_iterator()
        self._hatch_generator = self._init_hatch_generator()
        self._marker_generator = self._init_marker_generator()
        self.group_plot_settings = {}
        self.metric_plot_settings = {}                      
                              

    def set_plot_options_for_groups(self, settings):
        '''
        Sets plotting attributes for different groups/labels being plotted.
        The argument `settings' should be a dictionary that maps the 
        group name/label to a dictionary containing plotting settings. The 
        following keys are supported in this inner dictionary:
            color - marker color
            marker - type of marker
            markersize - size of the marker
            hatch - string with desired hatch pattern
            linecolor - color for line
        See the Matplotlib documentation for plot() and boxplot() for acceptable
        values for these.
        '''
        self.group_plot_settings = settings
            
    def set_plot_options_for_metrics(self, settings):
        '''
        Sets plotting attributes for the different 'y_metrics'.
        The argument `settings' should be a dictionary that maps the 
        metric to a dictionary containing plotting settings. The
        following keys are supported in this inner dictionary:
            color - marker color 
            marker - the type of marker
            markersize - size of the marker 
            hatch - string with the desired hatch pattern (e.g. for boxplots)
            linecolor - color for the line 
        See the Matplotlib documentation for plot() and boxplot() for acceptable
        values for these.
        '''         
        self.metric_plot_settings = settings
        
    def add_datapoint(self, modelrun, label):
        '''
        Adds a datapoint based on the passed in `modelrun' to the current Plotter.
        The passed in `label' will be applied to the Datapoint
        '''
        missing_attr_msg = lambda path, attr : \
            "The given ModelRun object at path %s does not contain " \
            "a member variable called %s - skipping this datapoint" \
            %(path, attr)
        y = {}
        try:
            x = getattr(modelrun, self.x_axis_metric)
        except:
            log.warn(missing_attr_msg(modelrun.path, self.x_axis_metric) )
            return
        for i, yMetric in enumerate(self.y_axis_metrics):
            try:
                #import pdb ; pdb.set_trace()
                y[yMetric] = getattr(modelrun, yMetric)
                log.debug("Adding y-metric '%s' for data point label '%s'. Value = '%s'"
                          %(yMetric, label, y[yMetric]))
                if self.y_axis_keys[i] is not None:
                    try:
                        key = self.y_axis_keys[i]
                        y[yMetric] = y[yMetric][key]
                    except:
                        log.error("Key %s does not exist for y-axis-metric %s"
                                  %(key, yMetric) ) 
            except:
                log.warn(missing_attr_msg(modelrun.path, self.x_axis_metric))
        # Note that this is only done if all values are present. Otherwise
        # we gotta deal with the possibility of NaNs
        self.datapoints.append(Datapoint(x, y, label))
        
    def create_line_plot(self):
        '''
        Create a scatter plot of the self.x_values and self.y_values for
        each metric in self.y_values
        '''
        unique_labels = set([dp.label for dp in self.datapoints])
        for label in unique_labels:
            for i,yMetric in enumerate(self.y_axis_metrics):
                #import pdb ; pdb.set_trace()
                #print self.y_vals[label]
                # self.y_vals[label][yMetric] is an array of arrays (whose
# length is the number of datapoints with the same label
                #yVals = self.y_axis_reducers[i](self.y_vals[label][yMetric])
                yVals = [ self.y_axis_reducers[i](v) for v in self.y_vals[label][yMetric]]
                #print 'x, y = ', self.x_vals[label], yVals
                # if plotting a single group, use the given yMetric color, 
                # otherwise use the group color
                if len(unique_labels) == 1:
                    plotColor = self.get_metric_color(yMetric)
                else:
                    plotColor = self.get_group_color(label)
                # if plotting a single metric, use the group marker, otherwise
                # use the metric's marker
                if len(self.y_axis_metrics) == 1:
                    plotMarker = self.get_group_marker(label)
                else:
                    plotMarker = self.get_metric_marker(yMetric)
                
                plot(self.x_vals[label], yVals, 
                     color=plotColor,
                     marker=plotMarker,
                     linestyle='',
                     label=label+" "+yMetric) 
        tight_layout(pad=1)
        xlim(min([min(vals) for k,vals in self.x_vals.iteritems()]) - 5,
             max([max(vals) for k,vals in self.x_vals.iteritems()]) + 5)

    def create_stacked_bar_plot(self, first_alpha=1.0, alpha_reduction_step=0.25):
        '''
        Creates a "stacked" bar chart, wherein each of the y_axis_metrics 
        are stacked on top of each other. Supports plotting of arbitrary 
        number of y-axis metrics.
        TODO : group experiments with the same label and x values and 
         distinguish them using the 'yerr' parameter to bar()
         This is done since bar plots are not good for showing too many
         different data sets
        The following decoration rules apply:
         - The hatch type used always corresponds to the y_axis_metric
         - The fill color used always corresponds to the group
          -> Alpha level will be modified to help destinguish different y 
             metrics within the group. 
        OPTIONAL ARGS
         - first_alpha - Alpha level for the first y-metric. Subsequent y metrics'
           alpha values will be reduced by `alpha_reduction_step'
         - alpha_reduction_step - fraction to reduce alpha level by for each
           subsequent y metric plotted.
        '''        
        unique_labels = set([dp.label for dp in self.datapoints])
        # yVals will be an array containing the y values for each y_axis_metric
        yVals = []
        for label in unique_labels:
            alpha = first_alpha
            for i,yMetric in enumerate(self.y_axis_metrics):
                yVals.append([self.y_axis_reducers[i](v) for v in self.y_vals[label][yMetric]])
                # Determine what yMetric was under this one
                if i == 0: 
                    bottom = None
                else: 
                    bottom = yVals[i-1]
                # Determine colors
                hatchType = self.get_metric_hatch(yMetric)
                fillColor = self.get_group_color(label)
                # Set the width of the bar according to x range
                # TODO : make this a function of the number of labels rather
                # than hardcoding to 20
                minX = min([min(vals) for k,vals in self.x_vals.iteritems()])
                maxX = max([max(vals) for k,vals in self.x_vals.iteritems()])
                barWidth = (maxX - minX) / 20
                
                legendLabel = yMetric + ", " + label 
                
                # progressively decrease the alpha level
                alpha = alpha - (alpha * alpha_reduction_step)
                
                bar(self.x_vals[label], yVals[i],
                    bottom = bottom,
                    width = barWidth,
                    color=fillColor,
                    edgecolor='black',
                    hatch=hatchType,
                    alpha=alpha,
                    label=legendLabel)
        #ax = gca()
        #ax.set_xticks(set([arr for arr in self._xvals]))
    
    @property
    def x_vals(self):
        ''' Dictionary mapping each label to a list of x-axis values '''
        if self._xvals is None:
            self._xvals = {}
            unique_labels = set([dp.label for dp in self.datapoints])
            for label in unique_labels:
                self._xvals[label] = [dp.x for dp in self.datapoints 
                                      if dp.label == label]
        return self._xvals
        
    @property
    def y_vals(self):
        ''' 
        Dictionary mapping each label to a dictionary that maps each 
        y_axis_metric to a list. 
        The list contains the y values corresponding to the x values.
        e.g. { mvapich : { commTimes : [...] , compTimes : [...] }
                 intel : { commTimes : [...] , compTimes : [...] } 
             }
        ''' 
        if self._yvals is None:
            self._yvals = {}
            #import pdb ; pdb.set_trace()
            unique_labels = set([dp.label for dp in self.datapoints])
            for label in unique_labels:
                self._yvals[label] = {}
                for metric in self.y_axis_metrics:
                    self._yvals[label][metric] = [dp.y[metric] for dp in self.datapoints
                                                  if dp.label == label]
                    # this should be a list of lists if the y_metric is a list, with
                    # the length of the outter list being the number of data points.
        return self._yvals

    
    def get_group_color(self, group):
        #import pdb ; pdb.set_trace()
        if not self.group_plot_settings.has_key(group):
            self.group_plot_settings[group] = {}
        return self._get_key_for_dict('color',self.group_plot_settings[group],
                                      self._color_generator)
    def get_metric_color(self, metric):
        if not self.metric_plot_settings.has_key(metric):
            self.metric_plot_settings[metric] = {}
        return self._get_key_for_dict('color',self.metric_plot_settings[metric],
                                      self._color_generator)
    def get_group_marker(self, group):
        if not self.group_plot_settings.has_key(group):
            self.group_plot_settings[group] = {}
        return self._get_key_for_dict('marker', self.group_plot_settings[group],
                                      self._marker_generator)
    def get_metric_marker(self, metric):
        if not self.metric_plot_settings.has_key(metric):
            self.metric_plot_settings[metric] = {}
        return self._get_key_for_dict('marker', self.metric_plot_settings[metric],
                                      self._marker_generator)
    def get_metric_hatch(self, metric):
        if not self.metric_plot_settings.has_key(metric):
            self.metric_plot_settings[metric] = {}
        return self._get_key_for_dict('hatch', self.metric_plot_settings[metric],
                                      self._hatch_generator)
                                      
    def _get_key_for_dict(self, key, el_dict, el_generator):
        '''
        Returns the value for the given `key' in dictionary `el_dict'. If
        the key does not exist, it is added using the value given by
        `el_generator''s next() method
        '''
        if not el_dict.has_key(key):
            el_dict[key] = el_generator.next() 
        return el_dict[key]
    
    def _init_color_iterator(self):
        '''
        Generator that returns a color from self._color_list sequentially
        and cyclically
        '''
        for i in range(len(self._color_list)):
            if i > len(self._color_list):
                i = 0
            yield self._color_list[i]
            
    def _init_hatch_generator(self):
        '''
        Generator that returns a color from self._hatch_list sequentially
        and cyclically
        '''
        for i in range(len(self._hatch_list)):
            if i > len(self._hatch_list):
                i = 0
            yield self._hatch_list[i]

    def _init_marker_generator(self):
        '''
        Generator that returns a color from self._marker_list sequentially
        and cyclically
        '''
        for i in range(len(self._marker_list)):
            if i > len(self._marker_list):
                i = 0
            yield self._marker_list[i]

class Datapoint(object):
    '''
    Encapsulates attributes of a datapoint to be plotted
    '''
    def __init__(self, x, y, label):
        '''
        PARAMS
         - x - The x-axis value for this datapoint
         - y - dictionary mapping y-axis labels to their values
         - label - The label for the dataset corresponding to this point
        '''
        self.x = x
        self.y = y
        self.label = label 
        
if __name__ ==  '__main__':
    from wrfrun import WrfRun
    from nmmbrun import NmmbRun
    
    log.basicConfig(level=log.DEBUG)
    
    domain_of_interest = 1
    fcst_duration = 9 * 3600
    plotter = ModelPerformancePlotter(x_axis_metric = "num_compute_tasks", 
                                      x_axis_label = "#cores",
                                      y_axis_metrics = ["halo_times", "integrate_times"], 
                                      y_axis_labels = ["halo times", 'integration times'],
                                      y_axis_reducers = [np.sum, np.sum],
                                      y_axis_keys = [1,1])
    plotter.set_plot_options_for_groups(
        { 'wrf-noNest' : 
            { 'marker' : 'x' , 'color' : 'blue' }  
        , 'wrf-nest' :
            { 'marker' : '+' , 'color' : 'red' }
        })
    wrfrun = WrfRun('./sample_1dom_hwrf_dir', duration=fcst_duration)
    plotter.add_datapoint(wrfrun, 'wrf-noNest')
    wrfrun = WrfRun('./sample_2dom_hwrf_dir', duration=fcst_duration)
    plotter.add_datapoint(wrfrun, 'wrf-nest')
    #plotter.add_datapoint(wrfrun, 'wrf-nest') # -> breaks the bar plot
    #plotter.create_line_plot()
    plotter.create_stacked_bar_plot()
    legend(loc='best')
    #show()
    savefig("scal.png")

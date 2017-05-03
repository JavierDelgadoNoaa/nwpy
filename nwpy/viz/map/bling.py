#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
'''
This module contains functions for making maps look nice by 
selecting reasonable default values and accepting optional
values.
'''

DEFAULT_WE_GRIDLINE_FREQ = 10
DEFAULT_NS_GRIDLINE_FREQ = 10

def decorate_map(m, **kwargs):
    '''
    Decorate a Basemap according to passed in arguments.
    
    If a keyword has a value of None, the decoration will not be applied. 
    If a keyword has a value that is not None, it will be used for the deocration
    If a keyword does not have a value, the default value defined in the module will be used.

    INPUT 
     m - The Basemap object to be decorated

    The following keyword arguments are supported:
     - ns_gridline_freq - Interval (in degrees) at which to show north->south
                          gridlines
     - we_gridline_freq - Interval (in degrees) at which to show west->east
                          gridlines
     - draw_coastlines - Draw lines separating water from land? (Default: True)
     - draw_countries - Show the countries? (Default: True)
     - water_color - Colorize bodies of water? Default: Just leave it
                     as background color (white)
     - lake_color - Color to use for lakes within continents. Default: Leave as
                    background color.
     - continent_color - Color to use for continents. Default: white
     - country_outlines_color = gray
     - country_outlines_width
     *NOTE: For b/w figure, only set draw_coastlines to True

    '''
    
    key_val = lambda d,k,default: d[k] if k in d.keys() else default
    # little func to set label color, since textcolor arg is new:
    # https://github.com/matplotlib/basemap/issues/145
    def setlabelcolor(x, color):
        for m in x:
            for t in x[m][1]:
                t.set_color(color)
    # aliases
    if kwargs.has_key("ocean_color"):
        kwargs["water_color"] = kwargs["ocean_color"]
    if kwargs.has_key("continents_fill_color"):
        kwargs["continents_color"] = kwargs["continents_fill_color"]
    if kwargs.has_key("latitude_line_freq"):
        kwargs["ns_gridline_freq"] = int(kwargs["latitude_line_freq"])
    if kwargs.has_key("longitude_line_freq"):
        kwargs["we_gridline_freq"] = int(kwargs["longitude_line_freq"])
    ## Defaults
    # set lat/lon lines and labels to gray if not set
    for arg in ("latitude_label_color","longitude_label_color", 
                "latitude_line_color", "longitude_line_color"):
        if not arg in kwargs:
            kwargs[arg] = "gray"

    # process args
    if kwargs.has_key('draw_coastlines') and not kwargs['draw_coastlines']:
        pass
    else:
        m.drawcoastlines()
    try:
        zorder = int(kwargs["zorder"])
    except KeyError:
        zorder = 0
    if kwargs.has_key('water_color'):
        m.drawmapboundary(fill_color=kwargs['water_color'])
    if kwargs.has_key('lake_color') or kwargs.has_key('continent_color'):
        if not kwargs.has_key('continent_color'):
            m.fillcontinents(lake_color=kwargs['lake_color'], color='white',
                             zorder=zorder)
        elif not kwargs.has_key('lake_color'):
            m.fillcontinents(color=kwargs['continent_color'], zorder=zorder)
        else:
            m.fillcontinents(color=kwargs['continent_color'],
                             lake_color=kwargs['lake_color'], zorder=zorder)

    if kwargs.has_key('ns_gridline_freq'):
        ns_gridline_freq = kwargs['ns_gridline_freq']
    else:
        ns_gridline_freq = DEFAULT_NS_GRIDLINE_FREQ
    lat_label_mask = key_val(kwargs, "latitude_label_mask", "1 0 0 0").split()
    lat_label_mask = [int(x) for x in lat_label_mask]
    assert len(lat_label_mask) == 4
    if ns_gridline_freq is not None:
        lowerLat = round(m.llcrnrlat)
        upperLat = round(m.urcrnrlat)
        foo = m.drawparallels(np.arange(lowerLat,upperLat, ns_gridline_freq),
                   labels=lat_label_mask, color=kwargs["latitude_line_color"])
        setlabelcolor(foo, kwargs["latitude_label_color"]) 
    if kwargs.has_key('we_gridline_freq'):
        we_gridline_freq = kwargs['we_gridline_freq']
    else:
        we_gridline_freq = DEFAULT_WE_GRIDLINE_FREQ
    lon_label_mask = key_val(kwargs, "longitude_label_mask", "0 0 0 1").split()
    lon_label_mask = [int(x) for x in lon_label_mask]
    assert len(lon_label_mask) == 4
    if we_gridline_freq is not None:
        loLon = round(m.llcrnrlon)
        hiLon = round(m.urcrnrlon)
        foo = m.drawmeridians(np.arange(loLon, hiLon, we_gridline_freq),
                    labels=lon_label_mask, color=kwargs["longitude_line_color"])#,
                        #textcolor=kwargs["longitude_label_color"])
        setlabelcolor(foo, kwargs["longitude_label_color"])
        #plt.gca().xaxis.label.set_color(kwargs["longitude_label_color"])
    try:
        country_outlines_color = kwargs["country_outlines_color"]
    except KeyError:
        country_outlines_color = 'k'
    try:
        country_outlines_width = kwargs["country_outlines_width"]
    except KeyError:
        country_outlines_width = 0.5
    if "draw_countries" in kwargs:
        m.drawcountries(linewidth=country_outlines_width, color=country_outlines_color)


def create_simple_legend(ax=None, skip_duplicates=True, alpha=0.9, 
                         position='best'):
    '''
    Creates  a simple legend using the handles and labels in gca() and 
    decorates it based on the cfg options.
    *Therefore, make sure the pass in a 'label' to the plotting function being
     used*
    
    OPTIONAL ARGUMENTS
    
    ax - Axes to show legend on. If not passed in, use gca()
    skip_duplicates - If True and there are elements with duplicate
                      labels plotted, only one of them is included in the 
                      legend. Otherwise, for example, if there are multiple 
                      plots for the same dataset (e.g. for different cycles) 
                      there would  be a bunch of duplicate entries.
    alpha - set alpha level. Default = 0.9
    position - argument passed as 'loc' argument of plt.legend()
    '''
    if ax is None:
        ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    if skip_duplicates:
        unique_handles = []
        unique_labels = []
        for i in range(len(labels)):
            if not labels[i] in unique_labels:
                unique_labels.append(labels[i])
                unique_handles.append(handles[i])
        handles = unique_handles
        labels = unique_labels
    leg = plt.legend(handles, labels, loc=position)
    leg.get_frame().set_alpha(alpha)


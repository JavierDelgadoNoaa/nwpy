'''
This module provides functions for plotting observations on a Basemap
'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca, savefig
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import ColorConverter

#
# Define some defaults
#
WEST_CORNER = -80.0
EAST_CORNER = -30.0
NORTH_CORNER = 40.0
SOUTH_CORNER = 10.0
# how often to put horizontal and vertical grid line, in degrees
GRIDLINE_WE = 10
GRIDLINE_NS = 10
# Coastlines and lakes with an area (in sq. km) smaller than this will not be plotted
AREA_THRESHOLD = 1000.
# Map resolution should be (c)rude, (l)ow, intermediate, high, full, or None
MAP_RESOLUTION = 'l'
# Map projection - only 'cyl' is tested
MAP_PROJECTION = 'cyl'
# Default size of marker used for plotting the obs
MARKER_SIZE = 10


def plot_obs(lats, lons, **kwargs):
    '''
    Given a list of latitudes and a list of longitudes corresponding to
    observation locations, create a scatter plot of the points.

    RETURNS - The basemap object used (if the optional basemap was passed in, 
               same basemap instace is returned)
               
    REQUIRED ARGUMENTS
     - lats - List of latitude points
     - lons - List of longitude points

    OPTIONAL KEYWORD ARGUMENTS (and their corresponding functionality):
     - label - Label to  use for the legend
     - marker - Type of marker to use
     - ax - Axes object to plot on. This overrides the 'basemap' argument.
     - basemap - If set, use it as the Basemap object for plotting.
                 Otherwise instantiate a new one.
     - extents - If set and basemap is not passed in, use these extents
                 for the map. Should be a 4-element list as follows:
                     [lowerLat, westLon, upperLat, eastLon]
                 If not passed in, use the global variables defined in
                 this module
      * If `basemap' and/or `extents' are omitted, global variables will be used
        to determine map extents and other options.
      * If None is passed in for `basemap', a new Basemap will be instantiated
      * If `basemap' is passed in and is not None, `extents' will be ignored
    '''
    # Set basic optional values

    # array of arguments to pass to the plot() call
    mpl_plot_kwargs = {}

    # Set the axes
    if kwargs.has_key('ax'):
        ax = kwargs['ax']
    else:
        ax = gca()

    # Instantiate/set the Basemap object
    if kwargs.has_key('basemap') and kwargs['basemap'] is not None:
        m = kwargs['basemap']
    else:
        if kwargs.has_key('extents'):
            lowerLat = kwargs['extents'][0]
            westLon = kwargs['extents'][1]
            upperLat = kwargs['extents'][2]
            eastLon = kwargs['extents'][3]
        else:
            lowerLat = SOUTH_CORNER
            westLon = WEST_CORNER
            upperLat = NORTH_CORNER
            eastLon = EAST_CORNER
        m = Basemap(llcrnrlon=westLon,llcrnrlat=lowerLat,
                    urcrnrlon=eastLon,urcrnrlat=upperLat,
                    projection=MAP_PROJECTION, resolution=MAP_RESOLUTION,
                    area_thresh=AREA_THRESHOLD,
                    ax=ax)
    if kwargs.has_key('label'):
        label = kwargs['label']
    else:
        label = None
    if kwargs.has_key('facecolor'):
        facecolor=kwargs['facecolor']
    else:
        facecolor = 'red'
    if kwargs.has_key('marker'):
        marker = kwargs['marker']
    else:
        marker = '.'
    #
    # Plot the data
    #
    m.scatter(lons, lats,  latlon=True,
           label=label,
           marker=marker,
           color=facecolor,
           s=MARKER_SIZE,
            )
    #          marker=markers[obGroup])

    return m
  

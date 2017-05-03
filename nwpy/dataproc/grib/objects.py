'''
This module provides classes to facilitate processing of GriB data using PyGrib

Javier.Delgado@noaa.gov
'''
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pygrib
import logging as log


DEFAULT_CONTOUR_LINE_PLOT_COLOR = 'black'
DEFAULT_COLOR_MAP = plt.cm.jet

class GribMessage:
   ''' 
   Very light-weight Encapsulation of a GriB message. Currently, it
   just stores it's name, lats, lons, and values. 
   The PyGrib grib message is optionally stored as well. This is 
   optional since it may take too much memory.
   The class provides methods for plotting the values onto an existing
   Basemap object as a contour or filled contour.
   '''
   def __init__(self, param_name, level, lats, lons, values, plotType, 
                grib_message=None, **kwargs):
      self.param_name = param_name # should match parameterName attribute of grib file
      self.level = level
      self.lats = lats
      self.lons = lons
      self.values = values
      self.plot_type = plotType
      if grib_message != None:
          self.grib_message = grib_message

   def plot_contour_lines(self, basemap, colors=DEFAULT_CONTOUR_LINE_PLOT_COLOR):
      '''
      Make a contour plot from the GriB message's values onto the given 
      Basemap object
      '''
      #small = np.min(self.values)
      #large = np.max(self.values)
      x,y = basemap(self.lons, self.lats)
      # use linspace arg to control how many bins to break the color map into
      #cs = map.contour(x, y, self.values, np.linspace(small,large,11), colors=colors )
      cs = basemap.contour(x, y, self.values, colors=colors ) 
      plt.clabel(cs, inline=1, fontsize=10) # put values in contour lines

   def plot_contours(self, basemap):
      '''
      Make a filled contour plot from the grib message's values onto the
      given basemap object
      '''  
      x,y = basemap(self.lons, self.lats)
      rangeMin = np.min(self.values)
      rangeMax = np.max(self.values)
      #cs = m.contourf(x,y,self.values, np.linspace(rangeMax,rangeMin,11), cmap=plt.cm.jet )
      cs = basemap.contourf(x,y,self.values, cmap=DEFAULT_COLOR_MAP )
      cb = basemap.colorbar()
      

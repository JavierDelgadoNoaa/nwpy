#!/usr/bin/env python
'''
Draw a basic map with some squares
'''

import matplotlib
#matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.patches as patches
from nwpy.viz.map.bling import decorate_map

loLon = 273 ; hiLon = 347
loLat = 3 ; hiLat = 43

PROJECTION = 'cyl'
RESOLUTION = 'l'

m = Basemap(llcrnrlon=loLon, llcrnrlat=loLat,
            urcrnrlon=hiLon, urcrnrlat=hiLat,
            projection=PROJECTION, 
            resolution=RESOLUTION)

decorate_map(m)

#Southwest points
grid_sw = {} ; grid_height = {} ; grid_width = {}
grid_sw['osse_hwrf_d01'] = (291.009, 10.889) 
grid_width['osse_hwrf_d01'] = 24.6
grid_height['osse_hwrf_d01'] = 42.2

ax = plt.gca()
for gridName in ['osse_hwrf_d01']:
    ax.add_patch(patches.Rectangle(grid_sw[gridName], 
                                   grid_width[gridName], 
                                   grid_height[gridName]))
plt.show()

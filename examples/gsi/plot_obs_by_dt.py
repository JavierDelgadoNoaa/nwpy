#!/usr/bin/env python
"""
Plot assimilated conv obs from a GsiRun, separating obs according to their dt

USAGE:
  $0 gsi_run_path [date]
  If `date' is not specified, try to guess from path. 
  If specified, it can be in seconds since epoch or YYYYMMDDhhmm

"""

import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca, savefig
from nwpy.io.metcomp.gsi_objects import GsiRun
from nwpy.viz.map.obs_plotter import plot_obs
from nwpy.viz.map.bling import decorate_map
from nwpy.dateutils import conversions


# Set path to diag file readers
host = os.uname()[1]
if host.startswith("tfe"):
    DIAG_READER_RAD = '/home/Javier.Delgado/apps/gsi/comgsi/3.3/2014/util/Analysis_Utilities/read_diag/read_diag_rad.exe'
    DIAG_READER_CONV = '/home/Javier.Delgado/apps/gsi/comgsi/3.3/2014/util/Analysis_Utilities/read_diag/read_diag_conv.exe'
elif host.startswith("fe"):
    DIAG_READER_RAD = '/home/Javier.Delgado/apps/gsi/comgsi/3.3/util/Analysis_Utilities/read_diag/read_diag_rad.exe'
    DIAG_READER_CONV = '/home/Javier.Delgado/apps/gsi/comgsi/3.3/util/Analysis_Utilities/read_diag/read_diag_conv.exe'
else:
    DIAG_READER_RAD = None
    DIAG_READER_CONV = None

#
# Set options
#
USE_PLAINTEXT_DIAG = False # will be overriden to False if the DIAG_READER_* paths do not exist
# If True, use the diagnostics from the first iteration of GSI 3dvar to plot. 
# If False, use the 2rd iteration (i.e. the analysis diags)
PLOT_GES = False #True 

#
# MAIN
#
if __name__ == '__main__':
    # parse commandline and check paths
    if len(sys.argv) < 2:
        sys.stderr.write("%s" %__doc__)
        sys.exit(1)
    path = sys.argv[1]
    if len(sys.argv) == 3:
        dateStr = sys.argv[2] # 1122897600 # 8/1/ 12z
        if len(dateStr) == 12:
            # reasonably safe to assume user is passing a date string
            date = conversions.yyyymmddHHMM_to_epoch(dateStr)
        else:
            date = float(dateStr)
    else:
       sys.exit(13)
       # TODO : Try to find a file with a date suffix (e.g. stdout.anl.200508011200)
    if not os.path.exists(path):
        sys.stderr.write("Given path does not exist")
        sys.exit(2)

    # set the diag_reader vars depending on whether the path set in global variable exists
    use_plaintext_diag = USE_PLAINTEXT_DIAG
    if not use_plaintext_diag:
        # override decision to use plaintext diag if readers are not present
        if os.path.exists(DIAG_READER_RAD):
            radReader = DIAG_READER_RAD
        else:
            radReader = None # use plaintext
            use_plaintext_diag = True
        if os.path.exists(DIAG_READER_CONV):
            convReader = DIAG_READER_CONV 
        else:
            convReader = None # use plaintext
            use_plaintext_diag = True
        
    gsirun = GsiRun(path, date=date, ob_types=['conv'],
                   diag_reader_rad=radReader,
                   diag_reader_conv=convReader,
                   use_plaintext_diag=use_plaintext_diag)
    
      
    # build dictionary to use as data to plot_obs() - obs will be clustered by dt
    ob_groups = {}
    ob_group_names = ['-3<=dt<=-2', '-2<dt<=-1', '-1<dt<=0', '0<dt<=1', '1<dt<=2', '2<dt<=3']
    for obGrp in ob_group_names:
        ob_groups[obGrp] = { 'lats':[], 'lons':[] }
    if PLOT_GES:
        diags = gsirun.firstguess_diags
    else:
        diags = gsirun.analysis_diags 
    if diags['conv'] is None:
        sys.stderr.write("No obs were assimilated")
        sys.exit(1)
    for ob in diags['conv'].non_skipped_obs: 
        dt = ob.time_delta
        if dt >= -3 and dt <= -2:
            ob_groups['-3<=dt<=-2']['lats'].append(ob.lat)
            ob_groups['-3<=dt<=-2']['lons'].append(ob.lon)
        elif dt > -2 and dt <= -1:
            ob_groups['-2<dt<=-1']['lats'].append(ob.lat)
            ob_groups['-2<dt<=-1']['lons'].append(ob.lon)
        elif dt > -1 and dt <= 0:
            ob_groups['-1<dt<=0']['lats'].append(ob.lat)
            ob_groups['-1<dt<=0']['lons'].append(ob.lon)
        elif dt > 0 and dt <= 1:
            ob_groups['0<dt<=1'] ['lats'].append(ob.lat)
            ob_groups['0<dt<=1'] ['lons'].append(ob.lon)
        elif dt > 1 and dt <= 2:
            ob_groups['1<dt<=2'] ['lats'].append(ob.lat)
            ob_groups['1<dt<=2'] ['lons'].append(ob.lon)
        elif dt > 2 and dt <= 3:
            ob_groups['2<dt<=3'] ['lats'].append(ob.lat)
            ob_groups['2<dt<=3'] ['lons'].append(ob.lon)

    # plot it
    facecolors = ['blue', 'green', 'red', 'purple', 'cyan', 'orange', 'gray']
    markers = ['s', 'o', '8', '*', '+', 'x', '.']
    m = None # This will hold the Basemap object after the first call to plot_obs
    for i,obGrp in enumerate(ob_group_names):
        #import pdb ; pdb.set_trace()
        # Note that this is calling all the decoration stuff every time
        print obGrp, ': ', len(ob_groups[obGrp]['lats'])
        if len(ob_groups[obGrp]['lats']) == 0: 
            continue
        m = plot_obs(ob_groups[obGrp]['lats'], ob_groups[obGrp]['lons'], 
                     label=obGrp,
                     basemap=m,
		     facecolor=facecolors[i],
		     marker=markers[i])
    decorate_map(m, 
		         ns_gridline_freq=10,
                 we_gridline_freq=10)
		 
    leg = plt.legend()
    leg.get_frame().set_alpha(0.85)
    savefig('obs.png')

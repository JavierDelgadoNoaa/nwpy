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

if __name__ == '__main__':
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

    gsirun = GsiRun(path, date=date, ob_types=['conv'],
                   diag_reader_rad=DIAG_READER_RAD,
                   diag_reader_conv=DIAG_READER_CONV)
    for obDiag in gsirun.firstguess_diags.itervalues():
       obDiag.make_plain_text_diag_file()
    if gsirun.is_analysis:
        for obDiag in gsirun.analysis_diags.itervalues():
           obDiag.make_plain_text_diag_file()


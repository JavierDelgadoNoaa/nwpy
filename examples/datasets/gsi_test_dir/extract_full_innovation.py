#!/usr/bin/env python
##
# Read GSI-generated diagnostic files for global o-f and o-a data
# and produce statistics.
# Reads from all GSI_ANALYSIS.<date> directories
# Two output files per cycle are generated in OUTPUT_DIRECTORY, 
#  o_minus_a : o-a stats            o_minus_f : o-f stats
# See A.2 (page 153) of GSI manual.
##

import logging as log
import os
import time
import subprocess

class GSIAnalysis(object):
   '''
   Encapsulates the attributes of a GSI Analysis
   '''
   
   def __init__(self, date):
      ''' Instantiate a `GSIAnalysis` with a given analysis date'''
      self.date = date # seconds since epoch
      self.conv_obs = []
      self.rad_obs = []
      self.oz_obs = []
   
   def append_conv_ob(self, ob):
      self.conv_obs.append(ob)

   @property
   def non_skipped_conv_obs(self):
      return [ob for ob in self.conv_obs if ob.skipped == False]
   @property
   def non_skipped_temp_obs(self):
      return [ob for ob in self.conv_obs if ob.skipped == False if ob.type == 't' ]
   @property
   def non_skipped_wind_obs(self):
      return [ob for ob in self.conv_obs if ob.skipped == False if ob.type == 'uv' ]
   @property
   def non_skipped_pres_obs(self):
      return [ob for ob in self.conv_obs if ob.skipped == False if ob.type == 'ps' ]
   @property
   def non_skipped_gps_obs(self):
      return [ob for ob in self.conv_obs if ob.skipped == False if ob.type == 'gps' ]
   @property
   def non_skipped_moisture_obs(self):
      return [ob for ob in self.conv_obs if ob.skipped == False if ob.type == 'q' ]
   

class GSIAbstractObservation(object):
   def __init__(self, type, station_id, time_delta, lat, lon, pres, iuse):
      self.type = type
      self.station_id = station_id
      self.time_delta = time_delta
      self.lat = lat
      self.lon = lon
      self.pres = pres
      if iuse < 1:
         self.skipped = True
      else:
         self.skipped = False

class GSIConventionalObservation(GSIAbstractObservation):
   def __init__(self, type, station_id, time_delta, lat, lon, pres, skipped, value, bias):
      super(GSIConventionalObservation, self).__init__(type, station_id, time_delta, lat, lon, pres, skipped)
      self.value = value
      self.bias = bias # o-g

class GSIWindObservation(GSIAbstractObservation):
   def __init__(self, type, station_id, time_delta, lat, lon, pres, skipped, u_value, v_value, u_bias, v_bias):
      super(GSIWindObservation, self).__init__(type, station_id, time_delta, lat, lon, pres, skipped)
      self.u_value = u_value ; self.v_value = v_value
      self.u_bias = u_bias ; self.v_bias = v_bias




TOP_DIR = os.curdir + "/../" # '/home/Javier.Delgado/projects/osse/model_rundir/new_hyperspectral/conv_obs_only/run'
GSI_DIR = '~/projects/apps/gsi/comgsi/3.3'
GSI_WORKDIR_PREFIX = 'GSI_ANALYSIS.' 

CONV_DIAG_EXTRACTOR = os.path.join(GSI_DIR, 'util/Analysis_Utilities/read_diag/read_diag_conv.exe')
RAD_DIAG_EXTRACTOR = os.path.join(GSI_DIR, 'util/Analysis_Utilities/read_diag/read_diag_rad.exe')

# output data files will be placed here
OUTPUT_DIRECTORY = os.path.join( os.getcwd() , 'postproc', 'gsi_stats')

start_time = time.strptime('8/1/2005 0600', '%m/%d/%Y %H%M')
end_time = time.strptime('8/1/2005 0600', '%m/%d/%Y %H%M') 
frequency = 3600 * 6

start_date = int( time.mktime(start_time) )
end_date = int( time.mktime(end_time) )
print start_date, end_date

def get_date_string(cycle_date):
   return time.strftime('%Y_%m_%d_%H_%M', time.gmtime(cycle_date) )

def get_gsi_dir(cycle_date):
  datestr = get_date_string(cycle_date)
  return GSI_WORKDIR_PREFIX + datestr

def get_diag_file_suffix(cycle_date):
  ''' 
  TODO : look for diag files that don't have the minutes
  '''
  #datestr = time.strftime('%Y%m%d%H', time.gmtime(cycle_date))
  #datestr = time.strftime('%Y%m%d%H', time.localtime(cycle_date))
  datestr = time.strftime('%Y%m%d%H%M', time.localtime(cycle_date))
  return datestr

def create_namelist(suffix, in_file_name, out_file_name):
  elfilename = 'namelist.' + suffix
  elfile = open(elfilename, 'w')
  elfile.write('&iosetup\n')
  elfile.write("  infilename='%s',\n" %in_file_name)
  elfile.write("  outfilename='%s',\n" %out_file_name)
  elfile.write("/\n")
  elfile.close()  

def process_data(infile, anl_date):
  '''
  Process the diag file. Assume the following write statements 
  were used (except the width; just tokenize)
  for non-uv obs:
     write (42,'(A3," @ ",A8," : ",I3,F10.2,F8.2,F8.2,F8.2,I5,2F10.2)') &
            var,stationID,itype,rdhr,rlat,rlon,rprs,iuse,robs1,ddiff
  for uv-obs:
      write (42,'(A3," @ ",A8," : ",I3,F10.2,F8.2,F8.2,F8.2,I5,4F10.2)') &
             var,stationID,itype,rdhr,rlat,rlon,rprs,iuse,robs1,ddiff,robs2, rdpt2
  * Note the '@' and ':' will affect indexing
  ddiff is the o-g  (ov u for uv obs)
  rlat/rlon are obs lat/lon ; rprs is ob pres
  rdhr is the time difference in hours, relative to analysis time
  iuse: was it used
  robs1 - ob value (of u in the case of uv)
  robs2 - ob value (of v)
  rdpt2 is the o-g of v
  UNITS: K,m/s, hPa, 
  
  RETURN : A GSIAnalysis object with the observations populated

  ***  NOTE: It's good to verify these against readconvobs.f90 (in enkf)  and $GSI_ROOT/util/Analysis_Utilities/read_diag/read_diag_conv
             if using a new version of GSI ***
  '''
  obs_vals = []
  diffs = []
  gsi_anl = GSIAnalysis(date=anl_date)
  elfile = open(infile)
  for line in elfile:
     data = line.strip().split()
     if len(data) < 10:
        print 'Unrecognized data:', line
        continue
     station_id = data[2]
     type = data[0]
     subtype = data[4]
     time_delta = float(data[5]) # in hours relative to analysis time
     lat = float(data[6])
     lon = float(data[7] )
     pres = float(data[8])
     iuse = int(data[9])
     value = float(data[10])
     diff = float(data[11]) # this is the o-g value
     if type != 'uv':
        ob = GSIConventionalObservation(type, station_id, time_delta, lat, lon, pres, iuse, value, diff)
     else:
        u = float(data[10])
        v = float(data[12])
        u_diff = float(data[11])
        v_diff = float(data[13])
        ob = GSIWindObservation(type, station_id, time_delta, lat, lon, pres, iuse, u, v, u_diff, v_diff)
     gsi_anl.append_conv_ob(ob)
  return gsi_anl

# __Main__
if not os.path.exists ( OUTPUT_DIRECTORY) : 
   os.makedirs( OUTPUT_DIRECTORY ) 
for cycleDate in range(start_date , end_date+1, frequency):
  print cycleDate
  gsi_dir = get_gsi_dir(cycleDate)
  datestr = get_date_string(cycleDate)
  curr_directory = os.path.join(TOP_DIR, gsi_dir)
  diag_file_suffix = get_diag_file_suffix(cycleDate)  
  os.chdir(curr_directory)
  # GSI will prodcue files with names like diag_conv_ges.<date> and diag_amsua_n15_anl.<date>
  for diagStep in [ 'anl', 'ges' ] :
    in_file_name = 'diag_conv_' + diagStep + '.' + diag_file_suffix
    out_file_name = 'results_conv_' + diagStep + '.' + diag_file_suffix
    create_namelist('conv', in_file_name, out_file_name)
    ret = subprocess.call([CONV_DIAG_EXTRACTOR], shell=True, stdout=open("/dev/null", 'w') )
    # TODO  ensure it returns with exit code 9999    AND    have option of sending output to a file instead of /dev/null
    gsi_anl = process_data(out_file_name, cycleDate)
    # output to a file
    if diagStep == 'anl': outfile = open( os.path.join(OUTPUT_DIRECTORY, 'o_minus_a.' + datestr + '.dat'), 'w')
    else: outfile = open( os.path.join(OUTPUT_DIRECTORY, 'o_minus_f.' + datestr + '.dat'), 'w')
    for ob in gsi_anl.non_skipped_conv_obs:
       if isinstance(ob, GSIWindObservation):
         outfile.write('%s, %f, %f, %f, %f, %f, %f, %f, %f\n' %(ob.type, ob.time_delta, ob.lat, ob.lon, ob.pres, ob.u_value, ob.v_value, ob.u_bias, ob.v_bias) )
       else:
         outfile.write('%s, %f, %f, %f, %f, %f, %f\n' %(ob.type, ob.time_delta, ob.lat, ob.lon, ob.pres, ob.value, ob.bias) )
    outfile.close()


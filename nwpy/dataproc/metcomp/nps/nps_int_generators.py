"""
This module contains functions to facilitate creating derived variables in
NPS intermediate format.
"""
import os
import subprocess
import logging
import shutil

from nps_int_utils import get_int_file_name

def concatenate_files(src_file_name, dest_file_name) :
   '''
   Appends contents of <src file name> to <dest file name>
   '''
   output_file = open( dest_file_name , 'ab' )
   source_file = open( src_file_name, 'rb' )
   shutil.copyfileobj( source_file, output_file, 65536 )
   source_file.close()
   output_file.close()

def create_ght_geos2wrf(kwargs):
    #inPrefix, outPath, currDate, createHGTexePath, inDir=".", modelTop=1.0):
    """
    Create geopotential height using using the createGHT tool of geos2wrf.
    Creates and runs in temp directory <inDir>.<currDate>, so that multiple 
    dates can be processed simultaneously.
    In the end, the nps_int file will either be copied to "inDir" or merged
    with the input file, depending on ``catOutput'' setting
    args expected in kwargs:
    @param inPrefix Prefix of file containing input data (part before the date)
    @param outPath Path to write the field to
    @param currDate datetime object with date being processed
    @param geos2wrf_utils Path in which to find the "createHGT" executable
    @param catOutput If True, concatenate the generated file (SOILHGT:<date>) 
           into the input file (<inPrefix>:<date>)
    @param inDir (optional) Directory containing input file. Default: '.'
    @param modelTop Pressure of model top in mb. Default: 1.0)
    @param log (optional) Logger to log with. If None, no logging is done
    """
    
    # required args
    inPrefix = kwargs["inPrefix"]
    outPath = kwargs["outPath"]
    currDate = kwargs["currDate"]
    geos2wrf_utils = kwargs["geos2wrf_utils"]
    createHGTexePath = os.path.join(geos2wrf_utils, "createHGT")
    assert os.path.exists(createHGTexePath)
    #optional args
    try:
        inDir = kwargs["inDir"]
    except KeyError:
        inDir = "."
    try:
        modelTop = float(kwargs["modelTop"])
    except KeyError:
        modelTop = 1.0
    try:
        log = kwargs["log"]
    except KeyError:
        log = None
    try:
        catOutput = kwargs["catOutput"]
    except KeyError:
        catOutput = False
    
    # Set variables
    input_file = get_int_file_name(inPrefix, currDate, includeMMSS=False)
    output_file = get_int_file_name("HGT_MODEL_LEVEL", currDate)
    
    # Create temporary run directory, since there may be multiple 
    # parallel executions
    #import pdb ; pdb.set_trace()
    tmpDir = os.path.join(inDir, "run.{d:%Y%m%d%H%M}".format(d=currDate))
    os.mkdir(tmpDir)
    orig_dir = os.getcwd()
    os.chdir(tmpDir)
    os.symlink( os.path.join(inDir, input_file), input_file)
    # Create namelist
    namelist = open("namelist.createHGT", "w")
    txt = []
    txt.append("&input")
    txt.append("directory='{indir}',".format(indir=tmpDir))
    txt.append("prefix={p}".format(p=inPrefix))
    txt.append("{d:year=%Y,\nmonth=%m,\nday=%d,\nhour=%H,}".format(d=currDate))
    txt.append("layerPressureThicknessName='DELP',")
    txt.append("layerTemperatureName='TT',")
    txt.append("layerSpecificHumidityName='SPECHUMD',")
    txt.append("soilHeightName='SOILHGT',")
    txt.append("modelTopPressure={f},".format(f=modelTop))
    txt.append("/%")
    txt = [ s+"\n" for s in txt]
    if  log:
        #fo = "\n".join(txt)
        log.debug("namelist.createHGT contents: {0}".format(txt))
    namelist.writelines(txt)
    namelist.close()    
    cmd = [createHGTexePath]
    
    # Run it
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output = p.communicate()[0].strip()
    if log:
        log.info("Output from createHGT: {0}".format(output))
        log.debug("Assuming units of 'm'. Could change in future geos2wrf")
    units="m"
    description = "Geopotential height. Computed with geos2wrf"
    outArgs = {"units": units, "description":description}

    # Merge output to main file or move to input file's directory
    if catOutput:
        [src,dest] = [os.path.join(tmpDir, p) for p in [output_file,input_file]]
        concatenate_files(src, dest)
        os.unlink(output_file)
    else:
        os.move(os.path.join(tmpDir, output_file), inDir)
    for fil in ["namelist.createHGT", input_file]:
        os.unlink(fil)
    os.chdir(orig_dir)
    os.rmdir(tmpDir)
    return outArgs

def test_ght():
    from datetime import datetime as dtime
    infile = "/home/Javier.Delgado/scratch/nems/g5nr/data/catera_pl/nps_int/G5NR:2006-09-10_00"
    in_dir = "/home/Javier.Delgado/scratch/nems/g5nr/data/catera_pl/nps_int"
    geos2wrf_path = "/home/Javier.Delgado/scratch/apps_tmp/nuwrf/dist/nu-wrf_v7lis7-3.5.1-p6/utils/geos2wrf_2"
    tim = dtime(year=2006, month=9, day=10, hour=0)
    log = logging.getLogger()
    #handler = logging.Handler(logging.DEBUG)
    handler = logging.Handler()
    ch = logging.StreamHandler()
    log.addHandler(ch)
    log.setLevel(logging.DEBUG)
    kwargs = {"inPrefix":"G5NR", "outPath":".", "currDate": tim,
              "geos2wrf_utils":geos2wrf_path, "inDir":in_dir, "log": log,
              "catOutput":True }
    create_ght_geos2wrf(kwargs)

if __name__ == "__main__":
    test_ght()

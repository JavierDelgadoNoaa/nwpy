import os
import sys
from optparse import OptionParser
import logging
from setuptools import find_packages

from numpy.distutils.core import Extension

#
# DEFAULTS
##
# Do not try to compile any natively-compiled stuff
g_skip_native = False
# Skip compilation of native-compiled fortran interpolation codes
g_skip_interp = False
# Skip compilation of native-compiled NPS/WPS intermediate format
# file utilities
g_skip_nps_int_native = False

def _parse_args():
    """
    Process command line arguments specific to nwpy, which  involves
    setting global fields and removing command line options so they
    can be passed to distutils::setup()
    """
    global g_skip_native, g_do_interp
    '''
    #import pdb ;pdb.set_trace()
    #usage="usage: %prog [options] <database file> <config file1 [config_file2 ...]> "\
    #      "[config_fileN ... config_fileN+M] [input_file1,...input_fileN]"\
    #      "".format() + __doc__
    usage = "TODO"
    parser = OptionParser(usage=usage)
    #parser.add_option("-f", "--file", dest="input_file", optional=False)
    #parser.add_option("-c", "--config", dest="config_file")
    parser.add_option("-n", "--skip-native", dest="skip_native",
                      action="store_true", default=False,
                      help="Compile all native-compiled things"
                      )
    parser.add_option("-l", "--log-level", dest="log_level", default=logging.INFO,
                      help="(0-100 or constant defined in the logging module")
    parser.add_option("-o", "--conf-override", dest="conf_overrides",
                      action="append",
                      help="Override config setting (section.item=value)")
    (options, args) = parser.parse_args()
    g_skip_native = False
    if options.skip_native: g_skip_native = True
    '''
    copy_args=sys.argv[1:]
    if '--skip-native' in copy_args:
      g_skip_native = True
      copy_args.remove('--skip-native')
    if "--skip-interp" in copy_args:
        g_skip_interp = True
        copy_args.remove("--skip-interp")
    if "--skip-nps-int-utils" in copy_args:
        g_skip_nps_int_native = True
        copy_args.remove("--skip-nps-int-utils")
    return copy_args

if __name__ == "__main__":
    from numpy.distutils.core import setup
    ext_modules = []
    copy_args = _parse_args()
    #f2py_opts = dict(f90exec="ifort")
    #f2py_opts = ["-m", "-c", "--f90exec", "ifort"]
    if not g_skip_native:
        # if g_do_interp:
            ext2 = Extension(name = 'interp_routines',
                             sources = ['src/interp/vinterp/interp_routines.f90'],
    #                        f2py_options = f2py_opts,
                    )
            ext_modules.append(ext2)
        # if g_do_nps_int
            ext1 = Extension(name="write_nps_int",
                             sources = ["src/io/metcomp/nps/write_nps_int.f90"],
                            )
            ext_modules.append(ext1)
    else:
        print "Skipping all modules that require native compile"

    os.environ["FFLAGS"] = "-convert big_endian -fPIC"
    setup(
          name              = "nwpy",
          ext_modules=ext_modules,
          script_args=copy_args,
          version           = "0.2",
          description       = "Plot data on map projections with matplotlib",
          long_description  = """TODO""",
          author            = "Javier Delgado",
          author_email      = "Javier.Delgado@noaa.gov",
          platforms         = ["any"],
          keywords          = ["python","numerical weather prediction", "interpolation", 
                               "io"],
                              # https://pypi.python.org/pypi?%3Aaction=list_classifiers
          classifiers       = ["Development Status :: 3 - Alpha",
                               "Intended Audience :: Science/Research",
                               "Programming Language :: Python",
                               "Topic :: Scientific/Engineering :: Visualization",
                               "Topic :: Scientific/Engineering :: Atmospheric Science",
                               "Operating System :: OS Independent"],
         # By specifying packages=find_packages(), it will include the .py files corresponding
         # to the .so files and create a .pth file in the install_prefix
         # ** REQUIRES you to have --prefix path in PYTHONPATH
         # ** Will __not__ autoinclude any other .py files
         # ...BUT...if there is a a directory matching "name" above it will
         # include the python files there and any subdirectories containing python files
         #
         # "If you do supply your own MANIFEST, you must specify everything: the default set of files described above does not apply in this case."
         # -> but not if the dir name matches project name,apparently
         # -> since this is from setuptools, maybe it doesn't even use the manifest
         #     (or maybe it does, just not for src code)
         packages           = find_packages(exclude=["era*"]), # recursively includes all subdirs with an __init__.py
         # ** setup.py install will not include package_data in the MANIFEST.in
         #    (setup.py bdist will), so non-source files must be included here
         package_data = { "nwpy.io.metcomp.nps":["etc/*.txt"] },
         #recursive-include nwpy/io/metcomp/nps/lib/etc *txt
         # https://github.com/alphaomega-technology/Equation
         depends=["Equation"],
         include_package_data = True,
         )
         # ** UPDATE: (TODO)**
         # According to https://python-packaging.readthedocs.io/en/latest/non-code-files.html,
         # gotta set include_package_data = True to include data files from MAnifest
         # Also, in Manifest, to include individual subdirs: include docs/*.txt
         # Also, see NOTE in the above link - data files should be in module directory (think nps default_units)

#f2py -m -c --f90exec=`which ifort` interp_routines interp_routines.f90

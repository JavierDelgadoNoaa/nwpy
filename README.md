# TODO
-Add README to Manifest (or use autoincluded file prefix)
-Fix issue with WrfRun (see hwrf screen on Jet)

# nwpy


To build the fortan modules  with f2py:
python setup.py build --fcompiler=intelem
* This will put natively-compiled libraries in a temporary build directory.
  If not doing a system install, you'll need to move them to the 
  corresponding lib folder

To do a system install:
python setup.py install [--prefix=/path/to/topdir]
* libs will go in: /path/to/topdir/lib/python2.7/site-packages   


# For SpecData
The stock input specification ("inspec") files are located in nwpy/io/specdata/conf
 -> Copy these somewhere to use them in apps that use the specdata library


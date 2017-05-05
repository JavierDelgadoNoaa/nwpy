# nwpy:
NWPY_DIST_PATH=$(readlink -f `dirname $0`/..) # path before the etc/ directory
export PYTHONPATH=$NWPY_DIST_PATH:$PYTHONPATH

# for Equation - may be used by specdata
echo $PYTHONPATH | grep ":/home/Javier.Delgado/local/lib/python2.7/site-packages:" &> /dev/null
[[ $? != 0 ]] && export PYTHONPATH=$PYTHONPATH:/home/Javier.Delgado/local/lib/python2.7/site-packages

scr=`readlink -f $0`
python -c "import nwpy"
[[ $? != 0 ]] && echo "nwpy import failed. Check $scr" && return 1

return 0

# nwpy:
export PYTHONPATH=$PYTHONPATH:/home/Javier.Delgado/libs/nwpy/nwpy

# for Equation - may be used by specdata
echo $PYTHONPATH | grep ":/home/Javier.Delgado/local/lib/python2.7/site-packages:" &> /dev/null
[[ $? != 0 ]] && export PYTHONPATH=$PYTHONPATH:/home/Javier.Delgado/local/lib/python2.7/site-packages


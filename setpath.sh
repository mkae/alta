# Bash command to set ALTA environment for ease of use.
# Invoke this script one GNU/Linux or OSX using the source
# command:
#
#    source setpath.sh
#
# ALTA command line programs and plugins will be available
# to the shell.
#
path=`pwd`/sources/build

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$path
export PATH=$PATH:$path
export PYTHONPATH=$PYTHONPATH:$path

# Bash command to set ALTA environment for ease of use.
# Invoke this script one GNU/Linux or OSX using the source
# command:
#
#    source setpath.sh
#
# ALTA command line programs and plugins will be available
# to the shell.
#
sources=`pwd`
if [ "$BASH_VERSION" ]; then
	sources=`dirname $BASH_SOURCE`
	sources=`builtin cd ${sources}; builtin pwd`
elif [ "$ZSH_VERSION" ]; then
	sources=$(dirname "$0:A")
fi
path=${sources}/sources/build/
scripts=${sources}/sources/scripts/

export ALTA=$sources
export ALTA_LIB=$path
export ALTA_DIR=${sources}/sources/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$path
export PATH=$PATH:$path:$scripts
export PYTHONPATH=$PYTHONPATH:$path

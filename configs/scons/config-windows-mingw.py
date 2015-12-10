import os, sys

##----------------------------------------------------------------##
## This file describes required and optional arguments to ALTA    ##
## compilation. If you want to manualy tune the use of an already ##
## present library, uncomment the according lines.                ##
##----------------------------------------------------------------##


## Compilators build flags
##
CXX            = 'g++'

# This is not the best place to past static linking since plugins need the -shared options 
#LINKFLAGS			 = [ '-static-libgcc', '-static-libstdc++']

## ALTA internal flags
##
CORE_LIB       = []
SOFT_LIB       = ['core']
PLUGIN_LIB     = ['core']


## Python and boost-python library
##
PYTHON_INC    = ['']
PYTHON_DIR    = []


## Eigen library
##
EIGEN_INC     = ['#external/build/include/Eigen']


## OpenMP flags
##
OPENMP_LIB    = ['gomp']


## OpenEXR library
##
OPENEXR_INC    = ['']
OPENEXR_DIR    = ['']


## FLANN library
##
FLANN_INC    = ['']
FLANN_DIR    = ['']
FLANN_LIB    = ['flann']


## QUADPROG library
##
## You have to specify the directory of the QuadProg library
##
QUADPROG_INC      = ['#external/quadprog++']
QUADPROG_DIR      = ['#external/build/lib']
QUADPROG_LIB      = ['quadprog++']


## CERES library
##
## You have to specify both the directory of the CERES library
## and the glog library
##
CERES_INC      = ['#external/build/include']
CERES_DIR      = ['#external/build/lib']
CERES_LIB      = ['ceres', 'miniglog']
CERES_OPT_LIB  = ['gomp', 'lapack', 'blas']


## NlOpt library
##
## You have to specify the directory of the NlOpt library
##
NLOPT_INC      = ['#external/build/include']
NLOPT_DIR      = ['#external/build/lib']
NLOPT_LIB      = ['nlopt']
NLOPT_OPT_LIB  = []


## coin IpOpt library
##
## You have to specify the directory of the IpOpt library
##
IPOPT_INC      = ['#external/build/include']
IPOPT_DIR      = ['#external/build/lib']
IPOPT_LIB      = ['ipopt']
IPOPT_OPT_LIB  = []


## MATLAB library and Engine
##
#MATLAB_INC  = ['']
#MATLAB_DIR  = ['']
#MATLAB_LIB  = ['eng', 'mex','mat']

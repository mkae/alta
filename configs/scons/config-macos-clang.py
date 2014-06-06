import os, sys

##----------------------------------------------------------------##
## This file describes required and optional arguments to ALTA    ##
## compilation. If you want to manualy tune the use of an already ##
## present library, update the according lines                    ##
##                                                                ##
## This configuration file is made of a MacOS X operating system  ##
## version 10.9, with clang and OpenEXR and FLANN libraries       ##
## installed using Ports.                                         ##
##----------------------------------------------------------------##


## Compilators build flags
##
CXX            = 'clang'
CCFLAGS        = ['-O3', '-Wall', '-m64']
LINKFLAGS      = []


## ALTA internal flags
##
CORE_LIB       = ['dl', 'stdc++']
SOFT_LIB       = ['core', 'dl', 'stdc++']
PLUGIN_LIB     = ['core']


## OpenMP flags
##
OPENMP_FLAG   = []
OPENMP_LIB    = []


## OpenEXR library
##
OPENEXR_INC    = ['/opt/local/include/OpenEXR']
OPENEXR_DIR    = ['/opt/local/lib']
OPENEXR_LIB    = ['Half', 'IlmImf', 'IlmThread']


## FLANN library
##
FLANN_INC    = ['/opt/local/include']
FLANN_DIR    = ['/opt/local/lib']
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
CERES_LIB      = ['ceres', 'glog']
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
MATLAB_INC  = ['/Applications/MATLAB_R2014a.app/extern/include/']
MATLAB_DIR  = ['/Applications/MATLAB_R2014a.app/bin/maci64/']
MATLAB_LIB  = ['eng', 'mex','mat']

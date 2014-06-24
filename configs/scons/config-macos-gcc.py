import os, sys

##----------------------------------------------------------------##
## This file describes required and optional arguments to ALTA    ##
## compilation. If you want to manualy tune the use of an already ##
## present library, uncomment the according lines.                ##
##----------------------------------------------------------------##


## Compilators build flags
##
CXX            = 'gcc'
CCFLAGS        = ['-O3', '-Wall', '-m64']
LINKFLAGS      = []


## ALTA internal flags
##
CORE_LIB       = ['dl']
SOFT_LIB       = ['core', 'dl']
PLUGIN_LIB     = ['core']


## OpenMP flags
##
OPENMP_FLAG   = ['-fopenmp']
OPENMP_LIB    = ['gomp']


## OpenEXR library
##
OPENEXR_INC    = ['/usr/include/OpenEXR']
OPENEXR_DIR    = ['/usr/lib']
OPENEXR_LIB    = ['Half', 'IlmImf', 'IlmThread']


## FLANN library
##
FLANN_INC    = ['/usr/include/flann']
FLANN_DIR    = ['/usr/lib/x86_64_linux-gnu']
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
CERES_OPT_LIB  = ['gomp', 'lapack', 'blas', 'amd', 'camd', 'ccolamd', 'colamd', 'cholmod', 'cxsparse']


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
MATLAB_LIB  = ['eng', 'mx','mat']

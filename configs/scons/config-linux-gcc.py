import os, sys

##----------------------------------------------------------------##
## This file describes required and optional arguments to ALTA    ##
## compilation. If you want to manualy tune the use of an already ##
## present library, uncomment the according lines.                ##
##----------------------------------------------------------------##


## Compilators build flags
##
CXX            = 'g++'


## ALTA internal flags
##
CORE_LIB       = ['rt', 'dl']
SOFT_LIB       = ['core', 'rt', 'dl']
PLUGIN_LIB     = ['core']


## Python and boost-python library
##
PYTHON_INC    = ['/usr/include/python2.7']
PYTHON_DIR    = []
PYTHON_LIB    = ['boost_python', 'python2.7']


## Eigen library
##
EIGEN_INC     = ['#external/build/include/Eigen']


## OpenMP flags
##
OPENMP_FLAG  = ['-fopenmp']
OPENMP_LIB   = ['gomp']


## OpenEXR library
##
OPENEXR_INC    = ['/usr/include/OpenEXR']
OPENEXR_DIR    = ['/usr/lib/x86_64_linux-gnu']
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
CERES_LIB      = ['ceres', 'miniglog']
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
MATLAB_INC  = ['']
MATLAB_DIR  = ['']
MATLAB_LIB  = ['eng', 'mex','mat']

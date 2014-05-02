import os, sys

##----------------------------------------------------------------##
## This file describes required and optional arguments to ALTA    ##
## compilation. If you want to manualy tune the use of an already ##
## present library, uncomment the according lines.                ##
##----------------------------------------------------------------##


## Compilators build flags
##
CXX            = 'g++'
CCFLAGS        = ['-O3', '-Wall', '-m64']
LINKFLAGS      = []


## OpenMP flags
##
OPENMP_FLAGS  = ['-fopenmp']
OPENMP_LIBS   = ['gomp']


## OpenEXR library
##
#OPENEXR_INC    = ['/usr/include/OpenEXR']
#OPENEXR_DIR    = ['/usr/lib']
#OPENEXR_LIB    = ['Half', 'IlmImf', 'IlmThread']


## QUADPROG library
##
## You have to specify the directory of the QuadProg library
##
QUADPROG_INC      = ['#external/quadprog++']
QUADPROG_DIR      = ['#external/build/lib']
QUADPROG_LIBS     = ['quadprog++']


## CERES library
##
## You have to specify both the directory of the CERES library
## and the glog library
##
CERES_INC      = ['#external/build/include']
CERES_DIR      = ['#external/build/lib']
CERES_LIBS     = ['ceres', 'glog']
CERES_OPT_LIBS = ['gomp', 'lapack', 'blas']


## NlOpt library
##
## You have to specify the directory of the NlOpt library
##
NLOPT_INC      = ['#external/build/include']
NLOPT_DIR      = ['#external/build/lib']
NLOPT_LIBS     = ['nlopt']
NLOPT_OPT_LIBS = []


## coin IpOpt library
##
## You have to specify the directory of the IpOpt library
##
IPOPT_INC      = ['#external/build/include']
IPOPT_DIR      = ['#external/build/lib']
IPOPT_LIBS     = ['ipopt']
IPOPT_OPT_LIBS = []


## MATLAB library and Engine
##
MATLAB_INC  = ['/Applications/MATLAB_R2014a.app/extern/include/']
MATLAB_DIR  = ['/Applications/MATLAB_R2014a.app/bin/maci64/']
MATLAB_LIBS = ['eng', 'mex','mat']

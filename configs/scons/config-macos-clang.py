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
CXX            = 'clang++'
CCFLAGS        = ['-O3', '-g', '-Wall', '-Xarch_x86_64', '-mmacosx-version-min=10.9']
LINKFLAGS      = ['-headerpad_max_install_names', '-Xarch_x86_64', '-mmacosx-version-min=10.9']


## ALTA internal flags
##
PLUGIN_LIB     = ['core', 'stdc++']


## Python and boost-python library
## Be sure to install boost-python and python
## You can do this by using Macports (http://www.macports.org/) 
## or Homebrew (http://brew.sh/)
PYTHON_INC    = ['/opt/local/include', '/usr/include/python2.7']
PYTHON_DIR    = ['/opt/local/lib']
PYTHON_LIB    = ['boost_python-mt', 'python2.7']


## Eigen library
##
EIGEN_INC     = ['#external/build/include/Eigen']


## OpenMP flags
##
OPENMP_FLAG   = []
OPENMP_LIB    = []


## OpenEXR library
##
#OPENEXR_INC    = ['/opt/local/include/OpenEXR']
#OPENEXR_DIR    = ['/opt/local/lib']


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
CERES_LIB      = ['ceres', 'miniglog']
CERES_OPT_LIB  = ['lapack', 'blas', 'amd', 'camd', 'ccolamd', 'colamd', 'cholmod', 'cxsparse']


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

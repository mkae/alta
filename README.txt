This is the repository for the ALTA library. ALTA is a set of cross-platform
command line executables and shared object plugins allowing to analyze, fit 
and understand BRDF data and models.

ALTA targets: people working on BRDFs and willing to benchmark new BRDF models
and compare them with state-of-the-art BRDF models and  data easily; people 
working on optical measurements and wanting to experiment different fitting 
procedures and models; or people wanting to perform statistical analysis on
your BRDF data.


1. Repository Organization
==========================

data/:      contains the data for which the fitting techniques are tested
            against. Files are separated by dimension of the input domain
            (e.g. 1d, 2d, 3d, ...).
configs/:   contains platform dependent configuration files for qmake and
            scons.
documents/: contains the documentation, should be build using doxygen into
            that directiory.
external/:  contains any third party library that needs to be used by ALTA.
            Contains a modified QuadProg++ library using Eigen. You can
            download and compile automaticaly some third party libraries using
            provided python scripts obtain_[libname].py.
sources/:   contains all the source files. A Makefile or VS project can be
            created there from the .pro file. Scons generation of the
            project is also supported.



2. Building
===========

ALTA uses SCons, a Python-based build system:

  http://scons.org/

To build ALTA, run a command like the following from the top-level
source directory (noted ${ALTA} in the remainder of this document):

  $ scons --cfg=config-file.py

Here, 'config-file.py' must be replaced with a suitable configuration
file for your platform.  For instance, when building with GCC on
GNU/Linux, you may run:

  $ scons --cfg=configs/scons/config-linux-gcc.py


3. Building Advises (for Qt enthousiasts)
==========================================

We use heavily the Qt profile functionality. To build some of the plugins you
will be required to create your own system dependant .prf for any used library.
For example ALTA core use the Eigen library. Therefore it is mandatory that
you provide a eigen.prf file and that this file is in your QMAKEFEATURES
directory. We provide example of such file in the ${ALTA}/configs/qt directory.
To ease the configuration for a first installation of ALTA, we advise to use
this directory as the main QMAKEFEATURES directory.


3.1. Dependencies:
------------------

 ALTA core: Eigen
 Plugin rational_eigen:     Eigen 3.x
 Plugin rational_quadprog:  Quadprog++
 Plugin rational_cgal:      The CGAL library
 Plugin rational_parallel:  The OpenMP library, Quadprog++ library and Eigen
 Plugin rational_matlab:    Matlab engine (matlab.prf required)
 Plugin rational_parsec_*:  PLASMA coreblas, PaRSEC runtime
 Plugin nonlinear_eigen:    Eigen
 Plugin nonlinear_ceres:    CERES library and its dependencies
 Plugin nonlinear_nlopt:    NLOpt library and its dependencies
 Plugin nonlinear_ipopt:    IpOpt library and its dependencies


3.2. Eigen Plugins
-------------------
 You must provide an eigen.prf file that contains

 INCLUDEPATH *= PATH_TO_EIGEN_DIRECTORY


3.3. Quadprog++ Plugins
-----------------------

 We provide our own version of quadprog++ which uses Eigen library.
 To compile it:

 Go to external/quadprog++/
 Use qmake (to generale the Makefile)
 make

 Then create a quadprog.prf file that will include
 LIBS *= PATH_TO_LIBQUADPROG/libquadprog++.a
 QMAKE_LIBDIR *= PATH_TO_LIBQUADPROG
 INCLUDEPATH *= PATH_TO_LIBQUADPROG_HEADERS


3.4. Parallel Plugin that requires OpenMP
-----------------------------------------

 Create an  openmp.prf file and add the following directives:
 QMAKE_CXXFLAGS *=-fopenmp


3.5. Matlab Plugin
------------------

 Create an matlab.prf file and add the following directives:
 The PATH_TO_MATLAB_INCLUDE_DIRECTORY must point to a directory that
 contains the file engine.h

 INCLUDEPATH *= PATH_TO_MATLAB_INCLUDE_DIRECTORY



4. Generate the documentation using Doxygen
===========================================

  $ cd ${ALTA}/documents/
  $ doxygen doxygen.conf

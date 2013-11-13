This is the repository for the ALTA Project

Repository Organization

data/:      contains the data for which the fitting techniques are tested 
            against. Files are separated by dimension of the input domain
            (e.g. 1d, 2d, 3d, ...).
documents/: contains the documentation, should be build using doxygen into 
            that directiory.
external/:  contains any third party library that needs to be used by ALTA.
            Contains a modified QuadProg++ library using Eigen. You can
            download and compile automaticaly some third party libraries using 
            the Scons scripts.
sources/:   contains all the source files. A Makefile or VS project can be 
            created there from the .pro file. Scons generation of the
            project is also supported.



Building Advises (for Qt enthousiasts)

We use heavily the Qt profile functionality. To build some of the plugins you 
will be required to create your own system dependant .prf for any used library.
For example ALTA core use the Eigen library. Therefore it is mandatory that 
you provide a eigen.prf file and that this file is in your QMAKEFEATURES 
directory.


Dependencies:
 ALTA core: Eigen
 Plugin rational_eigen:     Eigen 
 Plugin rational_quadprog:  Quadprog++ 
 Plugin rational_cgal:      The CGAL library
 Plugin rational_parallel:  The OpenMP library, Quadprog++ library and Eigen
 Plugin rational_matlab:    Matlab engine (matlab.prf required)
 Plugin nonlinear_eigen:    Eigen
 Plugin nonlinear_ceres:    CERES library and its dependancies
 Plugin nonlinear_nlopt:    NLOpt library and its dependancies
 Plugin nonlinear_ipopt:    IpOpt library and its dependancies

Eigen Plugins
 You must provide an eigen.prf file that contains

 INCLUDEPATH *= PATH_TO_EIGEN_DIRECTORY

Quadprog++ Plugins

 We provide our own version of quadprog++ which uses Eigen library.
 To compile it:
  
 Go to external/quadprog++/
 Use qmake (to generale the Makefile)
 make 
 
 Then create a quadprog.prf file that will include 
 LIBS *= PATH_TO_LIBQUADPROG/libquadprog++.a
 QMAKE_LIBDIR *= PATH_TO_LIBQUADPROG
 INCLUDEPATH *= PATH_TO_LIBQUADPROG_HEADERS
 

Parallel Plugin that requires OpenMP

 Create an  openmp.prf file and add the following directives:
 QMAKE_CXXFLAGS *=-fopenmp
 



Matlab Plugin

 Create an matlab.prf file and add the following directives:
 The PATH_TO_MATLAB_INCLUDE_DIRECTORY must point to a directory that
 contains the file engine.h

 INCLUDEPATH *= PATH_TO_MATLAB_INCLUDE_DIRECTORY
 




Generate the documentation using Doxygen
  cd ${ALTA}/documents/
  doxygen doxygen.conf

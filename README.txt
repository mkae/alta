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
configs/:   contains platform dependent configuration files for scons.
documents/: contains the documentation in doxygen format.
external/:  contains any third party library that needs to be used by ALTA.
            Contains a modified QuadProg++ library using Eigen. You can
            download and compile automaticaly some third party libraries using
            provided python scripts obtain_[libname].py.
sources/:   contains all the source files. 


2. Building
===========

ALTA uses SCons, a Python-based build system:

  http://scons.org/

To build ALTA, run a command like the following from the top-level source
directory (noted ${ALTA} in the remainder of this document):

  $ scons

If you want to use a specifi configuration (i.e. change the default compiler,
specify where some libaries should be found) you can use a configuration file
using the '--cfg [file]' option:

  $ scons --cfg=[config-file.py]

Here, 'config-file.py' must be replaced with a suitable configuration file for
your platform. For instance, when building with GCC on GNU/Linux, you may run:

  $ scons --cfg=configs/scons/config-linux-gcc.py

ALTA provide a Python interface for non command-line experts. However, this
interface is not build automatically. To build the Python interface, please
run the following command at the root of the repository:

  $ scons python

You can as well use a plateform specific configuration file using the '--cfg'
option.


3. Dependencies:
================

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


4. Use ALTA
===========

Once ALTA is compiled with its plugins, you can access the executables and
library in ${ALTA}/sources/build. To access them directly from the shell, you
can source the `setpath.sh` script at the root:

  $ source setpath.sh

This will expose the binary and plugins to the systeme and allow you to run
ALTA commands from anywhere.

For further use of ALTA, please refer to the documentation and tutorials.


5. Generate the documentation using Doxygen
===========================================

If you have doxygen installed on your system, you can build the documentation
(a static website) using the scons script:

  $ scons doc

The static website is then available at:

  ${ALTA}/documents/doxygen/html/index.html

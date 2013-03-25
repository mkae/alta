This is the repository for the ALTA Project

Repository Organization

data/:      contains the data for which the fitting techniques are tested against.
documents/: contains the documentation, should be build using doxygen into that
            directiory.
sources/:   Contains all the source files. A Makefile or VS project can be created
            there from the .pro file.



Building Advises

We use heavily the Qt profile functionality. To build some of the plugins you will
be required to create your own system dependant .prf for any used library. For example
all the rational BRDF fitters use the Eigen library. Therefore it is mandatory that you
provide a eigen.prf file and that this file is in your QMAKEFEATURES directory.

Generate the documentation using Doxygen
  cd ${ALTA}/documents/
  doxygen doxygen.conf

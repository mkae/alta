import obtain
import os
import sys
import shutil


# Download GLOG
obtain.obtain('GLOG v0.3.3', 'glog-0.3.3', 'http://google-glog.googlecode.com/files/glog-0.3.3.tar.gz', 'glog-0.3.3.tar.gz')
if os.name == 'nt':
	print '<<WARNING>> no automatic installation for this package'
else:
	print '<<INSTALL>> configure and build GLOG v0.3.3'
	obtain.configure_build('glog-0.3.3')
#end

# Download Eigen3
execfile('obtain_eigen.py')
 
# Download CERES
obtain.obtain('CERES v1.7.0', 'ceres-solver-1.7.0', 'http://ceres-solver.googlecode.com/files/ceres-solver-1.7.0.tar.gz', 'ceres-solver-1.7.0.tar.gz')
print '<<WARNING>> CERES installation requires CMake. You need to run it yourself'

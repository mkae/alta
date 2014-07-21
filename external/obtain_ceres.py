import obtain
import os
import sys
import shutil
import subprocess

# Download GLOG
obtain.obtain('GLOG v0.3.3', 'glog-0.3.3', 'http://google-glog.googlecode.com/files/glog-0.3.3.tar.gz', 'glog-0.3.3.tar.gz')

if not os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'glog'):
	if os.name == 'nt':
		print '<<WARNING>> no automatic installation for this package'
	else:
            if sys.platform == 'darwin':
                obtain.patch('glog-0.3.3/src/glog/stl_logging.h.in', 'glog.patch')
            #end
	    print '<<INSTALL>> configure and build GLOG v0.3.3'
            obtain.configure_build('glog-0.3.3', '--enable-static=no --enable-shared=true --with-pic=true')
	#end
else:
	print '<<INSTALL>> GLOG already installed'
#end


# Download Eigen3
execfile('obtain_eigen.py')
 
# Download CERES
version   = '1.9.0'
base      = 'ceres-solver'
name      = 'CERES v' + version
directory = base + '-' + version
url       = 'http://ceres-solver.org/ceres-solver-' + version + '.tar.gz'
filename  = 'ceres-solver-' + version + '.tar.gz'
obtain.obtain(name, directory, url, filename)

## Test for the presence of already compiled ceres version in
## the $ALTA/external/build directory. Then test for the
## presence of cmake.
compile_test = not os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'ceres')

with open(os.devnull, 'w') as fnull:
	res = subprocess.call(['cmake', '--version'], stdout = fnull, stderr = fnull, shell=True)
	if res != 0:
		compile_test = False
		print '<<ERROR>> cmake is not installed'
	#end
#end


if  compile_test:
	print '<<INSTALL>> configure and build CERES'
	os.chdir('.' + os.sep + 'ceres-solver-' + version)
	build_dir = os.pardir + os.sep + 'build' + os.sep

	libname = ''
	if os.name == 'posix':
		libname = 'libglog.so'
	elif os.name == 'nt':
		libname = 'glog.lib'
	else:
		libname = 'libglog.dylib'
	#end

	#cmake_cmd = 'cmake -DGLOG_LIB=' + build_dir + 'lib' + os.sep + libname + ' -DGLOG_INCLUDE=' + build_dir + 'include -DGFLAGS=OFF ' + '-DEIGEN_INCLUDE=' + build_dir + 'include -DCMAKE_INSTALL_PREFIX=' + build_dir + ' .' + ' -DDISABLE_TR1=ON -DBUILD_EXAMPLES=OFF ' + '-DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DMINIGLOG=OFF'
	cmake_cmd = 'cmake -DBUILD_SHARED_LIBS=ON -DGLOG_LIB=' + build_dir + 'lib' + ' -DGLOG_INCLUDE=' + build_dir + 'include -DGFLAGS=OFF ' + '-DEIGEN_INCLUDE=' + build_dir + 'include -DCMAKE_INSTALL_PREFIX=' + build_dir + ' .' + ' -DDISABLE_TR1=ON -DBUILD_EXAMPLES=OFF ' + '-DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DMINIGLOG=ON'
	
	if os.name == 'nt':
		ret = os.system(cmake_cmd + ' -G \"NMake Makefiles\"')
		ret = os.system('nmake install')
	else:
		cmake_cmd += ' -DCMAKE_CXX_FLAGS=\"-fPIC\"'
		cmake_cmd += ' -DCMAKE_C_FLAGS=\"-fPIC\"'
		ret = os.system(cmake_cmd)
		ret = os.system('make install')
	#end
	
	os.chdir(os.pardir)
else:
	print '<<INSTALL>> CERES already installed or cannot be installed automatically'
#end

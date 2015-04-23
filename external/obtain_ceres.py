import obtain
import os
import sys
import shutil
import subprocess
import SCons.Warnings as W
import SCons.SConf as C

# Download GLOG
obtain.obtain('GLOG v0.3.3', 'glog-0.3.3',
              'http://google-glog.googlecode.com/files/glog-0.3.3.tar.gz', 'glog-0.3.3.tar.gz',
              'fbf90c2285ba0561db7a40f8a4eefb9aa963e7d399bd450363e959929fe849d0')

if not os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'glog'):
	if os.name == 'nt':
		W.warn(obtain.AltaDependencyWarning, 'sorry, no automatic installation of GLOG')
	else:
		exists_archive = os.path.exists('.' + os.sep + 'glog-0.3.3')
		if sys.platform == 'darwin' and not exists_archive:
			obtain.patch('glog-0.3.3/src/glog/stl_logging.h.in', 'glog.patch')
		C.progress_display('configuring and building GLOG')
		obtain.configure_build('glog-0.3.3', '--enable-static=no --enable-shared=true --with-pic=true')
else:
	C.progress_display('GLOG is already installed')


# Download CERES.  Assume Eigen is already available.
version   = '1.9.0'
base      = 'ceres-solver'
name      = 'CERES v' + version
directory = base + '-' + version
url       = 'http://ceres-solver.org/ceres-solver-' + version + '.tar.gz'
filename  = 'ceres-solver-' + version + '.tar.gz'
sha256    = '30ac0729249f908afe80cb6fd06ae6d037f25a60d9fac54f61344389adab9c1a'
obtain.obtain(name, directory, url, filename, sha256)

## Test for the presence of already compiled ceres version in
## the $ALTA/external/build directory. Then test for the
## presence of cmake.
compile_test = not os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'ceres')

with open(os.devnull, 'w') as fnull:
	res = subprocess.call(['cmake', '--version'], stdout = fnull, stderr = fnull, shell=True)
	if res != 0:
		compile_test = False
		W.warn(obtain.AltaDependencyWarning,
					 'CMake is not installed but is needed to build CERES')


if compile_test:
	C.progress_display('configuring and building CERES')
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
	cmake_cmd = 'cmake -DBUILD_SHARED_LIBS=OFF -DGLOG_LIB=' + build_dir + 'lib' + ' -DGLOG_INCLUDE=' + build_dir + 'include -DGFLAGS=OFF ' + '-DEIGEN_INCLUDE=' + build_dir + 'include -DCMAKE_INSTALL_PREFIX=' + build_dir + ' .' + ' -DDISABLE_TR1=ON -DBUILD_EXAMPLES=OFF ' + '-DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DMINIGLOG=ON'
	
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
	W.warn(obtain.AltaDependencyWarning,
				 'CERES already installed or cannot be installed automatically')

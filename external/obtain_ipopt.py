import obtain
import os
import sys
import shutil
import subprocess
import SCons.Warnings as W
import SCons.SConf as C

# Download IpOpt.
obtain.obtain('IpOpt v3.11.8', 'Ipopt-3.11.8',
              'http://www.coin-or.org/download/source/Ipopt/Ipopt-3.11.8.tgz', 'Ipopt-3.11.8.tgz',
              '9f9b76075fbd9315286ea4d7c159c94cab4a4fb16122fb172b24910af5b5b75b')

## Test for the presence of already compiled ceres version in
## the $ALTA/external/build directory. Then test for the
## presence of cmake.
compile_test = not os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'Ipopt')


if os.path.exists('.' + os.sep + 'Ipopt-3.11.8.tgz') :
	C.progress_display('the IpOpt package is already downloaded or installed')
	C.progress_display('if the plugins using IpOpt fail to build, \
check this installation')

else:

	if not os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'coin/IpIpoptNLP.hpp'):
		if os.name == 'nt':
			W.warn(obtain.AltaDependencyWarning, 'sorry, no automatic installation of IpOpt')
		else:
			C.progress_display('configuring and building Ipopt and its dependencies')

			path  = os.getcwd()
			third = path + os.sep + 'Ipopt-3.11.8' + os.sep + 'ThirdParty' + os.sep
			os.chdir(third + 'Blas')
			ret = os.system('./get.Blas')
			os.chdir(third + 'Lapack')
			ret = os.system('./get.Lapack')
			os.chdir(third + 'ASL')
			ret = os.system('./get.ASL')
			os.chdir(third + 'Mumps')
			ret = os.system('./get.Mumps')

			os.chdir(path)
			obtain.configure_build('Ipopt-3.11.8', '--enable-dependency-linking')
		#end
	else:
		C.progress_display('IpOpt is already installed')
	#end
#end

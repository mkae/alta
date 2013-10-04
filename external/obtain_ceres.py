import urllib
import os
import sys
import shutil

# Download a file
def download(url, filename):
	urllib.urlretrieve(url, filename)
#end

# Uncompress the archive
def uncompress(filename):
	cmd = '';
	if os.name == 'posix':
		cmd = 'tar -xf '
	elif os.name == 'nt':
		cmd = '7zip '

	cmd += filename
	ret = os.system(cmd)
#end

# Obtain a package and uncompress it
def obtain(name, rep, url, filename):
	if not os.path.exists(rep):
		print '<<INSTALL>> Obtaining ' + name

		download(url, filename)
		uncompress(filename)
	else:
		print '<<INSTALL>> ' + name + ' already available'
	#end
#end

def configure_build(rep):
	os.chdir(rep)
	ret = os.system('./configure -q --prefix=' + os.getcwd() + os.sep + os.pardir + os.sep + 'build')
	if ret != 0:
		print '<<ERROR>> unable to configure package'
		exit
	#end
	ret = os.system('make -s && make -s install')
	if ret != 0:
		print '<<ERROR>> unable to build & install package'
		exit
	#end
	os.chdir(os.pardir)
#end

# Check if the build dir exists
if not os.path.exists('build'):
	os.mkdir('build')
#end

# Download GLOG
obtain('GLOG v0.3.3', 'glog-0.3.3', 'http://google-glog.googlecode.com/files/glog-0.3.3.tar.gz', 'glog-0.3.3.tar.gz')
if os.name == 'nt':
	print '<<WARNING>> no automatic installation for this package'
else:
	print '<<INSTALL>> configure and build GLOG v0.3.3'
	configure_build('glog-0.3.3')
#end

# Download Eigen3
obtain('Eigen v3.2.1', 'eigen-eigen-ffa86ffb5570', 'http://bitbucket.org/eigen/eigen/get/3.2.0.tar.gz', 'eigen-3.2.1.tar.gz')
shutil.copytree('eigen-eigen-ffa86ffb5570' + os.sep + 'Eigen', 'build' + os.sep + 'include' + os.sep + 'Eigen')
 
# Download CERES
obtain('CERES v1.7.0', 'ceres-solver-1.7.0', 'http://ceres-solver.googlecode.com/files/ceres-solver-1.7.0.tar.gz', 'ceres-solver-1.7.0.tar.gz')
print '<<WARNING>> CERES installation requires CMake. You need to run it yourself'

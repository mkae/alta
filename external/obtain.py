import urllib
import os
import sys
import shutil
import tarfile

# Check if the build dir exists
if not os.path.exists('build'):
	os.mkdir('build')
#end

# Download a file
def download(url, filename):
	urllib.urlretrieve(url, filename)
#end

# Uncompress the archive
def uncompress(filename):
	tfile = tarfile.open(filename, 'r:gz')
	tfile.extractall()

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


import urllib
import os
import sys
import shutil
import tarfile
import hashlib
import SCons.Errors

# Check if the build dir exists
if not os.path.exists('build'):
	os.mkdir('build')
#end

# Make sure that FILENAME has the given SHA256, and raise an error
# otherwise.
def check_integrity(filename, sha256):
        hasher = hashlib.sha256()
        BLOCK_SIZE = 65536
        with open(filename, 'rb') as input:
                buf = input.read(BLOCK_SIZE);
                while len(buf) > 0:
                        hasher.update(buf)
                        buf = input.read(BLOCK_SIZE)
        obtained = hasher.hexdigest()
        if obtained != sha256:
                # Keep the faulty file around but move it out of the
                # way so we don't end up using it by mistake.
                os.rename(filename, filename + ".wrong-hash")

                print >> sys.stderr, "error: downloaded file '" + \
                        filename + "' is inauthentic"
                print >> sys.stderr, \
                        "error: got sha256 hash {0} but expected {1}".format(obtained, sha256)
                raise SCons.Errors.BuildError("downloaded file is inauthentic",
                                              1, 1, filename)

# Download from URL to FILENAME.
def download(url, filename):
	urllib.urlretrieve(url, filename)

# Uncompress the archive
def uncompress(filename):
	tfile = tarfile.open(filename, 'r:gz')
	tfile.extractall()
#end

# Apply a patch to some file
def patch(filename, patch):
	ret = os.system('patch ' + filename + ' ' + patch)
#end

# Obtain a package, check its integrity, and uncompress it.
def obtain(name, rep, url, filename, sha256):
	if not os.path.exists(rep):
		print '<<INSTALL>> Obtaining ' + name

                if not os.path.exists(filename):
                        download(url, filename)
                check_integrity(filename, sha256)
		uncompress(filename)
	else:
		print '<<INSTALL>> ' + name + ' already available'
	#end
#end

def configure_build(rep, options = ''):
	os.chdir(rep)
	ret = os.system('./configure -q --prefix=' + os.getcwd() + os.sep + os.pardir + os.sep + 'build ' + options)
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


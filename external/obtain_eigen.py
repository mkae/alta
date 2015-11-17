import obtain
import os
import sys
import shutil
import SCons.Errors as Errors

# Download Eigen3
version   = '3.2.7'
base      = 'Eigen'
name      = 'Eigen v' + version
directory = 'eigen-eigen-b30b87236a1b'
url       = 'http://bitbucket.org/eigen/eigen/get/' + version + '.tar.bz2'
filename  = 'eigen-' + version + '.tar.bz2'
sha256    = 'e58e1a11b23cf2754e32b3c5990f318a8461a3613c7acbf6035870daa45c2f3e'
obtained  = obtain.obtain(name, directory, url, filename, sha256)

if obtained:
   rep = 'build' + os.sep + 'include' + os.sep + 'Eigen'
   if not os.path.exists(rep):
      shutil.copytree(directory + os.sep + 'Eigen', rep)

   unsup_rep = 'build' + os.sep + 'include' + os.sep + 'unsupported'
   if not os.path.exists(unsup_rep):
      shutil.copytree(directory + os.sep + 'unsupported', unsup_rep)

else:
   msg = "download from '" + url + "' failed"
   raise Errors.BuildError(errstr = msg,
                           filename = filename,
                           exitstatus = 1)

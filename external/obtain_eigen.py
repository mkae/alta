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
url       = 'http://bitbucket.org/eigen/eigen/get/' + version + '.tar.gz'
filename  = 'eigen-' + version + '.tar.gz'
sha256    = '5a50a006f83480a31f1f9beabec9e91dad95138df19363ee73ccf57676f10405'
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

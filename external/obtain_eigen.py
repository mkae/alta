import obtain
import os
import sys
import shutil
import SCons.Errors as Errors

# Download Eigen3
version   = '3.2.1'
base      = 'Eigen'
name      = 'Eigen v' + version
directory = 'eigen-eigen-6b38706d90a9'
url       = 'http://bitbucket.org/eigen/eigen/get/' + version + '.tar.gz'
filename  = 'eigen-' + version + '.tar.gz'
sha256    = 'fa9b1821608d8fd3b364ac9db62102b797364923fa0c74e6cbc4f9ba36c43e44'
obtained  = obtain.obtain(name, directory, url, filename, sha256)

if obtained:
   rep = 'build' + os.sep + 'include' + os.sep + 'Eigen'
   if not os.path.exists(rep):
      shutil.copytree(directory + os.sep + 'Eigen', rep)

   unsup_rep = 'build' + os.sep + 'include' + os.sep + 'unsupported'
   if not os.path.exists(unsup_rep):
      shutil.copytree(directory + os.sep + 'unsupported', unsup_rep)

else:
   msg = "Unable to obtain " + name
   raise Errors.BuildError(msg, 1, 1, filename)
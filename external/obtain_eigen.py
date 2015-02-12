import obtain
import os
import sys
import shutil

# Download Eigen3
directory = 'eigen-eigen-6b38706d90a9'
obtain.obtain('Eigen v3.2.1', directory,
              'http://bitbucket.org/eigen/eigen/get/3.2.1.tar.gz', 'eigen-3.2.1.tar.gz',
              'fa9b1821608d8fd3b364ac9db62102b797364923fa0c74e6cbc4f9ba36c43e44')

rep = 'build' + os.sep + 'include' + os.sep + 'Eigen'
if not os.path.exists(rep):
	shutil.copytree(directory + os.sep + 'Eigen', rep)
#end

unsup_rep =  rep + os.sep+ 'unsupported'
if not os.path.exists(unsup_rep):
	shutil.copytree(directory + os.sep + 'unsupported', unsup_rep)
#end

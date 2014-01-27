import obtain
import os
import sys
import shutil

# Download Eigen3
obtain.obtain('Eigen v3.2.1', 'eigen-eigen-ffa86ffb5570', 'http://bitbucket.org/eigen/eigen/get/3.2.0.tar.gz', 'eigen-3.2.1.tar.gz')

rep = 'build' + os.sep + 'include' + os.sep + 'Eigen'
if not os.path.exists(rep):
	shutil.copytree('eigen-eigen-ffa86ffb5570' + os.sep + 'Eigen', rep)
#end

if not os.path.exists(rep + os.sep + 'unsupported'):
	shutil.copytree('eigen-eigen-ffa86ffb5570' + os.sep + 'unsupported', rep + os.sep + 'Eigen' + os.sep+ 'unsupported')
#end

import obtain
import os
import sys
import shutil

# Download Eigen3
obtain.obtain('Eigen v3.2.1', 'eigen-eigen-ffa86ffb5570',
              'http://bitbucket.org/eigen/eigen/get/3.2.0.tar.gz', 'eigen-3.2.1.tar.gz',
              'ea084fa2517d33f3a3a99c1d4be9e114e91ebfd174cacbd897d90e7b9d6d91c3')

rep = 'build' + os.sep + 'include' + os.sep + 'Eigen'
if not os.path.exists(rep):
	shutil.copytree('eigen-eigen-ffa86ffb5570' + os.sep + 'Eigen', rep)
#end

unsup_rep =  rep + os.sep+ 'unsupported'
if not os.path.exists(unsup_rep):
	shutil.copytree('eigen-eigen-ffa86ffb5570' + os.sep + 'unsupported', unsup_rep)
#end

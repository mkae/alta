import obtain
import os
import sys
import shutil


# Download NlOpt
if not os.path.exists('.' + os.sep + 'nlopt-2.4.1.tar.gz'):
   obtain.obtain('NlOpt', 'nlopt-2.4.1',
              'http://ab-initio.mit.edu/nlopt/nlopt-2.4.1.tar.gz', 'nlopt-2.4.1.tar.gz',
              'fe9ade54ed79c87f682540c34ad4e610ff32c9a43c52c6ea78cef6adcd5c1319')

if not os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'nlopt.hpp'):
	if os.name == 'nt':
		print '<<WARNING>> no automatic installation for this package'
	else:
		print '<<INSTALL>> configure and build Nlopt v2.4.1'
		#obtain.configure_build('nlopt-2.4.1', '--enable-static --without-matlab --without-octave --without-python --without-guile')
		obtain.configure_build('nlopt-2.4.1', '--enable-shared --without-matlab --without-octave --without-python --without-guile  --with-pic')
	#end
else:
	print '<<INSTALL>> NlOpt already installed'
#end

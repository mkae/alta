import obtain
import os
import sys
import shutil
import subprocess
import SCons.Warnings as W
import SCons.SConf as C
from subprocess import Popen, PIPE

def getThirdParty(name) :
   path  = os.getcwd()
   third = path + os.sep + directory + os.sep + 'ThirdParty' + os.sep
   os.chdir(third + name)
   ret   = Popen(['./get.' + name]).wait()
   if(ret != 0):
      C.progress_display('Failed to get IpOpt third party \'' + name + '\'')
   os.chdir(path)

# Download IpOpt.
version   = '3.12.4'
base      = 'Ipopt'
name      = 'IpOpt v' + version
directory = base + '-' + version
url       = 'http://www.coin-or.org/download/source/Ipopt/Ipopt-' + version + '.tgz'
filename  = 'Ipopt-' + version + '.tgz'
sha256    = '292afd952c25ec9fe6225041683dcbd3cb76e15a128764671927dbaf881c2e89'
obtained  = obtain.obtain(name, directory, url, filename, sha256)

compiled  = os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'coin/IpIpoptNLP.hpp')

if obtained and not compiled:
   if os.name == 'nt':
      W.warn(obtain.AltaDependencyWarning, 'sorry, no automatic installation of IpOpt')
   else:
      C.progress_display('configuring and building ' + name + ' and its dependencies')

      getThirdParty('Blas')
      getThirdParty('Lapack')
      getThirdParty('ASL')
      getThirdParty('Mumps')

      obtain.configure_build(directory,
                             ['--enable-static', '--with-pic',
                              '--enable-dependency-linking'])
else:
   C.progress_display('IpOpt is already installed')
   C.progress_display('if the plugins using IpOpt fail to build, check its installation')

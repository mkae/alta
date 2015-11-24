import obtain
import os
import sys
import shutil
import SCons.SConf as C

# Download NlOpt.
version   = '2.4.1'
base      = 'nlopt'
name      = 'NlOpt v' + version
directory = base + '-' + version
baseurl   = 'http://ab-initio.mit.edu/nlopt/'
url       = baseurl + 'nlopt-' + version + '.tar.gz'
filename  = 'nlopt-' + version + '.tar.gz'
sha256    = 'fe9ade54ed79c87f682540c34ad4e610ff32c9a43c52c6ea78cef6adcd5c1319'
obtained  = obtain.obtain(name, directory, url, filename, sha256)

compiled  = os.path.exists('.' + os.sep + 'build' + os.sep + 'include' + os.sep + 'nlopt.hpp')
 
if obtained and not compiled:
   C.progress_display('configuring and building ' + name)
   if os.name == 'nt':

      # Get the CMakeList file
      obtain.download(baseurl + 'CMakeLists.txt', directory+os.sep+'CMakeLists.txt')

      # Get the Additional CMake script
      obtain.download(baseurl + 'config.cmake.h.in', directory+os.sep+'config.cmake.h.in')

      # Build cmake script
      obtain.cmake_build(directory)

   else:
      obtain.configure_build(directory,
                             ['--enable-static', '--with-pic',
                              '--without-matlab',
                              '--without-octave', '--without-python',
                              '--without-guile'])
else:
   print '<<INSTALL>> NlOpt already installed'

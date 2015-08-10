# ALTA --- Analysis of Bidirectional Reflectance Distribution Functions
#
# Copyright (C) 2014, 2015 CNRS
# Copyright (C) 2013, 2014, 2015 Inria
# Copyright (C) 2015 Universite de Montreal
#
# This file is part of ALTA.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0.  If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
import sys
import SCons.SConf as C

## Build the documentation. This is independant of everything else and
## should return afterwards, speeding-up the process of generating doc.
##
if 'doc' in COMMAND_LINE_TARGETS:
   os.chdir('documents')
   Alias('doc', os.system('doxygen doxygen.conf'))
   Exit(0)

## Add ALTA custom cmd configurations
##
AddOption('--cfg', help='Specify a configuration file')
AddOption('--no-externals', action="store_false", dest="obtain_externals", default=True, help='Do not download and build externals')


## Import configuration from a config file
##
configFile = GetOption('cfg')
if configFile == None:
	
	if sys.platform == 'win32':
		configFile = "./configs/scons/config-windows-cl.py"
	elif sys.platform == 'darwin':
		configFile = "./configs/scons/config-macos-clang.py"
	elif sys.platform == 'linux2':
		configFile = "./configs/scons/config-linux-gcc.py"
	else:
		print '<<ERROR>> You need to specify a configuration file using:'
		print '<<ERROR>>    scons --cfg=[filename]'
		print '<<ERROR>> Please find example of configuration files in ${ALTA}/configs/scons/'
		Exit(1)
#end

if not os.path.exists(configFile):
	print '<<ERROR>> the config file you specified \"' + configFile + '\" does not exists'
	Exit(1)
else:
	print '<<INFO>> Using config file \"' + configFile + '\"'
#end

# List of C++ compilers we look for by default.
cxx_compilers = [ 'g++', 'c++', 'clang++', 'cl' ]

def default_cxx_compiler():
	"Return the absolute file name of a C++ compiler."
	for cxx in cxx_compilers:
		cxx_full = WhereIs(cxx)
		if cxx_full != None:
			return cxx_full
	return cxx_compilers[0]

vars = Variables(configFile)
vars.Add('CXX',               'C++ compiler',
				 default = default_cxx_compiler())
vars.Add('CCFLAGS',           'Compiler\'s flags',
		 default = ['-g', '-O2', '-Wall'])
vars.Add('LINKFLAGS',         'Linker\'s flags',
		 default = [])
vars.Add('PLUGIN_LIB',        'Special links for ALTA plugin')
vars.Add('EIGEN_INC',         'Eigen include directory (mandatory)')
vars.Add('PYTHON_INC',        'Python and boost-python include directory')
vars.Add('PYTHON_DIR',        'Python and boost-python libraries directory')
vars.Add('PYTHON_LIB',        'Python and boost-python libraries', default = [])
vars.Add('OPENEXR_INC',       'OpenEXR include directory', default=[])
vars.Add('OPENEXR_LIB',       'OpenEXR libraries', default = [])
vars.Add('OPENEXR_DIR',       'OpenEXR libraries directory', default=[])
vars.Add('FLANN_INC',         'FLANN include directory')
vars.Add('FLANN_DIR',         'FLANN libraries directory')
vars.Add('FLANN_LIB',         'FLANN libraries')
vars.Add('OPENMP_FLAG',       'OpenMP required flags')
vars.Add('OPENMP_LIB',        'OpenMP libraries')
vars.Add('QUADPROG_INC',      'QUADPROG include directory')
vars.Add('QUADPROG_DIR',      'QUADPROG libraries directory')
vars.Add('QUADPROG_LIB',      'QUADPROG libraries')
vars.Add('CERES_INC',         'CERES include directory')
vars.Add('CERES_DIR',         'CERES libraries directory')
vars.Add('CERES_LIB',         'CERES libraries')
vars.Add('CERES_OPT_LIB',     'CERES optional libraries')
vars.Add('NLOPT_INC',         'NLOPT include directory')
vars.Add('NLOPT_DIR',         'NLOPT libraries directory')
vars.Add('NLOPT_LIB',         'NLOPT libraries', default = [])
vars.Add('NLOPT_OPT_LIB',     'NLOPT optional libraries', default = [])
vars.Add('IPOPT_INC',         'IPOPT include directory')
vars.Add('IPOPT_DIR',         'IPOPT libraries directory')
vars.Add('IPOPT_LIB',         'IPOPT libraries', default = [])
vars.Add('IPOPT_OPT_LIB',     'IPOPT optional libraries')
vars.Add('MATLAB_INC',        'MATLAB include directory')
vars.Add('MATLAB_DIR',        'MATLAB directory')
vars.Add('MATLAB_LIB',        'MATLAB libraries')


##
# Copy the system environment.
#
# Update the PKG_CONFIG_PATH variable to add the package configuration
# files contained is #external/build/lib/pkgconfig
##
envVars = {}
for var in [ 'PATH', 'CPATH', 'LIBRARY_PATH', 'PKG_CONFIG_PATH', 'TMP', 'TMPDIR' ]:
	if var in os.environ:
		envVars[var] = os.environ[var]
	else:
		envVars[var] = '';

if len(envVars['PKG_CONFIG_PATH']) > 0:
	envVars['PKG_CONFIG_PATH'] += ':'
envVars['PKG_CONFIG_PATH'] += os.path.abspath('external' + os.sep + 'build' + os.sep + 'lib' + os.sep + 'pkgconfig')
env = Environment(variables = vars, ENV = envVars )
env['DL_EXTERNALS'] = GetOption('obtain_externals')


# Generate help text for the build variables.
Help(vars.GenerateHelpText(env))

C.progress_display('the current platform is: ' + env['PLATFORM'])

## PLATFORM dependant section
##
if env['PLATFORM'] == 'darwin':

	# Adding the /usr/local/lib directory used to store libraries of
	# MacPorts or Brew.
	env.AppendUnique(LIBPATH = ['/usr/local/lib'])
	env.AppendUnique(CPPPATH = ['/usr/local/include'])

	env.AppendUnique(LIBPATH = ['/opt/local/lib/'])
	env.AppendUnique(CPPPATH = ['/opt/local/include/'])


def CheckPKG(context, name):
		"""Return True if package NAME can be found with 'pkg-config'."""
		context.Message('Checking for %s using pkg-config... ' % name)
		ret = context.TryAction('pkg-config --exists \'%s\'' % name)[0]
		context.Result(ret)
		return ret


"""
Return True if the given library is available. First look for the LIB_VAR and
INC_VAR construction variables, honoring them if they are defined. Then look
for PKGSPEC using pkg-config. Last, try to build LIB with HEADER. Configure
ENV accordingly.
"""
def library_available(env, pkgspec='', lib='', header='',
					  language='c++', inc_var='', lib_var=''):

	conf = Configure(env, custom_tests = { 'CheckPKG' : CheckPKG })

	result = False
	# If a XXX_LIB is specified in the environment, add the various path
	# and link flags. Check if the library is correctly compiling and
	# linking with the header.
	if (lib in env) and (len(env[lib]) > 0):
		env.AppendUnique(LIBPATH = env[lib_var])
		env.AppendUnique(CPPPATH = env[inc_var])
		env.AppendUnique(LIBS = env[lib])

		# Check whether the library is usable.
		result = conf.CheckLibWithHeader(env[lib], header, language)

	elif conf.CheckPKG(pkgspec):
		env.ParseConfig('pkg-config --cflags --libs "' + pkgspec + '"')
		result = True

	conf.Finish()
	return result

def openexr_available(env):
		"""Return True if OpenEXR is available."""
		return library_available(env, pkgspec='OpenEXR',
								 inc_var='OPENEXR_INC',
								 lib_var='OPENEXR_DIR',
								 lib='OPENEXR_LIB',
								 header='ImfRgbaFile.h')

# Export these for use in SConscripts.
Export('CheckPKG', 'library_available', 'openexr_available')

conf = Configure(env)

# Determine the extra libraries that libcore (and thus everything
# else) depends on.  Plugins need to specify it in addition to -lcore
# because libcore is not a shared library.
ALTA_LIBS = []

# Libcore's uses 'dlopen', which is in libdl in GNU libc.
if conf.CheckLibWithHeader('dl', 'dlfcn.h', 'c++'):
		ALTA_LIBS = ['dl']

# Libcore uses 'clock_gettime', which is in librt in GNU libc.
if conf.CheckLibWithHeader('rt', 'sched.h', 'c++'):
		ALTA_LIBS = ALTA_LIBS + ['rt']

Export('ALTA_LIBS')

conf.Finish()

## Load the configuration file if it exists. The configuration file
## is a python script that updates the env variable with different
## paths.
##
env.AppendUnique(LIBPATH = ['#external/build/lib'])
env.AppendUnique(LIBPATH = ['#sources/build'])
env.AppendUnique(CPPPATH = ['#external/build/include'])
env.AppendUnique(CPPPATH = ['#sources'])

# Consider files changed as soon as their modification time changes.
env.Decider('timestamp-newer')

## Launch the compilations
##l
Export('env')

external = env.SConscript('external/SConscript')
core     = env.SConscript('sources/core/SConscript')
plugins  = env.SConscript('sources/plugins/SConscript')
softs    = env.SConscript('sources/softs/SConscript')

if 'python' in COMMAND_LINE_TARGETS:
	python = env.SConscript('sources/python/SConscript')
	env.Depends(python, core)

if 'tests' in COMMAND_LINE_TARGETS:
	tests = env.SConscript('sources/tests/SConscript')
	env.Depends(tests, core)
	if 'python' in COMMAND_LINE_TARGETS:
		env.Depends(tests, python)

env.Depends(plugins, core)
env.Depends(softs, core)

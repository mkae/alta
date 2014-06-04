import os
import sys

## Add ALTA custom cmd configurations
##
AddOption('--cfg', help='Specify a configuration file (see config.example')


## Import configuration from a config file
##
configFile = GetOption('cfg')
if configFile == None:
	print '<<ERROR>> You need to specify a configuration file using:'
	print '<<ERROR>>    scons --cfg [filename]'
	print '<<ERROR>> Please find example of configuration files in ${ALTA}/configs/scons/'
	Exit(1)
#end

if not os.path.exists(configFile):
	print '<<ERROR>> the config file you specified does not exists'
	Exit(1)
#end

vars = Variables(configFile)
vars.Add('CXX',               'Compiler')
vars.Add('CCFLAGS',           'Compiler\'s flags')
vars.Add('LINKFLAGS',         'Linker\'s flags')
vars.Add('OPENEXR_INC',       'OpenEXR include directory')
vars.Add('OPENEXR_DIR',       'OpenEXR libraries directory')
vars.Add('OPENEXR_LIBS',      'OpenEXR libraries')
vars.Add('OPENMP_FLAGS',      'OpenMP required flags')
vars.Add('OPENMP_LIBS',       'OpenMP libraries')
vars.Add('QUADPROG_INC',      'QUADPROG include directory')
vars.Add('QUADPROG_DIR',      'QUADPROG libraries directory')
vars.Add('QUADPROG_LIBS',     'QUADPROG libraries')
vars.Add('CERES_INC',         'CERES include directory')
vars.Add('CERES_DIR',         'CERES libraries directory')
vars.Add('CERES_LIBS',        'CERES libraries')
vars.Add('CERES_OPT_LIBS',    'CERES optional libraries')
vars.Add('NLOPT_INC',         'NLOPT include directory')
vars.Add('NLOPT_DIR',         'NLOPT libraries directory')
vars.Add('NLOPT_LIBS',        'NLOPT libraries')
vars.Add('NLOPT_OPT_LIBS',    'NLOPT optional libraries')
vars.Add('IPOPT_INC',         'IPOPT include directory')
vars.Add('IPOPT_DIR',         'IPOPT libraries directory')
vars.Add('IPOPT_LIBS',        'IPOPT libraries')
vars.Add('IPOPT_OPT_LIBS',    'IPOPT optional libraries')
vars.Add('MATLAB_INC',        'MATLAB include directory')
vars.Add('MATLAB_DIR',        'MATLAB directory')
vars.Add('MATLAB_LIBS',       'MATLAB libraries')


env = Environment(variables = vars)


## PLATFORM dependant section
##
if sys.platform == 'darwin':

	# Adding the /usr/local/lib directory used to store libraries of
	# MacPorts or Brew.
	env.AppendUnique(LIBPATH = ['/usr/local/lib'])
	env.AppendUnique(CPPPATH = ['/usr/local/include'])

	env.AppendUnique(LIBPATH = ['/opt/local/lib/'])
	env.AppendUnique(CPPPATH = ['/opt/local/include/'])

#end


## COMPILER dependant section
##
if env['CC'] in ['gcc', 'clang'] or env['CXX'] in ['g++', 'clang']:
	print '<<INFO>> add -fPIC option to the compiler'
	env.AppendUnique(CCFLAGS = '-fPIC')
#end


## Load the configuration file if it exists. The configuration file
## is a python script that updates the env variable with different
## paths.
##
env.AppendUnique(LIBPATH = ['#external/build/lib'])
env.AppendUnique(LIBPATH = ['#sources/build'])
env.AppendUnique(CPPPATH = ['#external/build/include'])
env.AppendUnique(CPPPATH = ['#sources'])



## Launch the compilations
##
Export('env')

external = env.SConscript('external/SConscript')

core    = env.SConscript('sources/core/SConscript')
plugins = env.SConscript('sources/plugins/SConscript')
softs   = env.SConscript('sources/softs/SConscript')
#sources  = env.SConscript('sources/SConscript')
#env.Depends(sources, external)

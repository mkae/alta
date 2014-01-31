import sys

env = Environment()

## PLATFORM dependant section
##

if sys.platform == 'darwin':

	# Adding the /usr/local/lib directory used to store libraries of
	# MacPorts or Brew.
	env.AppendUnique(LIBPATH = '/usr/local/lib')
	env.AppendUnique(CPPPATH = '/usr/local/include')

	env.AppendUnique(LIBPATH =  '/opt/local/lib/')
	env.AppendUnique(CPPPATH =  '/opt/local/include/')

#end

## COMPILER dependant section
##

if env['CC'] == ['gcc', 'g++', 'clang']:
	env.AppendUnique(CCFLAGS = '-fPIC')
#end




Export('env')

external = env.SConscript('external/SConscript')
sources  = env.SConscript('sources/SConscript')

env.Depends(sources, external)

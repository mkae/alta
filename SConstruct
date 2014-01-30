env = Environment()

if env['CC'] == 'gcc' or env['CC'] == 'clang':
	env.Append(CCFLAGS = '-fPIC')
#end

Export('env')

external = env.SConscript('external/SConscript')
sources  = env.SConscript('sources/SConscript')

env.Depends(sources, external)

env = Environment()

external = env.SConscript('external/SConscript')
sources  = env.SConscript('sources/SConscript')

env.Depends(sources, external)

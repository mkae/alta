TEMPLATE = subdirs
SUBDIRS  = core    \
           plugins \
           tests

plugins.depends = core
tests.depends   = core
tests.depends   = plugins

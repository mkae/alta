TEMPLATE = subdirs
SUBDIRS  = core    \
           plugins \
           softs

plugins.depends = core
softs.depends   = core
softs.depends   = plugins

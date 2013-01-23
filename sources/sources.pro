TEMPLATE = subdirs
SUBDIRS  = plugins \
           tests

tests.depends = plugins

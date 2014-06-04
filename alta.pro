TEMPLATE = subdirs
SUBDIRS  = sources             \
           external/quadprog++-v2

sources.depends = external/quadprog++-v2

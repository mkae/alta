TEMPLATE = subdirs
SUBDIRS  = sources             \
           external/quadprog++

sources.depends = external/quadprog++

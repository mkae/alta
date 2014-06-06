CONFIG -= qt

TEMPLATE = subdirs
SUBDIRS  = sources                \
           external/quadprog++    \
           external/quadprog++-v2

sources.depends = external/quadprog++
sources.depends = external/quadprog++-v2

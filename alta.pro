CONFIG -= qt

TEMPLATE = subdirs
SUBDIRS  = sources                \
           external               \

sources.depends = external

#!/bin/sh

# Copy all the file to the ALTA-alpha directory
hg archive ../alta-alpha
cd ../alta-alpha

# Remove unecessary parts (data, retro plugins, ...)
rm -rf data
rm -rf papers
rm -rf sources/tests
rm -rf sources/matlab
rm -rf sources/xml
rm -rf sources/scripts/unitary*
rm -rf sources/plugins/*retro*
rm -rf sources/plugins/shifted_gamma_function
rm -rf sources/softs/tests
rm -rf sources/softs/rational_1d
rm -rf sources/softs/generate_data
rm -rf sources/softs/fourieranalysis
rm -rf sources/softs/data2diff

# Creating the C++ LICENSE header
echo "/*\n$(cat LICENSE.txt | while read LINE; do echo " * $LINE"; done;)\n */\n\n" > LICENSE.C

# Add headers to all .cpp .h files
find ./sources/ -regex ".*\.\(h\|cpp\)$" | while read FILE; do
	cat LICENSE.C $FILE > $FILE.temp;
	mv $FILE.temp $FILE;
done
rm -rf LICENSE.C

# Create the archive
cd ../
tar -cf alta-alpha.tar alta-alpha
gz alta-alpha.tar 

# Upload the archive
scp alta-alpha.tar.gz belcour@scm.gforge.inria.fr:/home/groups/alta/htdocs/downloads

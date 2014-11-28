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
rm -rf sources/plugins/shifted_gamma_function
rm -rf sources/softs/tests
rm -rf sources/softs/rational_1d
rm -rf sources/softs/generate_data
rm -rf sources/softs/fourieranalysis
rm -rf sources/softs/data2diff
rm -rf sources/scripts/generate-alpha.sh

# Creating the C++ LICENSE header
echo "/*$(cat LICENSE.txt | while read LINE; do echo " * $LINE"; done;) */" > LICENSE.C

# Add headers to all .cpp .h files
find ./sources/ -regex ".*\.\(h\|cpp\)$" | while read FILE; do
	cat LICENSE.C $FILE > $FILE.temp;
	mv $FILE.temp $FILE;
done
#rm -rf LICENSE.C


# Create the archive
cd ../
tar -cf alta-alpha.tar alta-alpha
gzip alta-alpha.tar 

# Copy file that contains not reference to retro
#cp ../alta-alpha-temp/sources/plugins/SConscript ./sources/plugins/
#cp ../alta-alpha-temp/sources/core/params.* ./sources/core/
#cp ../alta-alpha-temp/sources/plugins/nonlinear_function_beckmann/function.* ./sources/plugins/nonlinear_function_beckmann/
#cp ../alta-alpha-temp/sources/plugins/nonlinear_function_spherical_gaussian//function.* ./sources/plugins/nonlinear_function_spherical_gaussian/

# Upload the archive
#scp alta-alpha.tar.gz belcour@scm.gforge.inria.fr:/home/groups/alta/htdocs/downloads

/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */


/* Test whether the 'load' and 'save' methods of 'data' work as expected.  We
 * do that by bitwise comparison of the files produced by 'save', which is
 * not ideal (we may want to compare the in-memory representations
 * instead.)  */

#include <core/args.h>
#include <core/data.h>
#include <core/data_storage.h>
#include <core/vertical_segment.h>

#include <string>
#include <iostream>

#include <cstring>
#include <cstdlib>

// Get the 'unlink' declaration.
#ifdef _WIN32
# include <io.h>
# define unlink _unlink
#else
# include <unistd.h>
#endif

// Files that are automatically deleted upon destruction.
class temporary_file
{
public:
		temporary_file()
		{
				static std::string suffix;
				_name = "t-data-io" + suffix;
				suffix += "-";
		}

		~temporary_file()
		{
				::unlink(_name.c_str());
		}

		const std::string& name() { return _name; };
		operator const std::string& () { return _name; }

private:
		std::string _name;
};

// Try hard to read N bytes from INPUT into BUF.
static std::streamsize read_exactly(std::istream &input, char *buf, std::streamsize n)
{
		std::streamsize total;

		for (total = 0;
				 total < n && !input.eof();
				 total += input.gcount())
		{
				input.read(&buf[total], n - total);
		}

		return total;
}

// Return true if the contents of FILE1 are identical to the contents of
// FILE2.
static bool files_are_equal(const std::string &file1, const std::string &file2)
{
		std::ifstream i1, i2;

		i1.exceptions(std::ios_base::failbit);
		i1.open(file1.c_str());
		i2.open(file2.c_str());
		i1.exceptions(std::ios_base::goodbit);

		while (!i1.eof() && !i2.eof())
		{
				char buf1[16384], buf2[16384];
				std::streamsize len1 = read_exactly(i1, buf1, sizeof buf1);
				std::streamsize len2 = read_exactly(i2, buf2, sizeof buf2);
				if (len1 != len2 || memcmp(buf1, buf2, len1) != 0) {
						return false;
				}
		}

		return i1.eof() && (i1.eof() == i2.eof());
}

int main(int argc, char** argv)
{
		std::string input_file;
		temporary_file temp_file1, temp_file2, temp_file3;
		vertical_segment sample1, sample2, sample3;

		if (argc > 1)
				// Process the user-specified file.
				input_file = argv[1];
		else
		{
				// Process the default test file.
				static const std::string data_file = "Kirby2.dat";
				std::string data_dir = getenv("TEST_DATA_DIRECTORY") != NULL
						? getenv("TEST_DATA_DIRECTORY") : ".";
				input_file = data_dir + "/" + data_file;
		}

		try
		{
				// Use the standard load/save methods.
				sample1.load(input_file);
				sample1.save(temp_file1);

				sample2.load(temp_file1);
				sample2.save(temp_file2);

				// Now use the binary output format.
				std::ofstream out;
				out.open(temp_file3.name().c_str());
				save_data_as_binary(out, sample2);
				out.close();

				std::ifstream in;
				in.open(temp_file3.name().c_str());
				load_data_from_binary(in, sample3);
				in.close();
		}
		CATCH_FILE_IO_ERROR(input_file);

		return (sample1.equals(sample2)
						&& files_are_equal(temp_file1, temp_file2)
				    && sample2.equals(sample3))
				? EXIT_SUCCESS : EXIT_FAILURE;
}

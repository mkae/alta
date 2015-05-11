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
#include <core/vertical_segment.h>

#include <string>
#include <iostream>

#include <cstring>
#include <cstdlib>

static void make_temp_file_name(std::string &result)
{
		static std::string suffix;
		result = "t-data-io" + suffix;
		suffix += "-";
}

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
		std::string input_file, temp_file1, temp_file2;
		vertical_segment sample1, sample2;

		make_temp_file_name(temp_file1);
		make_temp_file_name(temp_file2);

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
		}
		CATCH_FILE_IO_ERROR(input_file);

		return files_are_equal(temp_file1, temp_file2)
				? EXIT_SUCCESS : EXIT_FAILURE;
}

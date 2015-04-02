/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <iostream>
#include "data.h"
#include "data_storage.h"

void data::save(const std::string& filename) const
{
		std::ofstream file;

		file.exceptions(std::ios_base::failbit);
		file.open(filename.c_str(), std::ios_base::trunc);
		file.exceptions(std::ios_base::goodbit);

		save_data_as_text(file, *this);

		file.close();
}

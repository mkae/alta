/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

#include <iostream>
#include "data.h"
#include "common.h"
#include "args.h"

namespace alta
{
	// Write DATA to OUT in ALTA's text format.
	void save_data_as_text(std::ostream& out, const alta::data &data);

	// Write DATA to OUT in a compact binary format.
	void save_data_as_binary(std::ostream& out, const alta::data& data);

	// Initialize DATA from the binary-formatted stream IN.
	void load_data_from_binary(std::istream& in, const alta::arguments& header,
	                           alta::data &data);	
}


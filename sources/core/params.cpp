#include "params.h"

// Assing the input params map
const std::map<params::input, const params::param_info> params::input_map = {
	{COS_TH,                {"COS_TH",                1, "Cosine of the Half angle"}},
	{RUSIN_TH_TD,           {"RUSIN_TH_TD",           2, "Radialy symmetric Half angle parametrization"}},
	{RUSIN_TH_TD_PD,        {"RUSIN_TH_TD_PD",        3, "Isotropic Half angle parametrization"}},
	{RUSIN_TH_PH_TD_PD,     {"RUSIN_TH_PH_TD_PD",     4, "Complete Half angle parametrization"}},
	{SPHERICAL_TL_PL_TV_PV, {"SPHERICAL_TL_PL_TV_PV", 4, "Complete classical parametrization"}},
	{CARTESIAN,             {"CARTESIAN",             6, "Complete vector parametrization"}}
};

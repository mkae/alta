#include "params.h"

struct param_info
{
	param_info(std::string n, int d, std::string i) : 
		name(n), dimension(d), info(i) { };

	std::string name;
	int dimension;
	std::string info;
};

// Assing the input params map
static const std::map<params::input, const param_info> input_map = {
	{params::COS_TH,                {"COS_TH",                1, "Cosine of the Half angle"}},
	{params::RUSIN_TH_TD,           {"RUSIN_TH_TD",           2, "Radialy symmetric Half angle parametrization"}},
	{params::RUSIN_TH_TD_PD,        {"RUSIN_TH_TD_PD",        3, "Isotropic Half angle parametrization"}},
	{params::RUSIN_TH_PH_TD_PD,     {"RUSIN_TH_PH_TD_PD",     4, "Complete Half angle parametrization"}},
	{params::SPHERICAL_TL_PL_TV_PV, {"SPHERICAL_TL_PL_TV_PV", 4, "Complete classical parametrization"}},
	{params::CARTESIAN,             {"CARTESIAN",             6, "Complete vector parametrization"}}
};

params::input params::parse_input(const std::string& txt)
{
	for(std::map<params::input, const param_info>::const_iterator it=input_map.begin(); it != input_map.end(); ++it)
	{
		if(txt.compare(it->second.name) == 0)
		{
			return it->first;
		}
	}

	return params::UNKNOWN_INPUT;

	/*
		if(txt == std::string("COS_TH"))
		{
		return params::COS_TH;
		}
		else if(txt == std::string("RUSIN_TH_TD"))
		{
		return params::RUSIN_TH_TD;
		}
		else if(txt == std::string("RUSIN_TH_PH_TD_PD"))
		{
		return params::RUSIN_TH_PH_TD_PD;
		}
		else
		{
		return params::UNKNOWN_INPUT;
		}
		*/
}
		  
std::string params::get_name(const params::input param)
{
	std::map<params::input, const param_info>::const_iterator it = input_map.find(param);
	if(it != input_map.end())
	{
		return it->second.name;
	}

	return std::string();
}

int  params::dimension(params::input t)
{
	std::map<params::input, const param_info>::const_iterator it = input_map.find(t);
	if(it != input_map.end())
	{
		return it->second.dimension;
	}
	else
	{
		return -1;
	}
	/*
		switch(t)
		{
	// 1D Parametrizations
	case params::COS_TH:
	return 1;
	break;

	// 2D Parametrizations
	case params::ISOTROPIC_TD_PD:
	case params::RUSIN_TH_TD:
	case params::ROMEIRO_TH_TD:
	case params::COS_TH_TD:
	return 2;
	break;

	// 3D Parametrization
	case params::RUSIN_TH_PH_TD:
	case params::RUSIN_TH_TD_PD:
	case params::ISOTROPIC_TV_TL_DPHI:
	return 3;
	break;

	// 4D Parametrization
	case params::RUSIN_TH_PH_TD_PD:
	case params::SPHERICAL_TL_PL_TV_PV:
	return 4;
	break;

	// 6D Parametrization
	case params::CARTESIAN:
	return 6;
	break;

	default:
	assert(false);
	return -1;
	break;
	}
	*/
}

void params::print_input_params()
{
	for(std::map<params::input, const param_info>::const_iterator it=input_map.begin(); it != input_map.end(); ++it)
	{
		std::cout << it->second.name << std::endl;
	}
}

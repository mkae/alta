#include "params.h"
#include "common.h"

struct param_info
{
	param_info(std::string n, int d, std::string i) : 
        name(n), dimension(d), info(i) { }

	std::string name;
	int dimension;
	std::string info;
};

//#ifdef WIN32
std::map<params::input, const param_info> create_map()
{
	std::map<params::input, const param_info> _map;
	/* 1D Params */
	_map.insert(std::make_pair<params::input, const param_info>(params::COS_TH, param_info("COS_TH", 1, "Cosine of the Half angle")));
	_map.insert(std::make_pair<params::input, const param_info>(params::COS_TK, param_info("COS_TK", 1, "Cosine of the Back angle")));

	/* 2D Params */
	_map.insert(std::make_pair<params::input, const param_info>(params::RUSIN_TH_TD, param_info("RUSIN_TH_TD", 2, "Radialy symmetric Half angle parametrization")));
	_map.insert(std::make_pair<params::input, const param_info>(params::ISOTROPIC_TV_PROJ_DPHI, param_info("ISOTROPIC_TV_PROJ_DPHI", 2, "Isoptropic projected phi parametrization without a light direction.")));

	/* 3D Params */
	_map.insert(std::make_pair<params::input, const param_info>(params::RUSIN_TH_TD_PD, param_info("RUSIN_TH_TD_PD", 3, "Isotropic Half angle parametrization")));
	_map.insert(std::make_pair<params::input, const param_info>(params::ISOTROPIC_TV_TL_DPHI, param_info("ISOTROPIC_TV_TL_DPHI", 3, "Isotropic Light/View angle parametrization")));
	_map.insert(std::make_pair<params::input, const param_info>(params::RUSIN_VH, param_info("RUSIN_VH", 3, "Vector representation of the Half angle only")));
    _map.insert(std::make_pair<params::input, const param_info>(params::SCHLICK_VK, param_info("SCHLICK_VH", 3, "Vector representation of the Back angle only")));

	/* 4D Params */
	_map.insert(std::make_pair<params::input, const param_info>(params::RUSIN_TH_PH_TD_PD, param_info("RUSIN_TH_PH_TD_PD", 4, "Complete Half angle parametrization")));
	_map.insert(std::make_pair<params::input, const param_info>(params::SPHERICAL_TL_PL_TV_PV, param_info("SPHERICAL_TL_PL_TV_PV", 4, "Complete classical parametrization")));
	_map.insert(std::make_pair<params::input, const param_info>(params::STEREOGRAPHIC, param_info("STEREOGRAPHIC", 4, "Light/View vector in stereographic projection"))),

	/* 6D Param */
	_map.insert(std::make_pair<params::input, const param_info>(params::CARTESIAN, param_info("CARTESIAN", 6, "Complete vector parametrization")));

	return _map;
}
static const std::map<params::input, const param_info> input_map = create_map();
//#else

// Assing the input params map
//static const std::map<params::input, const param_info> input_map = {
	/* 1D Params */
//	{params::COS_TH,                {"COS_TH",                1, "Cosine of the Half angle"}},

	/* 2D Params */
//	{params::RUSIN_TH_TD,           {"RUSIN_TH_TD",           2, "Radialy symmetric Half angle parametrization"}},

	/* 3D Params */
//	{params::RUSIN_TH_TD_PD,        {"RUSIN_TH_TD_PD",        3, "Isotropic Half angle parametrization"}},
//	{params::ISOTROPIC_TV_TL_DPHI,  {"ISOTROPIC_TV_TL_DPHI",  3, "Isotropic Light/View angle parametrization"}},
	
	/* 4D Params */
//	{params::RUSIN_TH_PH_TD_PD,     {"RUSIN_TH_PH_TD_PD",     4, "Complete Half angle parametrization"}},
//	{params::SPHERICAL_TL_PL_TV_PV, {"SPHERICAL_TL_PL_TV_PV", 4, "Complete classical parametrization"}},
//    {params::STEREOGRAPHIC,         {"STEREOGRAPHIC",         4, "Light/View vector in stereographic projection"}},

	/* 6D Params */
//	{params::CARTESIAN,             {"CARTESIAN",             6, "Complete vector parametrization"}}
//};
//#endif

void params::to_cartesian(const double* invec, params::input intype,
		double* outvec)
{
	switch(intype)
	{
		// 1D Parametrizations
		case params::COS_TH:
			half_to_cartesian(acos(invec[0]), 0.0, 0.0, 0.0, outvec);
			break;
		case params::COS_TK:
			outvec[0] = invec[0]*invec[0]-1.0;
			outvec[1] = 0.0;
			outvec[2] = invec[0];
			outvec[3] = 1.0-invec[0]*invec[0];
			outvec[4] = 0.0;
			outvec[5] = invec[0];
			break;

			// 2D Parametrizations
		case params::COS_TH_TD:
			half_to_cartesian(acos(invec[0]), 0.0, acos(invec[1]), 0.0, outvec);
			break;
		case params::RUSIN_TH_TD:
         half_to_cartesian(invec[0], 0.0, invec[1], 0.5*M_PI, outvec);
			break;
		case ISOTROPIC_TV_PROJ_DPHI:
		{
			const double theta = 0.5*sqrt(invec[0]*invec[0] + invec[1]*invec[1]);
			outvec[0] = invec[0]/theta*sin(theta);
			outvec[1] = invec[1]/theta*sin(theta);
			outvec[2] = cos(theta);
			outvec[3] = 0.0;
			outvec[4] = 0.0;
			outvec[5] = 1.0;
		}
			break;

			// 3D Parametrization
		case params::RUSIN_TH_PH_TD:
			half_to_cartesian(invec[0], invec[1], invec[2], 0.0, outvec);
			break;
		case params::RUSIN_TH_TD_PD:
			half_to_cartesian(invec[0], 0.0, invec[1], invec[2], outvec);
#ifdef DEBUG
			std::cout << outvec[2] << std::endl;
#endif
			break;
		case params::ISOTROPIC_TV_TL_DPHI:
			classical_to_cartesian(invec[0], 0.0, invec[1], invec[2], outvec);
			break;
		case params::RUSIN_VH:
            half_to_cartesian(acos(invec[2]), atan2(invec[1], invec[0]), 0.0, 0.0, outvec);
			break;
        // \todo I should handle the phi_k in the conversion to CARTESIAN
        case params::SCHLICK_VK:
            outvec[0] = invec[2]*invec[2]-1.0;
            outvec[1] = 0.0;
            outvec[2] = invec[2];
            outvec[3] = 1.0-invec[2]*invec[2];
            outvec[4] = 0.0;
            outvec[5] = invec[2];
            break;

			// 4D Parametrization
		case params::RUSIN_TH_PH_TD_PD:
			half_to_cartesian(invec[0], invec[1], invec[2], invec[3], outvec);
			break;

		case params::SPHERICAL_TL_PL_TV_PV:
			outvec[0] = cos(invec[1])*sin(invec[0]);
			outvec[1] = sin(invec[1])*sin(invec[0]);
			outvec[2] = cos(invec[0]);
			outvec[3] = cos(invec[3])*sin(invec[2]);
			outvec[4] = sin(invec[3])*sin(invec[2]);
			outvec[5] = cos(invec[2]);
			break;

		case params::STEREOGRAPHIC:
      {
			// Project the 2D direction on the surface invec[0,1] to the View
         // vector outvec[0,1,2]
         const double normL = invec[0]*invec[0] + invec[1]*invec[1];
         const double alphL = (1.0 - normL) / (1.0 + normL);
         outvec[0] = invec[0] + alphL;
         outvec[1] = invec[1] + alphL;
         outvec[2] = alphL;

         // Project the 2D direction on the surface invec[2,3] to the Light
         // vector outvec[3,4,5]
         const double normV = invec[2]*invec[2] + invec[3]*invec[3];
         const double alphV = (1.0 - normV) / (1.0 + normV);
         outvec[3] = invec[2] + alphV;
         outvec[4] = invec[3] + alphV;
         outvec[5] = alphV;
		}
         break;

		// 6D Parametrization
      case params::CARTESIAN:
			memcpy(outvec, invec, 6*sizeof(double));
			break;

		default:
			std::cerr << "<<ERROR>> Transformation not implemented, " << get_name(intype) << " " << __FILE__ << ":" << __LINE__ << std::endl;
			throw;
			break;
	}

}

void params::from_cartesian(const double* invec, params::input outtype,
		double* outvec)
{
	// Compute the half vector
	double half[3] ;
	half[0] = invec[0] + invec[3];
	half[1] = invec[1] + invec[4];
	half[2] = invec[2] + invec[5];
	double half_norm = sqrt(half[0]*half[0] + half[1]*half[1] + half[2]*half[2]);
	half[0] /= half_norm;
	half[1] /= half_norm;
	half[2] /= half_norm;

	// Difference vector 
	double diff[3];

	switch(outtype)
	{
		// 1D Parametrizations
		case params::COS_TH:
			outvec[0] = half[2];
			break;
		case params::COS_TK:
		{
			const double Kx = invec[0]-invec[3];
			const double Ky = invec[1]-invec[4];
			const double Kz = invec[2]+invec[5];
			outvec[0] = (invec[2] + invec[5]) / sqrt(Kx*Kx + Ky*Ky + Kz*Kz);
		}
			break;

			// 2D Parametrizations
		case params::COS_TH_TD:
			outvec[0] = half[2];
			outvec[1] = half[0]*invec[0] + half[1]*invec[1] + half[2]*invec[2];
			break;
		case params::RUSIN_TH_TD:
			outvec[0] = acos(half[2]);
			outvec[1] = acos(clamp(half[0]*invec[0] + half[1]*invec[1] + half[2]*invec[2], -1.0, 1.0));
			break;
		case ISOTROPIC_TV_PROJ_DPHI:
		{
			const double theta = acos(invec[2]);
			const double dphi  = atan2(invec[1], invec[0]) - atan2(invec[4], invec[3]);
			outvec[0] = theta * cos(dphi);
			outvec[1] = theta * sin(dphi);
		}
			break;

			// 3D Parametrization
		case params::RUSIN_TH_PH_TD:
			outvec[0] = acos(half[2]);
			outvec[1] = atan2(half[1], half[0]);
			outvec[2] = acos(half[0]*invec[0] + half[1]*invec[1] + half[2]*invec[2]);
			break;
		case params::RUSIN_TH_TD_PD:
			outvec[0] = acos(half[2]);

			// Compute the diff vector
			diff[0] = invec[0];
			diff[1] = invec[1];
			diff[2] = invec[2];

			rotate_normal(diff, -atan2(half[1], half[0]));
			rotate_binormal(diff, -outvec[0]);

			outvec[1] = acos(diff[2]);
			outvec[2] = atan2(diff[1], diff[0]);
			break;
		case params::ISOTROPIC_TV_TL_DPHI:
			outvec[0] = acos(invec[2]);
			outvec[1] = acos(invec[5]);
			outvec[2] = atan2(invec[1], invec[0]) - atan2(invec[4], invec[3]);
			break;
		case params::RUSIN_VH:
			outvec[0] = half[0];  
			outvec[1] = half[1];  
			outvec[2] = half[2];  
			break;
      case params::SCHLICK_VK:
      {
			const double Kx = invec[0]-invec[3];
         const double Ky = invec[1]-invec[4];
         const double Kz = invec[2]-invec[5];

         const double norm =  sqrt(Kx*Kx + Ky*Ky + Kz*Kz);
			if(norm > 1.0E-10)
			{
				outvec[0] = Kx / norm;
		      outvec[1] = Ky / norm;
			   outvec[2] = Kz / norm;
			}
			else
			{
				outvec[0] = 0.0;
		      outvec[1] = 0.0;
			   outvec[2] = 1.0;
			}
      }
			break;

			// 4D Parametrization
		case params::RUSIN_TH_PH_TD_PD:
			outvec[0] = acos(half[2]);
			outvec[1] = atan2(half[0], half[1]);

			// Compute the diff vector
			diff[0] = invec[0];
			diff[1] = invec[1];
			diff[2] = invec[2];
			rotate_normal(diff, -atan2(half[1], half[0]));
			rotate_binormal(diff, -outvec[0]);

			outvec[2] = acos(diff[2]);
			outvec[3] = atan2(diff[1], diff[0]);
			break;

		case params::SPHERICAL_TL_PL_TV_PV:
			outvec[0] = acos(invec[2]);
			outvec[1] = atan2(invec[1], invec[0]);
			outvec[2] = acos(invec[5]);
			outvec[3] = atan2(invec[4], invec[3]);
#ifdef DEBUG
			std::cout << invec[2] << " - acos -> " << outvec[0] << std::endl;
#endif
			break;

      case params::STEREOGRAPHIC:
      {
			// Project the View vector invec[0,1,2] on a 2D direction on the
         // surface outvec[0,1]
         const double dotVN = invec[2];
         outvec[0] = invec[0] / (1.0+dotVN);
         outvec[1] = invec[1] / (1.0+dotVN);

         // Project the Light vector invec[0,1,2] on a 2D direction on the
         // surface outvec[2,3]
         const double dotLN = invec[5];
         outvec[2] = invec[3] / (1.0+dotLN);
         outvec[3] = invec[4] / (1.0+dotLN);

		}
         break;

			// 6D Parametrization
		case params::CARTESIAN:
			memcpy(outvec, invec, 6*sizeof(double));
			break;

		default:
			std::cerr << "<<ERROR>> Transformation not implemented, " << get_name(outtype) << " " << __FILE__ << ":" << __LINE__ << std::endl;
			assert(false);
			break;
	}
}
params::input params::parse_input(const std::string& txt)
{
	for(std::map<params::input, const param_info>::const_iterator it=input_map.begin(); it != input_map.end(); ++it)
	{
		if(txt.compare(it->second.name) == 0)
		{
			return it->first;
		}
	}

	std::cout << "<<INFO>> the input parametrization is UNKNOWN_INPUT" << std::endl;
	return params::UNKNOWN_INPUT;
}
		  
std::string params::get_name(const params::input param)
{
	std::map<params::input, const param_info>::const_iterator it = input_map.find(param);
	if(it != input_map.end())
	{
		return it->second.name;
	}

	std::cerr << "<<ERROR>> Unknown parametrization, " << get_name(param) << " " << __FILE__ << ":" << __LINE__ << std::endl;
	return std::string("UNKNOWN_INPUT");
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
}

void params::print_input_params()
{
	for(std::map<params::input, const param_info>::const_iterator it=input_map.begin(); it != input_map.end(); ++it)
	{
		std::cout << it->second.name << std::endl;
	}
}

// Acces to the domain of definition of the function
void parametrized::setMin(const vec& min)
{
#ifdef DEBUG
    assert(min.size() == _nX) ;
#endif
    _min = min ;
}
void parametrized::setMax(const vec& max)
{
#ifdef DEBUG
    assert(max.size() == _nX) ;
#endif
    _max = max ;
}
vec parametrized::min() const
{
    return _min ;
}
vec parametrized::max() const
{
    return _max ;
}

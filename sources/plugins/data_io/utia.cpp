/***************************************************************************

  This plugin provides an input/output data plugin for the ALTA library to
  handle UTIA anisotropic BRDF dataset. This dataset can be obtain at:

       http://btf.utia.cas.cz/

  This is an adaptation of the TBDRF.cpp file from Filip Jiri to the ALTA
  library. Conversion was made by Laurent Belcour on February 2015.

  Authors: 
*****************************************************************************/

#include <core/common.h>
#include <core/data.h>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "EXR_IO.h"

class UTIA : public data {
private:

	float   step_t,step_p;
	int     nti, ntv, npi, npv, planes;
	int     nPerPlane;

	// Database values
	double* Bd;

public:
	UTIA() : data() {
		this->step_t = 15;
		this->step_p = 7.5;
		this->nti = 6;
		this->ntv = 6;
		this->npi = (int)(360.f/step_p);
		this->npv = (int)(360.f/step_p);
		this->planes = 3;
		this->nPerPlane = nti*npi*ntv*npv;
		this->Bd = new double[planes*nti*npi*ntv*npv];

		// Set the input and output parametrization
	    _in_param  = params::SPHERICAL_TL_PL_TV_PV;
	    _out_param = params::RGB_COLOR;
	    _nX = 4;
	    _nY = 3;
	}

	virtual ~UTIA() {
		delete[] this->Bd;
	}

	// Load data from a file
	virtual void load(const std::string& filename) {

		/* If the file is an OpenEXR image */
		if(filename.substr(filename.find_last_of(".") + 1) == "exr") {
			double* temp;
			int W, H;
			if(!t_EXR_IO<double>::LoadEXR(filename.c_str(), W, H, temp, 3) || W != npi*nti || H != ntv*npv) {
				std::cerr << "<<ERROR>> Unable to open file '" << filename << "'" << std::endl;
				throw;

			} else {
				/* Data copy */
				for(int i=0; i<H; ++i)
					for(int j=0; j<W; ++j){
						int index = i*W+j;
						// TODO: The indexing is not the same here
		    			Bd[index + 0*nPerPlane] = temp[3*index + 0];
		    			Bd[index + 1*nPerPlane] = temp[3*index + 1];
		    			Bd[index + 2*nPerPlane] = temp[3*index + 2];
		  			}

				delete[] temp;
			}

		/* If the file is a binary */
		} else {

			std::ifstream stream(filename, std::ios_base::in | std::ios_base::binary);
			/*
			if (!stream.is_open()) {
			  std::cerr << "<<ERROR>> Unable to open file '" << filename << "' !" << std::endl;
			  throw;
			}*/

			int count = 0;
			for(int isp=0;isp<planes;isp++)	{
			  for(int ni=0;ni<nti*npi;ni++)
			    for(int nv=0;nv<ntv*npv;nv++) {
			    	stream >> Bd[count++];
			    }
			}
		}

		std::cout << "<<INFO>> Reading BRDF from file '" << filename << "' ...done" << std::endl;
	}
	virtual void load(const std::string& filename, const arguments&) {
		load(filename);
	}

	virtual void save(const std::string& filename) const {
		/* If the file is an OpenEXR image */
		if(filename.substr(filename.find_last_of(".") + 1) == "exr") {
			int W = npi*nti, H = npv*ntv;
			double* temp = new double[W*H*3];

			/* Data copy */
			for(int i=0; i<H; ++i)
				for(int j=0; j<W; ++j){
					int index = i*W+j;
					// TODO: The indexing is not the same here
	    			temp[3*index + 0] = Bd[index + 0*nPerPlane];
	    			temp[3*index + 1] = Bd[index + 1*nPerPlane];
	    			temp[3*index + 2] = Bd[index + 2*nPerPlane];
	  			}
	  		t_EXR_IO<double>::SaveEXR(filename.c_str(), W, H, temp, 3);
			delete[] temp;
		} else {
			std::ofstream stream(filename, std::ios_base::out | std::ios_base::binary);
			int count = 0;
			for(int isp=0;isp<planes;isp++)	{
			  for(int ni=0;ni<nti*npi;ni++)
			    for(int nv=0;nv<ntv*npv;nv++) {
			    	stream << Bd[count++];
			    }
			}
		}
	}

	// Acces to data
	virtual vec get(int i) const {
		vec res(7);
		res[3] = 2.0*M_PI*double(i%npv)/double(npv); i /= npv;
		res[2] = 0.5*M_PI*double(i%ntv)/double(ntv); i /= ntv;
		res[1] = 2.0*M_PI*double(i%npi)/double(npi); i /= npi;
		res[0] = 0.5*M_PI*double(i%nti)/double(nti);
		res[4] = Bd[i + 0*nPerPlane];
		res[5] = Bd[i + 1*nPerPlane];
		res[6] = Bd[i + 2*nPerPlane];
		return res;
	}
	virtual vec operator[](int i) const {
		vec res(7);
		res[3] = 2.0*M_PI*double(i%npv)/double(npv); i /= npv;
		res[2] = 0.5*M_PI*double(i%ntv)/double(ntv); i /= ntv;
		res[1] = 2.0*M_PI*double(i%npi)/double(npi); i /= npi;
		res[0] = 0.5*M_PI*double(i%nti)/double(nti);
		res[4] = Bd[i + 0*nPerPlane];
		res[5] = Bd[i + 1*nPerPlane];
		res[6] = Bd[i + 2*nPerPlane];
		return res;
	}

	virtual vec value(const vec& in) const {
		// Input and Ouput parameters
		vec RGB(3);
		double theta_i = in[0];
		double phi_i   = in[1];
		double theta_v = in[2];
		double phi_v   = in[3];

		double PI2 = M_PI*0.5;
		if(theta_i>PI2 || theta_v>PI2) {
			RGB[0] = 0.f;
			RGB[1] = 0.f;
			RGB[2] = 0.f;
			return RGB;
		}

		float d2r = 180.f/M_PI;
		theta_i *= d2r;
		theta_v *= d2r;
		phi_i *= d2r;
		phi_v *= d2r;
		if(phi_i>=360.f)
		phi_i = 0.f;
		if(phi_v>=360.f)
		phi_v = 0.f;

		int iti[2],itv[2],ipi[2],ipv[2];
		iti[0] = (int)(floor(theta_i/step_t));
		iti[1] = iti[0]+1;
		if(iti[0]>nti-2) {
		  iti[0] = nti-2;
		  iti[1] = nti-1;
		}
		itv[0] = (int)(floor(theta_v/step_t));
		itv[1] = itv[0]+1;
		if(itv[0]>ntv-2) {
		  itv[0] = ntv-2;
		  itv[1] = ntv-1;
		}

		ipi[0] = (int)(floor(phi_i/step_p));
		ipi[1] = ipi[0]+1;
		ipv[0] = (int)(floor(phi_v/step_p));
		ipv[1] = ipv[0]+1;

		float sum;
		float wti[2],wtv[2],wpi[2],wpv[2];
		wti[1] = theta_i - (float)(step_t*iti[0]);
		wti[0] = (float)(step_t*iti[1]) - theta_i;
		sum = wti[0]+wti[1];
		wti[0] /= sum;
		wti[1] /= sum;
		wtv[1] = theta_v - (float)(step_t*itv[0]);
		wtv[0] = (float)(step_t*itv[1]) - theta_v;
		sum = wtv[0]+wtv[1];
		wtv[0] /= sum;
		wtv[1] /= sum;

		wpi[1] = phi_i - (float)(step_p*ipi[0]);
		wpi[0] = (float)(step_p*ipi[1]) - phi_i;
		sum = wpi[0]+wpi[1];
		wpi[0] /= sum;
		wpi[1] /= sum;
		wpv[1] = phi_v - (float)(step_p*ipv[0]);
		wpv[0] = (float)(step_p*ipv[1]) - phi_v;
		sum = wpv[0]+wpv[1];
		wpv[0] /= sum;
		wpv[1] /= sum;

		if(ipi[1]==npi) 
		ipi[1] = 0;
		if(ipv[1]==npv) 
		ipv[1] = 0;

		int nc = npv*ntv;
		int nr = npi*nti;
		for(int isp=0;isp<planes;isp++)	{
		  RGB[isp] = 0.f;
		  for(int i=0;i<2;i++)
		    for(int j=0;j<2;j++)
		      for(int k=0;k<2;k++)
		        for(int l=0;l<2;l++)
		          RGB[isp] += Bd[isp*nr*nc + nc*(npi*iti[i]+ipi[k]) + npv*itv[j]+ipv[l]] * wti[i] * wtv[j] * wpi[k] * wpv[l];
		  //      RGB[isp] *= cos(theta_i/d2r);
		}
		return RGB;
	}

	// Set data
	virtual void set(const vec& x) {
		assert(x.size() == dimX()+dimY());
	}

	virtual void set(int i, const vec& x) {
		assert(x.size() == dimY());
		for(int isp=0; isp<planes; ++isp) {
			Bd[isp*nPerPlane + i] = x[isp];
		}
	}

	// Get data size, e.g. the number of samples to fit
	virtual int size() const {
		return nPerPlane;
	}
};

ALTA_DLL_EXPORT data* provide_data()
{
    return new UTIA();
}

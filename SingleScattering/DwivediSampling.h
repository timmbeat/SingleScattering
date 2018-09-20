/*
 * Created by Tim Mend on 12.09.2018
 */



#pragma once
#include "Sampling.h"
class DwivediSampling : public Sampling
{
	public:
	DwivediSampling(double absorption, double scattering, double anisotropy, double diameter, double delr, std::size_t runs);
	~DwivediSampling();
	//Sampling Specific functions
	double dirdis(double const wz) const;
	double pathdis(double const wz, double const t) const;
	double samplePathDistribution(double wz) const;
	double sampleDirDistribution() const ;
	double calculateLr() override;


	
	//Getter
	double V0() const;

	private:
	double const v0;

};


/*
 * Created by Tim Mend on 12.09.2018
 */



#pragma once
#include "Sampling.h"
class DwivediSampling : public Sampling
{
	public:
	DwivediSampling(Real absorption, Real scattering, Real anisotropy, double diameter, double delr, std::size_t runs);
	~DwivediSampling();


	//Sampling Specific functions
	Real dirdis(Real const wz) const;
	Real pathdis(Real const wz, Real const t) const;
	Real samplePathDistribution(Real wz) const;
	Real sampleDirDistribution() const;
	Real calculateLr() override;


	
	//Getter
	Real V0() const;

	private:
	Real const v0;

};


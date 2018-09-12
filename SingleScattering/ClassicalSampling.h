/*
* Created by Tim Mend on 12.09.2018
*/


#pragma once
#include "constants.h"
#include "Sampling.h"
class ClassicalSampling : public Sampling
{
	public:
	ClassicalSampling(Real absorption, Real scattering, Real anisotropy, double diameter, double delr, std::size_t runs);
	~ClassicalSampling();


	Real samplePathDistribution() const;
	Real sampleDirDistribution() const;
	Real calculateLr() override;
};


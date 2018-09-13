/*
* Created by Tim Mend on 12.09.2018
*/


#pragma once
#include "constants.h"
#include "Sampling.h"
class ClassicalSampling : public Sampling
{
	public:
	ClassicalSampling(double absorption, double scattering, double anisotropy, double diameter, double delr, std::size_t runs);
	~ClassicalSampling();


	double samplePathDistribution() const;
	double sampleDirDistribution() const;
	double calculateLr() override;
};


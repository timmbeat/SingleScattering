#pragma once
#include "constants.h"
#include "Sampling.h"
class ClassicalSampling
{
	public:
	ClassicalSampling(float absorption, float scattering, float anisotropy, size_t binsr, float delr);
	~ClassicalSampling();

	double sdirectionaldistribution()
	{
		if (anisotropy == 0) return 2 * random() - 1;
		auto const qani = anisotropy * anisotropy;

		return (1 / (2 * anisotropy)) * (1 + qani - pow(((1 - qani) / (1 - anisotropy + 2 * anisotropy*random())), 2));

	}

	double spathdistribution()
	{
		return -log(random()) / (absorption + scattering);
	}

	double calculatelr();


	decltype(auto) getbinsr() const
	{
		return binsr;
	}
	/*
	 * Returns a reference to the bin
	 */
	decltype(auto) getbin() const
	{
		return &bins;
	}

private:
float absorption;
float scattering;
float anisotropy;
int binsr;
float delr;
std::vector<double> bins;

};


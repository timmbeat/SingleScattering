#pragma once
#include "Sampling.h"
class DwivediSampling
{
	public:
	DwivediSampling(float const absorption, float const scattering, float const anisotropy, float size, float const delr);
	~DwivediSampling();


	double sampledirdis(double const v0) const 
	{
		return v0 - (v0 + 1)*pow((v0 - 1) / (v0 + 1), random());
	}

	double samplepathdis(double const wz, double const v0, double const mut) const
	{
		return (-log(random())) / ((1 - v0 / wz)*mut);
	}

	double dirdis(double const v0, double const wz) const 
	{
		return (1 / log((v0 + 1) / (v0 - 1))) * (1 / (v0 - wz));
	}

	double pathdis(double const v0, double const wz, double const mut, double const t) const
	{
		return (1 - wz / v0)*mut*exp(-(1 - wz / v0)*mut*t);
	}

	double calculateonelr(size_t runs);

	double getv0()
	{
		return v0;
	}


	/*
	 * Returns a reference to the Bins
	 */
	decltype(auto) getBins() const 
	{
		return &bins;
	}

	decltype(auto) getDelr()const
	{
		return delr;
	}
	
	private:
	double const v0;
	float const binsr;
	float const delr;
	float const absorption;
	float const scattering;
	float const anisotropy;
	std::vector<double> bins;

};


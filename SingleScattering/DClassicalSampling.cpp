#include "DClassicalSampling.h"
#include <glm.hpp>
#include <iostream>

DClassicalSampling::DClassicalSampling(glm::dvec3 position, glm::dvec3 direction, double absorption, double scattering, double anisotropy, double diameter, double delr, std::size_t runs, int times, bool forcescattering) :
	position(position), direction(direction), Lr(1.0), direction_org(direction), position_org(position), times(times), forcescattering(forcescattering), scatterevent(0),
	ClassicalSampling{ absorption, scattering, anisotropy, diameter, delr, runs }
{

}


DClassicalSampling::~DClassicalSampling()
{
}

double DClassicalSampling::calculateLr()
{
	//SCHRITT 1
	double const s1 = samplePathDistribution();

	//Check if Boundary is hit.
	if (boundary(s1))
	{
		return Lr;
	}


	auto const cos_theta = ClassicalSampling::sampleDirDistribution();
	updatePosition(s1);
	updateDirection(cos_theta);
	scatterevent++;
	Lr *= (Sampling::tauo_di(s1, Absorption(), Scattering()) / taudi_r(s1, Absorption(), Scattering()))*(henyey_greenstein(acos(cos_theta), Anisotropy()) / henyey_greenstein_norm(acos(cos_theta), Anisotropy()));
}

void DClassicalSampling::reset()
{
	position = position_org;
	direction = direction_org;
	Lr = 1.0;
}

void DClassicalSampling::updateDirection(double const cos_theta)
{
	auto const sint = sqrt(1 - cos_theta * cos_theta);


	auto const uz = direction.z;
	auto const ux = direction.x;
	auto const uzq = uz * uz;

	if (abs(uz) > cos1<double>())
	{
		direction.x = sint;
		direction.y = 0.0;
		direction.z = glm::sign(uz)*cos_theta;
	}
	else
	{
		direction.x = sint / sqrt(1 - uzq)*(ux*uz) + ux * cos_theta;
		direction.y = 0.0;
		direction.z = -sint * sqrt(1 - uzq) + uz * cos_theta;
	}
}

void DClassicalSampling::updatePosition(double const stepsize)
{
	position.x = position.x + direction.x*stepsize;
	position.y = 0.0;
	position.z = position.z + direction.z*stepsize;
}

void DClassicalSampling::out()
{
	if (direction.z <= 0.0)

	{
		Lr = 0.0;
	}

	else
	{

		auto const stepleft = -position.z / direction.z;

		updatePosition(stepleft);


		Lr *= taudi_r(abs(stepleft), Absorption(), Scattering());

		auto binnum = static_cast<size_t>(abs(position.x) / Delr());
		if (binnum >= Binsr()) binnum = Binsr() - 1;


		auto const Lr_sized = static_cast<double>(Lr / Runs());

		(*Bins())[binnum] += Lr_sized;
	}
}

bool DClassicalSampling::boundary(double const stepsize)
{
	if (-position.z / direction.z <= stepsize && direction.z > 0.0)
	{
		if (scatterevent < times && forcescattering)
		{
			Lr = 0.0;
			return true;
		}
		Lr *= taudi_r(abs(-position.z / direction.z), Absorption(), Scattering());
		auto const Lr_sized = static_cast<double>(Lr / Runs());
		updatePosition(-position.z / direction.z);

		auto binnumber = static_cast<size_t>(abs(position.x) / Delr());
		if (binnumber >= Binsr()) binnumber = Binsr() - 1;


		(*Bins())[binnumber] += Lr_sized;
		return true;

	}

	return false;

}

double DClassicalSampling::run()
{
	reset();
	for (auto i = 0; i < times; i++)
	{
		calculateLr();
	}

	out();

	return Lr / Runs();
}

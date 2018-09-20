#include "DDwivediSampling.h"
#include <iostream>


DDwivediSampling::DDwivediSampling(glm::dvec3 position, glm::dvec3 direction, double absorption, double scattering,
								   double anisotropy, double diameter, double delr, std::size_t runs, double wz, int times, bool forcescattering) :
	position(position), direction(direction), wz(wz), wz_old(wz), Lr(1.0), Lr_org(1.0), wz_org(wz), position_org(position), direction_org(direction), times(times), forcescattering(forcescattering), scatterevent(0),
	DwivediSampling{ absorption, scattering, anisotropy, diameter, delr, runs }
{
}

DDwivediSampling::~DDwivediSampling()
{
}

void DDwivediSampling::updateDirection()
{

	auto const sint = sqrt(1 - wz * wz);
	auto const uz = direction.z;
	direction.x = sint;
	direction.y = 0.0;
	direction.z = glm::sign(uz) * wz;
}


void DDwivediSampling::updatePosition(double stepsize)
{
	position.x = position.x + direction.x*stepsize;
	position.y = 0.0;
	position.z = position.z + direction.z*stepsize;
}

double DDwivediSampling::calculateLr()
{
	if (scatterevent == times)
	{
		dead = true;
		return Lr;
	}


	double const s1 = samplePathDistribution(wz_old);


	//Check if Boundary is hit.
	if (boundary(s1))
	{
		return Lr;
	}


	double wz1 = sampleDirDistribution();
	wz = wz1; //Set new wz;
	updatePosition(s1);
	auto direction_old = direction;
	updateDirection();
	scatterevent++;
	auto direction_new = direction;
	auto dot = glm::dot(direction_new, direction_old);
	auto dot_theta = dot / (glm::length(direction_old) * glm::length(direction_new));
	auto dot_norm = glm::dot(direction_new, glm::dvec3(0.0, 0.0, 1.0)) / (glm::length(direction_new)*glm::length(glm::dvec3(0.0, 0.0, 1.0)));

	auto cos_theta = dot_theta < 0.0 ? - 1 - dot_theta : 1 - dot_theta;
	if (dot_theta == 0.0)
	{
		if(direction.z < 0.0)
		{
			cos_theta = -1.0;
		}
		else
		{
			cos_theta = 1.0;
		}
	}
	auto tmp = (henyey_greenstein_norm(cos_theta, Anisotropy()) / dirdis(wz));
	Lr *= (Sampling::tauo_di(s1, Absorption(), Scattering()) / DwivediSampling::pathdis(wz_old, s1))*tmp;

	wz_old = wz;
}



void DDwivediSampling::out()
{
	dead = true;
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

bool DDwivediSampling::boundary(double stepsize)
{
	if (-position.z / direction.z <= stepsize && direction.z > 0.0)
	{
		dead = true;
		if(scatterevent < times && forcescattering)
		{
			Lr = 0.0;
			return true;
		}
		else
		{

			Lr *= taudi_r(abs(-position.z / direction.z), Absorption(), Scattering());
			double const Lr_sized = static_cast<double>(Lr / Runs());
			updatePosition(-position.z / direction.z);

			auto binnumber = static_cast<size_t>(abs(position.x) / Delr());
			if (binnumber >= Binsr()) binnumber = Binsr() - 1;


			(*Bins())[binnumber] += Lr_sized;
			return true;
		}
	}

	return false;
}

double DDwivediSampling::run()
{
	reset();

	while(!dead)
	{
		calculateLr();
	}

	out();

	return Lr/Runs();
}

void DDwivediSampling::reset()
{
	Lr = Lr_org;
	position = position_org;
	direction = direction_org;
	wz = wz_org;
	wz_old = wz_org;
	scatterevent = 0;
	dead = false;
}

double DDwivediSampling::sampleDirDistribution_hen() const
{

	if (Anisotropy() == 0) return 2 * random() - 1;
	auto const qani = Anisotropy() * Anisotropy();

	return (1 / (2 * Anisotropy())) * (1 + qani - pow(((1 - qani) / (1 - Anisotropy() + 2 * Anisotropy()*random())), 2));
}

glm::dvec3 DDwivediSampling::Position() const
{
	return position;
}

glm::dvec3 DDwivediSampling::Direction() const
{
	return direction;
}



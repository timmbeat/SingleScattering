#include "DDwivediSampling.h"
#include <iostream>


DDwivediSampling::DDwivediSampling(glm::dvec3 position, glm::dvec3 direction, double absorption, double scattering,
								   double anisotropy, double diameter, double delr, std::size_t runs) :
	position(position), direction(direction), DwivediSampling{absorption, scattering, anisotropy, diameter, delr, runs}
{
}

DDwivediSampling::~DDwivediSampling()
{
}

void DDwivediSampling::updateDirection(double wz)
{

	//auto const theta = acos(-wz);
	auto const y = 2 * mcss::pi<double>() * random();
	//auto const sint = sqrt(1 - wz * wz);
	auto const sint = sqrt(1 - wz * wz);
	auto const cosp = cos(y);
	double sinp;
	auto const uz = direction.z;
	if (y < mcss::pi<double>())
	{
		sinp = sqrt(1.0 - cosp * cosp);
	}
	else
	{
		sinp = -sqrt(1.0 - cosp * cosp);
	}


	direction.x = sint * sinp;
	direction.y = wz;
}


void DDwivediSampling::updatePosition(double stepsize)
{
	position.x = position.x + direction.x*stepsize;
	position.y = position.y + direction.y*stepsize;
	//position.z = position.z + direction.z*stepsize;
}

double DDwivediSampling::calculateLr()
{
	double Lr = 1.0;


	//SCHRITT 1
	double const wz1 = sampleDirDistribution();
	double const s1 = samplePathDistribution(-1.0);
	updatePosition(s1);
	auto direction_old = direction;
	updateDirection(wz1);
	auto direction_new = direction;
	double dot = glm::dot(direction_new, direction_old);
	double dot_theta = dot / (glm::length(direction_old) * glm::length(direction_new));


	Lr *= (Sampling::tauo_di(s1, Absorption(), Scattering())/DwivediSampling::pathdis(-1.0, s1))*(henyey_greenstein(acos(dot_theta), Anisotropy())/dirdis(wz1));
		

	//SCHRITT 2
	double const s2 = samplePathDistribution(wz1);
	//Check if Boundary is hit
	if (-position.y / direction.y <= s2)
	{
		double const teiler = static_cast<double>(Runs());
		Lr *= taudi_r(abs(-position.y / direction.y), Absorption(), Scattering());
		double const Lr_sized = Lr / teiler;
		updatePosition(-position.y / direction.y);

		auto binnumber = static_cast<size_t>(abs(position.x) / Delr());
		if(binnumber >= Binsr()) binnumber = Binsr() - 1;
		

		(*Bins())[binnumber] += Lr_sized;
		return Lr_sized;

	}
	updatePosition(s2);
	double const wz2 = sampleDirDistribution();
	direction_old = direction;
	updateDirection(wz2);
	direction_new = direction;
	dot = glm::dot(direction_old, direction_new);
	dot_theta = dot / (glm::length(direction_old) * glm::length(direction_new));

	
	
	

	Lr *= (Sampling::tauo_di(s2, Absorption(), Scattering()) / DwivediSampling::pathdis(wz1, s2))*(henyey_greenstein(acos(dot_theta), Anisotropy()) / dirdis(wz2));


	

	if (direction.y < 0) return 0.0;

	double const rest = -position.y / direction.y;


	updatePosition(rest);

	Lr *= taudi_r(abs(rest), Absorption(), Scattering());
	auto test = abs(position.x) / Delr();
	auto binnumber = static_cast<size_t>(abs(position.x) / Delr());
	if (binnumber >= Binsr()) binnumber = Binsr() - 1;

	double const teiler = static_cast<double>(Runs());
	double const Lr_sized = Lr / teiler;
	(*Bins())[binnumber] += Lr_sized;


	return Lr_sized;

}

glm::dvec3 DDwivediSampling::Position() const
{
	return position;
}

glm::dvec3 DDwivediSampling::Direction() const
{
	return direction;
}



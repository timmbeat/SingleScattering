#pragma once
#include "DwivediSampling.h"
#include "glm.hpp"
class DDwivediSampling :
	public DwivediSampling
{




	public:
	DDwivediSampling(glm::dvec3 position, glm::dvec3 direction, double absorption, double scattering, double anisotropy, double diameter, double delr, std::size_t runs, double wz, int times, bool forcescattering);
	~DDwivediSampling();


	void updateDirection();
	void updatePosition(double stepsize);
	double calculateLr() override;
	void out();
	bool boundary(double stepsize);
	double run();
	void reset();
	double sampleDirDistribution_hen() const;

	glm::dvec3 Position() const;
	glm::dvec3 Direction() const;



private:
glm::dvec3 position;
glm::dvec3 direction;
double wz;
double wz_old;
double Lr;



glm::dvec3 position_org;
glm::dvec3 direction_org;
double wz_org;
double Lr_org;
int times;
bool forcescattering;
int scatterevent;
bool dead;
};


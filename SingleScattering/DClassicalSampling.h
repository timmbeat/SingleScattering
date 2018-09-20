#pragma once
#include "ClassicalSampling.h"
#include <vec3.hpp>

class DClassicalSampling : public ClassicalSampling
{
	public:
	DClassicalSampling(glm::dvec3 position, glm::dvec3 direction, double absorption, double scattering, double anisotropy, double diameter, double delr, std::size_t runs, int times, bool forcescattering);
	~DClassicalSampling();
	double calculateLr() override;
	void reset();
	void updateDirection(double theta);
	void updatePosition(double stepsize);
	void out();
	bool boundary(double stepsize);
	double run();



private:

	
double Lr;
glm::dvec3 position;
glm::dvec3 direction;

glm::dvec3 position_org;
glm::dvec3 direction_org;
int times;
bool forcescattering;
int scatterevent;
bool dead;
};


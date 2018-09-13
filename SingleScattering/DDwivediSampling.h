#pragma once
#include "DwivediSampling.h"
#include "glm.hpp"
class DDwivediSampling :
	public DwivediSampling
{




	public:
	DDwivediSampling(glm::dvec3 position, glm::dvec3 direction, double absorption, double scattering, double anisotropy, double diameter, double delr, std::size_t runs);
	~DDwivediSampling();


	void updateDirection(double wz);
	void updatePosition(double stepsize);
	double calculateLr() override;

	glm::dvec3 Position() const;
	glm::dvec3 Direction() const;



private:
glm::dvec3 position;
glm::dvec3 direction;

};


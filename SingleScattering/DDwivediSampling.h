#pragma once
#include "DwivediSampling.h"
#include "glm.hpp"
class DDwivediSampling :
	public DwivediSampling
{




	public:
	DDwivediSampling(glm::dvec3 position, glm::dvec3 direction, Real absorption, Real scattering, Real anisotropy, double diameter, double delr, std::size_t runs);
	~DDwivediSampling();


	void updateDirection(Real wz);
	void updatePosition(Real stepsize);


	glm::dvec3 Position() const;
	glm::dvec3 Direction() const;

private:
glm::dvec3 position;
glm::dvec3 direction;

};


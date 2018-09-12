#include "DDwivediSampling.h"





DDwivediSampling::DDwivediSampling(glm::dvec3 position, glm::dvec3 direction, Real absorption, Real scattering,
	Real anisotropy, double diameter, double delr, std::size_t runs) :
	position(position), direction(direction), DwivediSampling{absorption, scattering, anisotropy, diameter, delr, runs}
{
}

DDwivediSampling::~DDwivediSampling()
{
}

void DDwivediSampling::updateDirection(Real wz)
{

	//auto const theta = acos(-wz);
	auto const y = 2 * mcss::pi<double>() * random();
	//auto const sint = sqrt(1 - wz * wz);
	auto const sint = sqrt(1 - wz * wz);
	auto const cosp = cos(y);
	double sinp;
	auto const ux = direction.x;
	auto const uy = direction.y;
	auto const uz = direction.z;
	auto const direction_old = direction;





	if (y < mcss::pi<double>())
	{
		sinp = sqrt(1.0 - cosp * cosp);
	}
	else
	{
		sinp = -sqrt(1.0 - cosp * cosp);
	}


	direction.x = sint * cosp;
	direction.y = sint * sinp;
	direction.z = glm::sign(uz)*-wz;

	auto const direction_new = direction;
	auto const dot = glm::dot(direction_new, direction_old);
	auto const dot_theta = dot / (glm::length(direction_old) * glm::length(direction_new));
	photon->directiont = dot_theta;
	auto const hg = directional_distribution_hg(mat->matproperties->anisotropy, acos(dot_theta)); //Bis hier hin ist es richtig!

	auto dwi = directional_distribution_dwi(mat->matproperties->v0, photon->wz_new);




	auto new_weight = (photon->weight * hg) / dwi;


	photon->weight = new_weight;

}

void DDwivediSampling::updatePosition(Real stepsize)
{
	position.x = position.x + direction.x*stepsize;
	position.y = position.y + direction.y*stepsize;
	position.z = position.z + direction.z*stepsize;
}

glm::dvec3 DDwivediSampling::Position() const
{
	return position;
}

glm::dvec3 DDwivediSampling::Direction() const
{
	return direction;
}



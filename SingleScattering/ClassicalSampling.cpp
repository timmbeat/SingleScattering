#include "ClassicalSampling.h"



ClassicalSampling::ClassicalSampling(double const absorption, double const scattering, double const anisotropy, double const diameter, double const delr, std::size_t const runs)
	: Sampling{ absorption, scattering, anisotropy, diameter, delr, runs }
{
}


ClassicalSampling::~ClassicalSampling() = default;


double ClassicalSampling::samplePathDistribution() const
{
	return -log(random()) / (Absorption() + Scattering());
}

double ClassicalSampling::sampleDirDistribution() const
{
	if (Anisotropy() == 0) return 2 * random() - 1;
	auto const qani = Anisotropy() * Anisotropy();

	return (1 / (2 * Anisotropy())) * (1 + qani - pow(((1 - qani) / (1 - Anisotropy() + 2 * Anisotropy()*random())), 2));
}

double ClassicalSampling::calculateLr()
{

	auto const costheta = sampleDirDistribution();
	auto const theta = acos(costheta);

	auto const di = samplePathDistribution();
	auto const r = -di * tan(theta);


	if (r <= 0)
	{
		return 0.0;
	}

	auto const li = Sampling::li(r, di);
	auto const tdi_r = Sampling::taudi_r(li, Absorption(), Scattering());
	auto const to_di = Sampling::tauo_di(di, Absorption(), Scattering());
	auto const pdfp = Sampling::henyey_greenstein_norm(theta, Anisotropy());
	auto const pdi_r = Sampling::henyey_greenstein(theta, Anisotropy());
	auto const pdft = Sampling::taudi_r(di, Absorption(), Scattering());
	auto const lr = Sampling::lr(tdi_r, pdi_r, to_di, pdfp, pdft);

	auto binnumber = static_cast<size_t>(r / Delr());

	if (binnumber > Binsr() - 1) binnumber = Binsr() - 1;

	(*Bins())[binnumber] += lr / Runs();

	return lr / Runs();
}

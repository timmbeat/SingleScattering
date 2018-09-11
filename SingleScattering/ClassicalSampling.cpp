#include "ClassicalSampling.h"



ClassicalSampling::ClassicalSampling(float absorption, float scattering, float anisotropy, float size, float delr)
	: absorption(absorption), scattering(scattering), anisotropy(anisotropy), binsr(size/delr), delr(delr)
{
	bins.resize(binsr);
}


ClassicalSampling::~ClassicalSampling() = default;

double ClassicalSampling::calculatelr(size_t runs)
{
	
		double const costheta = sdirectionaldistribution();
		auto const theta = acos(costheta);
		if (theta <= mcss::pi<double>() / 2)
		{
			return 0.0;
		}
		double const di = spathdistribution();
		auto const r = -di * tan(theta);
		
		double const li = Sampling::li(r, di);

		double const tdi_r = Sampling::taudi_r(li, absorption, scattering);
		double const to_di = Sampling::tauo_di(di, absorption, scattering);
		double const pdfp = Sampling::henyey_greenstein_norm(theta, anisotropy);
		double const pdi_r = Sampling::henyey_greenstein(theta, anisotropy);
		double const pdft = Sampling::taudi_r(di, absorption, scattering);

		double const lr = Sampling::lr(tdi_r, pdi_r, to_di, pdfp, pdft);

		auto binnumber = static_cast<size_t>(r / delr);

		if (binnumber > binsr - 1) binnumber = binsr - 1;
		
		bins[binnumber] += lr/runs;

		return lr/runs;

}

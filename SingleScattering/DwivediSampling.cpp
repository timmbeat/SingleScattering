#include "DwivediSampling.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "ClassicalSampling.h"


double td_ir(double mu_t, double li, double wz, double v0);
double to_di(double mu_t, double di, double wz, double v0, double absorption);
DwivediSampling::DwivediSampling(float const absorption, float const scattering, float const anisotropy, int const binsr, float const delr) :
	absorption(absorption), scattering(scattering), anisotropy(anisotropy), binsr(binsr), delr(delr), v0(mcss::v0(scattering / (absorption + scattering)))
{
	bins.resize(binsr);
}

double DwivediSampling::calculateonelr(size_t const runs)
{

	auto const wz = sampledirdis(v0);
	auto const theta = acos(-wz);

	if (theta >= mcss::pi<double>() / 2)
	{
		return 0.0;
	}

		
		auto const di = samplepathdis(-1.0, v0, absorption + scattering);
		auto const r = -di * tan(theta);
		double const li = Sampling::li(r, di);

		//double const tdi_r = Sampling::taudi_r(li, absorption, scattering);
		double const tdi_r = Sampling::taudi_r(li, absorption, scattering);

		//double const to_di = Sampling::tauo_di(di, absorption, scattering);
		double const to_di = Sampling::tauo_di(di, absorption, scattering);
		double const pdi_r = Sampling::henyey_greenstein(theta, anisotropy);

		auto const pdfp = dirdis(v0, wz);
		auto const pdft = pathdis(v0, wz, absorption + scattering, di);   


		auto binnumber = static_cast<size_t>(di*tan(acos(-wz)) / delr);

		if (binnumber > binsr - 1) binnumber = binsr - 1;
		double const lr = Sampling::lr(tdi_r, pdi_r, to_di, pdfp, pdft);
		bins[binnumber] += lr/runs;

		return lr;




}


double td_ir(double mu_t, double li, double wz, double v0)
{
	return (1 - wz / v0)*mu_t*exp(-(1 - wz / v0)*mu_t*li);
}

double to_di(double mu_t, double di, double wz, double v0, double absorption)
{
	return (1 - wz / v0)*mu_t*exp(-(1 - wz / v0)*mu_t*di)*(1 - absorption * di);
}
int main()
{


	std::ofstream ccout("./plot/Manyplots.csv", std::ofstream::trunc);

	std::stringstream csvout;
	csvout << std::setw(15) << std::left << "DWIVEDI" << std::setw(15) << std::left << "CLASSICAL " << " RUN";
	ccout << csvout.str() << std::endl;
	csvout.str("");


	for (auto i = 0; i < 100; i++)
	{
		DwivediSampling dwi(1.0, 1.0, 0.99, 100, 0.005);
		ClassicalSampling cla(1.0, 1.0, 0.99, 100, 0.005);
		auto sumdwi = 0.0;
		auto sumclas = 0.0;
		auto runs = 100000;

		for (auto j = 0; j < runs; j++)
		{

			sumdwi += dwi.calculateonelr(runs);
			sumclas += cla.calculatelr(runs);
		}
		csvout << std::setw(15) << std::left << sumdwi/runs << std::setw(15) << std::left << sumclas/runs << i;
		ccout << csvout.str() << std::endl;
		csvout.str("");
		Sampling::createPlotFile(dwi.getBins(), cla.getbin(), "./plot/plotfile.csv");

	}


	
}
DwivediSampling::~DwivediSampling() = default;

#include "DwivediSampling.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "ClassicalSampling.h"


DwivediSampling::DwivediSampling(double const absorption, double const scattering, double const anisotropy, float size, float const delr) :
	absorption(absorption), scattering(scattering), anisotropy(anisotropy), binsr(size/delr), delr(delr), v0(mcss::v0(scattering / (absorption + scattering)))
{
	bins.resize(binsr);
}

double DwivediSampling::calculateonelr(size_t const runs)
{

	auto const wz = sampledirdis(v0);
	auto const theta = mcss::pi<double>() - acos(-wz);

	if (theta <= mcss::pi<double>() / 2)
	{
		return 0.0;
	}

		auto const di = samplepathdis(-1.0, v0, absorption + scattering);
		auto const r = -di * tan(theta);
		auto const li = Sampling::li(r, di);

		//double const tdi_r = Sampling::taudi_r(li, absorption, scattering);
		auto const tdi_r = Sampling::taudi_r(li, absorption, scattering);

		//double const to_di = Sampling::tauo_di(di, absorption, scattering);
		auto const to_di = Sampling::tauo_di(di, absorption, scattering);
		auto const pdi_r = Sampling::henyey_greenstein(theta, anisotropy);

		auto const pdfp = dirdis(v0, wz);
		auto const pdft = pathdis(v0, -1.0, absorption + scattering, di);   


		auto binnumber = static_cast<size_t>(r / delr);

		if (binnumber > binsr - 1) binnumber = binsr - 1;
		auto const lr = Sampling::lr(tdi_r, pdi_r, to_di, pdfp, pdft);
		bins[binnumber] += lr/runs;

		return lr/runs;




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
		DwivediSampling dwi(1.0, 1.0, 0.99, 1.0, 0.005);
		ClassicalSampling cla(1.0, 1.0, 0.99, 1.0, 0.005);
		auto sumdwi = 0.0;
		auto sumclas = 0.0;
		auto runs = 1000000;

		for (auto j = 0; j < runs; j++)
		{

			sumdwi += dwi.calculateonelr(runs);
			sumclas += cla.calculatelr(runs);
		}

		std::cout << sumdwi << std::endl;
		std::cout << sumclas << std::endl;
		std::cout << i << std::endl;
		csvout << std::setw(15) << std::left << sumdwi/runs << std::setw(15) << std::left << sumclas/runs << i;
		ccout << csvout.str() << std::endl;
		csvout.str("");
		
		Sampling::createPlotFile(dwi.getBins(), cla.getbin(), dwi.getDelr(),  "./plot/plotfile.csv");

	}


	
}
DwivediSampling::~DwivediSampling() = default;

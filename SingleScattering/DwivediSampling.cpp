#include "DwivediSampling.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "ClassicalSampling.h"
#include "DDwivediSampling.h"


DwivediSampling::DwivediSampling(double absorption, double scattering, double anisotropy, double diameter, double delr, std::size_t runs) :
	Sampling{absorption, scattering, anisotropy, diameter, delr, runs}, v0(mcss::v0(scattering/(absorption+scattering)))
{
}
DwivediSampling::~DwivediSampling() = default;



double DwivediSampling::dirdis(double const wz) const
{
		return (1 / log((v0 + 1) / (v0 - 1))) * (1 / (v0 - wz));
}

double DwivediSampling::pathdis(double const wz, double const t) const
{
		return (1 - wz / v0)*MU_T()*exp(-(1 - wz / v0)*MU_T()*t);	
}

double DwivediSampling::V0() const
{
	return v0;
}

double DwivediSampling::samplePathDistribution(double wz) const
{
	return static_cast<double>(-log(random()) / ((1 - v0 / wz)*MU_T()));
}

double DwivediSampling::sampleDirDistribution() const
{
	return v0 - (v0 + 1)*pow((v0 - 1) / (v0 + 1), random());
}





double DwivediSampling::calculateLr()
{

	auto const wz = sampleDirDistribution();
	auto const theta = acos(-wz);

	auto const di = samplePathDistribution(-1.0);
	auto const r = -di * tan(theta);

	if (r <= 0)
	{
		return 0.0;
	}



	auto const li = Sampling::li(r, di);

	//double const tdi_r = Sampling::taudi_r(li, absorption, scattering);
	auto const tdi_r = Sampling::taudi_r(li, Absorption(), Scattering());

	//double const to_di = Sampling::tauo_di(di, absorption, scattering);
	auto const to_di = Sampling::tauo_di(di, Absorption(), Scattering());
	auto const pdi_r = Sampling::henyey_greenstein(theta, Anisotropy());

	auto const pdfp = dirdis(wz);
	auto const pdft = pathdis(-1.0, di);


	auto binnumber = static_cast<size_t>(r / Delr());

	if (binnumber > Binsr() - 1) binnumber = Binsr() - 1;
	auto const lr = Sampling::lr(tdi_r, pdi_r, to_di, pdfp, pdft);
	(*Bins())[binnumber] += lr / Runs();

	return lr / Runs();

}


int main()
{

	std::size_t runs = 100000;
	double const absorption = 1.0;
	double const scattering = 1.0;
	double const anisotropy = 0.99;

	std::ofstream ccout("./plot/Manyplots.csv", std::ofstream::trunc);

	std::stringstream csvout;
	csvout << std::setw(15) << std::left << "DWIVEDI" << std::setw(15) << std::left << "CLASSICAL " << " RUN";
	ccout << csvout.str() << std::endl;
	csvout.str("");
	for (auto j = 0; j < 100; j++)
	{

		std::cout << "RUN " << j << std::endl;
		DDwivediSampling ddwivedi{ glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.0, 0.0, -1.0) , absorption, scattering, anisotropy, 1.0, 0.005, runs };


		auto sum = 0.0;
		for (auto i = 0; i < runs; i++)
		{
			sum += ddwivedi.calculateLr();
		}

		std::cout << sum;
		Sampling::createPlotFile(ddwivedi.Bins(), ddwivedi.Bins(), 0.005, "./plot/ddwivedi.csv");

		csvout << std::setw(15) << std::left << sum << std::setw(15) << std::left << sum << j;
		ccout << csvout.str() << std::endl;
		csvout.str("");
	}

}

	//int main()
//{
//
//	std::ofstream ccout("./plot/Manyplots.csv", std::ofstream::trunc);
//
//	std::stringstream csvout;
//	csvout << std::setw(15) << std::left << "DWIVEDI" << std::setw(15) << std::left << "CLASSICAL " << " RUN";
//	ccout << csvout.str() << std::endl;
//	csvout.str("");
//
//	std::size_t runs = 100000;
//	Real const absorption = 1.0;
//	Real const scattering = 1.0;
//	Real const anisotropy = 0.99;
//
//	for (auto i = 0; i < 1000; i++)
//	{
//		DwivediSampling dwivedi{ absorption, scattering, anisotropy, 1.0, 0.005, runs };
//		ClassicalSampling classical{ absorption, scattering, anisotropy, 1.0, 0.005, runs };
//		auto sumdwi = 0.0;
//		auto sumclas = 0.0;
//		
//		for (auto j = 0; j < runs; j++)
//		{
//
//			sumdwi += dwivedi.calculateLr();
//			sumclas += classical.calculateLr();
//		}
//
//		std::cout << sumdwi << std::endl;
//		std::cout << sumclas << std::endl;
//		std::cout << i << std::endl;
//		csvout << std::setw(15) << std::left << sumdwi/runs << std::setw(15) << std::left << sumclas/runs << i;
//		ccout << csvout.str() << std::endl;
//		csvout.str("");
//		
//		Sampling::createPlotFile(dwivedi.Bins(), classical.Bins(), dwivedi.Delr(),  "./plot/plotfile.csv");
//
//	}
//
//
//	
//}

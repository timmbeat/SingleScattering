#include "DwivediSampling.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "ClassicalSampling.h"


DwivediSampling::DwivediSampling(float const absorption, float const scattering, float const anisotropy, int const binsr, float const delr) :
	absorption(absorption), scattering(scattering), anisotropy(anisotropy), binsr(binsr), delr(delr), v0(mcss::v0(scattering / (absorption + scattering)))
{
	bins.resize(binsr);
}

double DwivediSampling::calculateonelr()
{

	auto const wz = sampledirdis(v0);
	auto const theta = acos(-wz);

	if (theta >= mcss::pi<double>() / 2)
	{
		return 0.0;
	}

		
		auto const di = samplepathdis(wz, v0, absorption + scattering);
		auto const r = -di * tan(theta);
		double const li = Sampling::li(r, di);

		double const tdi_r = Sampling::taudi_r(li, absorption, scattering);
		double const to_di = Sampling::tauo_di(di, absorption, scattering);

		double const pdi_r = Sampling::henvey_greenstein(theta, anisotropy);

		auto const pdfp = dirdis(v0, wz);
		auto const pdft = pathdis(v0, wz, absorption + scattering, di);   


		auto binnumber = static_cast<size_t>(di/tan(acos(-wz)) / delr);

		if (binnumber > binsr - 1) binnumber = binsr - 1;
		double const lr = Sampling::lr(tdi_r, pdi_r, to_di, pdfp, pdft);
		bins[binnumber] += lr;

		return lr;




}




int main()
{
	DwivediSampling dwi(1.0, 1.0, 0.99, 100, 0.005) ;
	ClassicalSampling cla(1.0, 1.0, 0.99, 100, 0.005);
	auto sumdwi = 0.0;
	auto sumclas = 0.0;
	auto runs = 1000;
	for(auto i = 0; i < runs; i++)
	{

		sumdwi += dwi.calculateonelr();
		sumclas += cla.calculatelr();
	}
	

	std::cout << sumdwi / runs << "DWIVEDI";
	std::cout << sumclas / runs << "CLASSICAL";
	std::cin.get();
	Sampling::createPlotFile(dwi.getBins(), cla.getbin(), "./plot/plotfile.csv");

	
}
DwivediSampling::~DwivediSampling() = default;

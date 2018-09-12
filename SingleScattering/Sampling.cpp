#include "Sampling.h"
#include <fstream>
#include <sstream>
#include <iomanip>


Sampling::Sampling(Real absorption, Real scattering, Real anisotropy, double diameter, double delr, std::size_t runs) :
	absorption(absorption), scattering(scattering), anisotropy(anisotropy), diameter(diameter), delr(delr), binsr(diameter / delr), runs(runs), mu_t(absorption + scattering)
{
	bins.resize(binsr);
}


Sampling::~Sampling()
{
}



/*
 * SAMPLING FUNCTIONS FOR BOTH SAMPLING TECHNIQUES
 */
Real Sampling::tauo_di(Real const di, Real const absorption, Real const scattering)
{
	return (absorption + scattering)*exp(-(absorption + scattering)*di)*(1 - absorption * di);
}

Real Sampling::taudi_r(Real const li, Real const absorption, Real const scattering)
{
	return (absorption + scattering)*exp(-(absorption + scattering)*li);
}

Real Sampling::henyey_greenstein(Real const theta, Real const g)
{
	return (1 - g * g) / (4 * pi<double>()*pow(1 + g * g - 2 * g*cos(theta), 3.0 / 2.0));
}

Real Sampling::henyey_greenstein_norm(Real const theta, Real const g)
{
	return (1 - g * g) / (2 * pow(1 + g * g - 2 * g*cos(theta), 3.0 / 2.0));
}

Real Sampling::lr(Real const taudi_r, Real const pdi_r, Real const tauo_di, Real const pdfp, Real const pdftau)
{
	return taudi_r * (pdi_r / pdfp)*(tauo_di / pdftau);
}
Real Sampling::li(Real r, Real di)
{

	return sqrt(di*di + r * r);
}


/*
 * Creating the File for the TeX - Plot
 */
void Sampling::createPlotFile(const std::vector<Real> *  binsA, const std::vector<Real> * binsB, const float delr, std::string filename)
{



	std::ofstream ccout(filename, std::ofstream::trunc);

	std::stringstream csvout;
	csvout << std::setw(15) << std::left << "weight_A" << std::setw(15) << std::left << "weight_B " << " ir";
	ccout << csvout.str() << std::endl;

	csvout.str("");

	for (size_t i = 0; i < binsA->size(); i++)
	{
		csvout << std::setw(15) << std::left << (*binsA)[i] << std::setw(15) << std::left << (*binsB)[i] << " " << i * delr;
		ccout << csvout.str() << std::endl;
		csvout.str("");
	}


}




/*
 * GETTER
 */
float Sampling::Diameter() const
{
	return diameter;
}

std::size_t Sampling::Runs() const
{
	return runs;
}

Real Sampling::MU_T() const
{
	return mu_t;
}

std::vector<Real> * Sampling::Bins()
{
	return &bins;
}

float Sampling::Delr() const
{
	return delr;
}

Real Sampling::Anisotropy() const
{
	return anisotropy;
}

Real Sampling::Scattering() const
{
	return scattering;
}

Real Sampling::Absorption() const
{
	return absorption;
}

std::size_t Sampling::Binsr() const
{
	return binsr;
}

#include "Sampling.h"
#include <fstream>
#include <sstream>
#include <iomanip>


Sampling::Sampling(double absorption, double scattering, double anisotropy, double diameter, double delr, std::size_t runs) :
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
double Sampling::tauo_di(double const di, double const absorption, double const scattering)
{
	return (absorption + scattering)*exp(-(absorption + scattering)*di)*(1 - absorption * di);
}

double Sampling::taudi_r(double const li, double const absorption, double const scattering)
{
	return (absorption + scattering)*exp(-(absorption + scattering)*li);
}

double Sampling::henyey_greenstein(double const theta, double const g)
{
	return (1 - g * g) / (4 * pi<double>()*pow(1 + g * g - 2 * g*cos(theta), 3.0 / 2.0));
}

double Sampling::henyey_greenstein_norm(double const cos_theta, double const g)
{
	return (1 - g * g) / (2 * pow(1 + g * g - 2 * g* cos_theta, 3.0 / 2.0));
}

double Sampling::lr(double const taudi_r, double const pdi_r, double const tauo_di, double const pdfp, double const pdftau)
{
	return taudi_r * (pdi_r / pdfp)*(tauo_di / pdftau);
}
double Sampling::li(double r, double di)
{

	return sqrt(di*di + r * r);
}


/*
 * Creating the File for the TeX - Plot
 */
void Sampling::createPlotFile(const std::vector<double> *  binsA, const std::vector<double> * binsB, const float delr, std::string filename)
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

double Sampling::MU_T() const
{
	return mu_t;
}

std::vector<double> * Sampling::Bins()
{
	return &bins;
}

float Sampling::Delr() const
{
	return delr;
}

double Sampling::Anisotropy() const
{
	return anisotropy;
}

double Sampling::Scattering() const
{
	return scattering;
}

double Sampling::Absorption() const
{
	return absorption;
}

std::size_t Sampling::Binsr() const
{
	return binsr;
}

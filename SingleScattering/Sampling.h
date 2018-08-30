#pragma once
#include <cmath>
#include <random>
#include "constants.h"

static std::random_device rd;
static std::mt19937 gen(rd());
static std::uniform_real_distribution<double> dis(0.0, 1.0);

using namespace mcss;
class Sampling
{
	public:
	Sampling();
	~Sampling();




	static double tauo_di(double const di, double const absorption, double const scattering) 
	{
		return (absorption + scattering)*exp(-(absorption + scattering)*di)*(1 - absorption * di);
	}

	static double taudi_r(double const li, double const absorption, double const scattering)
	{
		return (absorption + scattering)*exp(-(absorption + scattering)*li);
	}
	static double henvey_greenstein(double const theta, double const g)
	{
		return (1 - g * g) / (4 * pi<double>()*pow(1 + g * g - 2 * g*cos(theta), 3.0 / 2.0));
	}

	static double henvey_greenstein_norm(double const theta, double const g)
	{
		return (1 - g * g) / (2 * pow(1 + g * g - 2 * g*cos(theta), 3.0 / 2.0));
	}

	static double lr(double const taudi_r, double const pdi_r, double const tauo_di, double const pdfp, double const pdftau)
	{
		return taudi_r * (pdi_r / pdfp)*(tauo_di / pdftau);
	}

	static double li(double r, double di)
	{
		return sqrt(di*di + r * r);
	}


	static void createPlotFile(const std::vector<double> *  binsA, const std::vector<double> * binsB, std::string filename);


	
};

static double random()
{
	return dis(gen);
}
/*
* Created by Tim Mend on 12.09.2018
*/


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
	Sampling(double absorption, double scattering, double  anisotropy, double diameter, double delr, std::size_t runs);
	virtual ~Sampling();



public:
//Function which are needed by both Sampling Method
double tauo_di(double const di, double const absorption, double const scattering);
double taudi_r(double const li, double const absorption, double const scattering);
double henyey_greenstein(double const theta, double const g);
double henyey_greenstein_norm(double const cos_theta, double const g);
double lr(double const taudi_r, double const pdi_r, double const tauo_di, double const pdfp, double const pdftau);
double li(double r, double di);
static void createPlotFile(const std::vector<double> *  binsA, const std::vector<double> * binsB, const float delr,  std::string filename);
virtual double calculateLr() = 0;


//Getter
std::size_t Binsr()	const;
double Absorption()	const;
double Scattering()	const;
double Anisotropy()	const;
float Delr()		const;
std::vector<double> *  Bins();
float Diameter()	const;
std::size_t Runs()	const;
double MU_T() const;

private:
double absorption;
double scattering;
double anisotropy;
double mu_t;
std::size_t binsr;
std::size_t runs;
float diameter;
float delr;
std::vector<double> bins;
	
};

static double random()
{
	return dis(gen);
}
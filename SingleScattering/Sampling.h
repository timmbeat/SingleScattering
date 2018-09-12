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
	Sampling(Real absorption, Real scattering, Real anisotropy, double diameter, double delr, std::size_t runs);
	virtual ~Sampling();



public:
//Function which are needed by both Sampling Method
Real tauo_di(Real const di, Real const absorption, Real const scattering);
Real taudi_r(Real const li, Real const absorption, Real const scattering);
Real henyey_greenstein(Real const theta, Real const g);
Real henyey_greenstein_norm(Real const theta, Real const g);
Real lr(Real const taudi_r, Real const pdi_r, Real const tauo_di, Real const pdfp, Real const pdftau);
Real li(Real r, Real di);
static void createPlotFile(const std::vector<Real> *  binsA, const std::vector<Real> * binsB, const float delr,  std::string filename);
virtual Real calculateLr() = 0;


//Getter
std::size_t Binsr()	const;
Real Absorption()	const;
Real Scattering()	const;
Real Anisotropy()	const;
float Delr()		const;
std::vector<Real> *  Bins();
float Diameter()	const;
std::size_t Runs()	const;
Real MU_T() const;

private:
Real absorption;
Real scattering;
Real anisotropy;
Real mu_t;
std::size_t binsr;
std::size_t runs;
float diameter;
float delr;
std::vector<Real> bins;
	
};

static double random()
{
	return dis(gen);
}
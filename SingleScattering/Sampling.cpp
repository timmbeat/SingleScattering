#include "Sampling.h"
#include <fstream>
#include <sstream>
#include <iomanip>


Sampling::Sampling()
{
}


Sampling::~Sampling()
{
}

void Sampling::createPlotFile(const std::vector<double>  * binsA, const std::vector<double> * binsB, std::string filename)
{

	
	
		std::ofstream ccout(filename, std::ofstream::trunc);

		std::stringstream csvout;
		csvout << std::setw(15) << std::left << "weight_A" << std::setw(15) << std::left << "weight_B "<< " ir";
		ccout << csvout.str() << std::endl;

		csvout.str("");

		for (size_t i = 0; i < binsA->size(); i++)
		{
			csvout << std::setw(15) << std::left << (*binsA)[i] << std::setw(15) << std::left << (*binsB)[i] << " " << i;
			ccout << csvout.str() << std::endl;
			csvout.str("");
		}
	

}

#include "experiment.h"

Experiment::Experiment(
		std::string output_file,
		const std::function<double(unsigned long)>& pH,
		const std::function<double(unsigned long, double)>& dP
		)
	: system(ChemicalSystem::defaultEquilibrium()), csv_file(output_file), pH(pH), dP(dP) {
	csv_file << "T,H,P,C,S,SH,SP,SC" << std::endl;
}

void Experiment::writeToCsv() {
	csv_file << t << ","
		<< system.cH() << ","
		<< system.cP() << ","
		<< system.cC() << ","
		<< system.cS() << ","
		<< system.cSH() << ","
		<< system.cSP() << ","
		<< system.cSC() << std::endl;

	//csv_file << t << ","
		//<< system.nH() << ","
		//<< system.nP() << ","
		//<< system.nC() << ","
		//<< system.nS() << ","
		//<< system.nSH() << ","
		//<< system.nSP() << ","
		//<< system.nSC() << std::endl;
}

void Experiment::run(unsigned long time) {
	for(t = 0; t < time; t++) {
		std::cout << "Exp pH: " << pH(t) << std::endl;
		system.incrementP(dP(time, system.nP()));
		system.setPH(pH(t));
		system.distanceToEquilibrium();
		writeToCsv();
	}
}

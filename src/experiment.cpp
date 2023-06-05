#include "experiment.h"

namespace mineral {
	Experiment::Experiment(
			std::string output_file,
			const ChemicalSystem& chemical_system,
			const std::function<double(unsigned long)>& pH,
			const std::function<double(unsigned long, double)>& dP
			)
		: system(chemical_system), csv_file(output_file), pH(pH), dP(dP) {
			csv_file << "T,H,P,C,S,SH,SP,SC,pH" << std::endl;
		}

	void Experiment::writeToCsv() {
		//csv_file << t << ","
			//<< system.cH() << ","
			//<< system.cP() << ","
			//<< system.cC() << ","
			//<< system.cS() << ","
			//<< system.cSH() << ","
			//<< system.cSP() << ","
			//<< system.cSC() << std::endl;

		csv_file << t << ","
		<< system.nH() << ","
		<< system.nP() << ","
		<< system.nC() << ","
		<< system.nS() << ","
		<< system.nSH() << ","
		<< system.nSP() << ","
		<< system.nSC() << ","
		<< system.pH() << std::endl;
	}

	void Experiment::run(unsigned long time) {
		t = 0;
		while(t < time) {
			std::cout << "Exp pH: " << pH(t) << std::endl;
			writeToCsv();
			system.incrementP(dP(time, system.nP()));
			//writeToCsv();
			//++t;
			system.solveEquilibrium();
			//system.setPH(pH(t));
			++t;
		}
	}
}

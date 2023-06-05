#include "experiment.h"

using namespace mineral;

int main(int, char *[])
{
	{
		auto system = ChemicalSystem::Devau2011Control();
		double end_pH1 = 2;
		double end_pH2 = 11;
		unsigned long n = 1000;
		Experiment ph_experiment(
				"ph_output.csv", system,
				[&system, end_pH1, end_pH2, n] (unsigned long time) {
				if(time < n)
					return system.initPH() + ((double) time)*(end_pH1 - system.initPH())/((double) n);
				else
					return end_pH1 + ((double) time-n)*(end_pH2 - end_pH1)/((double) 2*n);

				},
				[] (unsigned long, double) {return 0;}
				);
		ph_experiment.run(3*n);
	}
	{
		//double init_P = 5.4e-8;
		//double min_P = 1e-8;
		unsigned long n = 1000;
		auto system = ChemicalSystem::Devau2011Control();
		auto init_P = system.nP();
		Experiment P_experiment(
				"P_output.csv", system,
				[] (unsigned long) {
				return 7.5;
				},
				//[init_P, min_P, n] (unsigned long, double) {
				//return (min_P-init_P)/(double) n;
				//}
				[init_P] (unsigned long, double) {
				return +0.5*init_P;
				}
				);
		P_experiment.run(n);
	}
	
	return 0;
}

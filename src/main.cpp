#include "experiment.h"

using namespace mineral;

int main(int, char *[])
{
	/*
	 *ChemicalSystem system;
	 *system.distanceToEquilibrium();
	 *system.setPH(7.6);
	 *system.distanceToEquilibrium();
	 *system.setPH(8);
	 *system.distanceToEquilibrium();
	 */

	{
		double init_pH = 7.5;
		double end_pH = 2;
		unsigned long n = 100;
		Experiment ph_experiment(
				"ph_output.csv",
				[init_pH, end_pH, n] (unsigned long time) {
				return init_pH + ((double) time)*(end_pH - init_pH)/((double) n);
				},
				[] (unsigned long, double) {return 0;}
				);
		ph_experiment.run(n);
	}
	/*
	 *{
	 *    //double init_P = 5.4e-8;
	 *    //double min_P = 1e-8;
	 *    unsigned long n = 1000;
	 *    Experiment P_experiment(
	 *            "P_output.csv",
	 *            [] (unsigned long) {
	 *            return 7.5;
	 *            },
	 *            //[init_P, min_P, n] (unsigned long, double) {
	 *            //return (min_P-init_P)/(double) n;
	 *            //}
	 *            [] (unsigned long, double current_P) {
	 *            return -current_P/10.0;
	 *            }
	 *            );
	 *    P_experiment.run(n);
	 *}
	 */
	
	return 0;
}

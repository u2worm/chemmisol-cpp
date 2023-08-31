#include "chemmisol.h"

using namespace chemmisol;

int main(int, char *[])
{
	ChemicalSystem chemical_system;
	// Defines the reaction H2O <-> OH- + H+ (log K = -13.997)
	chemical_system.addReaction("OH-", -13.997, {
			{"OH-", -1},
			{"H+", -1},
			{"H2O", 1}
			});
	// Defines the reaction Na+ + Cl- <-> NaCl (log K = -0.3)
	chemical_system.addReaction("NaCl", -0.3, {
			{"NaCl", -1},
			{"Na+", AQUEOUS, 1}, // AQUEOUS specified here only for
								 // demonstration purpose, since it should be
								 // the default phase
			{"Cl-", AQUEOUS, 1}
			});
	// Defines the reaction H2O + Na+ <-> NaOH + H+ (log K = -13.897)
	chemical_system.addReaction("NaOH", -13.897, {
			{"NaOH", -1},
			{"H+", -1},
			{"Na+", 1},
			{"H2O", 1}
			});

	// Defines the Na+ component and sets its total concentration to 0.1 mol/l
	chemical_system.addComponent("Na+", 0.1*mol/l);
	// Defines the Cl- component and sets its total concentration to 0.1 mol/l
	chemical_system.addComponent("Cl-", 0.1*mol/l);
	// Defines the H2O component as a solvent
	chemical_system.addSolvent("H2O");

	// Sets up a basic CSV output
	std::ofstream csv_file("output.csv");
	csv_file << "i, pH, Na+, Cl-, NaCl, NaOH"
		<< std::endl;

	for(std::size_t i = 0; i <= 10; i++) {
		// Makes the ph vary from 3 to 8,
		chemical_system.fixPH(3 + i*0.5);
		// Reduces the total quantity of Na+ from 0.1 mol/l to 0.05 mol/l
		chemical_system.setTotalConcentration(
				chemical_system.getComponent("Na+"),
				0.1*mol/l - i*0.005*mol/l);

		// Solves the new equilibrium state, from the previous equilibrium
		chemical_system.solveEquilibrium();

		// Concentration output (notice getSpecies() is used, not getComponent())
		csv_file << i << "," <<
			chemical_system.getPH() << "," <<
			chemical_system.getSpecies("Na+").concentration() << "," <<
			chemical_system.getSpecies("Cl-").concentration() << "," <<
			chemical_system.getSpecies("NaCl").concentration() << "," <<
			chemical_system.getSpecies("NaOH").concentration() << std::endl;
	} 
}

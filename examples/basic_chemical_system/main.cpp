#include "chemmisol.h"

using namespace chemmisol;

int main(int, char *[])
{
	ChemicalSystem chemical_system;
	chemical_system.addReaction("OH-", -13.997, {
			{"OH-", -1},
			{"H+", -1},
			{"H2O", 1}
			});
	chemical_system.addReaction("NaCl", -0.3, {
			{"NaCl", -1},
			{"Na+", AQUEOUS, 1}, // AQUEOUS specified here only for
								 // demonstration purpose, since it should be
								 // the default phase
			{"Cl-", AQUEOUS, 1}
			});
	chemical_system.addReaction("NaOH", -13.897, {
			{"NaOH", -1},
			{"H+", -1},
			{"Na+", 1},
			{"H2O", 1}
			});

	chemical_system.addComponent("Na+", 0.1*mol/l);
	chemical_system.addComponent("Cl-", 0.1*mol/l);
	chemical_system.addSolvent("H2O");

	// Automatically adds the H+ component
	chemical_system.fixPH(7);

	chemical_system.solveEquilibrium();
}

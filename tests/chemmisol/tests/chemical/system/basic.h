#include "gmock/gmock.h"
#include "chemmisol/chemical/system.h"

using namespace testing;
using namespace chemmisol;

/**
 * Tests basic ChemicalSystem features for a simple reaction system with only
 * unitary stoichiometric coefficients. Equilibrium solving features are not
 * tested yet.
 */
class BasicAqueousChemicalSystemTest : public Test {
	protected:
		ChemicalSystem chemical_system;
		const Reaction* oh_reaction;
		const Reaction* nacl_reaction;
		const Reaction* naoh_reaction;

		void SetUp() override {
			oh_reaction = &chemical_system.addReaction("OH-", -13.997, {
					{"OH-", -1},
					{"H+", -1},
					{"H2O", 1}
					});
			nacl_reaction = &chemical_system.addReaction("NaCl", -0.3, {
					{"NaCl", -1},
					{"Na+", AQUEOUS, 1}, // AQUEOUS specified here only for test
										 // purpose, since it should be the default
										 // phase
					{"Cl-", AQUEOUS, 1}
					});
			naoh_reaction = &chemical_system.addReaction("NaOH", -13.897, {
					{"NaOH", -1},
					{"H+", -1},
					{"Na+", 1},
					{"H2O", 1}
					});


			chemical_system.addComponent("Na+", 0.1*mol/l);
			chemical_system.addComponent("Cl-", 0.1*mol/l);
			chemical_system.addSolvent("H2O");

			// Automatically adds the H+ component
			chemical_system.initPH(7);

			chemical_system.setUp();
		}

};

class BasicMineralChemicalSystemTest : public Test {
	protected:
		ChemicalSystem chemical_system {
			2.5 * g/l,
			24.2 * m2/g,
			0.8 * entities/nm2
		};

		void SetUp() override {
			chemical_system.addReaction("OH-", -13.997, {
					{"OH-", -1},
					{"H+", -1},
					{"H2O", 1}
					});
			chemical_system.addReaction("=SOH2", 3.46, {
					{"=SOH2", MINERAL, -1},
					{"=SOH", MINERAL, 1}, // AQUEOUS specified here only for test
										  // purpose, since it should be the default
										  // phase
					{"H+", AQUEOUS, 1}
					});
			chemical_system.addComponent("=SOH", MINERAL, 1.0);
			chemical_system.addSolvent("H2O");

			chemical_system.fixPH(7);

			chemical_system.setUp();
		}
};

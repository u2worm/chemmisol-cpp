#include "chemmisol/chemical/system.h"
#include <gmock/gmock.h>

using namespace chemmisol;
using namespace testing;

TEST(InvalidReaction, empty_reagents_in_reaction_set_up) {
	ChemicalSystem system;
	system.addReaction("O2", -10, {});

	EXPECT_THROW(system.setUp(), EmptyReagents);
	try {
		system.setUp();
	} catch(const EmptyReagents& e) {
		ASSERT_THAT(e.getChemicalSystem(), Ref(system));
		ASSERT_THAT(e.getInvalidReaction(), Ref(system.getReaction("O2")));
		CHEM_LOGV(5) << "Thrown EmptyReagentsMissingProducedSpeciesInReaction exception:" << std::endl << e.what();
	}

	EXPECT_THROW(system.solveEquilibrium(), EmptyReagents);
	try {
		system.solveEquilibrium();
	} catch(const EmptyReagents& e) {
		ASSERT_THAT(e.getChemicalSystem(), Ref(system));
		ASSERT_THAT(e.getInvalidReaction(), Ref(system.getReaction("O2")));
		CHEM_LOGV(5) << "Thrown EmptyReagents exception:" << std::endl << e.what();
	}
}

TEST(InvalidReaction, missing_produced_species_in_reaction_set_up) {
	ChemicalSystem system;
	system.addComponent("H2O", AQUEOUS);
	system.addComponent("H2", AQUEOUS);
	system.addComponent("O2", AQUEOUS);
	system.addReaction("O2", -10, {
			{"H2O", 2},
			{"H2", -2},
			{"O2", -1}
			});

	EXPECT_THROW(system.setUp(), MissingProducedSpeciesInReaction);
	try {
		system.setUp();
	} catch(const MissingProducedSpeciesInReaction& e) {
		ASSERT_THAT(e.getChemicalSystem(), Ref(system));
		ASSERT_THAT(e.getInvalidReaction(), Ref(system.getReaction("O2")));
		CHEM_LOGV(5) << "Thrown MissingProducedSpeciesInReaction exception:" << std::endl << e.what();
	}

	EXPECT_THROW(system.solveEquilibrium(), MissingProducedSpeciesInReaction);
	try {
		system.solveEquilibrium();
	} catch(const MissingProducedSpeciesInReaction& e) {
		ASSERT_THAT(e.getChemicalSystem(), Ref(system));
		ASSERT_THAT(e.getInvalidReaction(), Ref(system.getReaction("O2")));
		CHEM_LOGV(5) << "Thrown MissingProducedSpeciesInReaction exception:" << std::endl << e.what();
	}
}

TEST(InvalidReaction, invalid_species_in_reaction) {
	ChemicalSystem system;
	system.addComponent("H2O", AQUEOUS);
	system.addReaction("O2", -10, {
			{"H2O", 2},
			{"H2", -2},
			{"O2", -1}
			});

	EXPECT_THROW(system.setUp(), TooManyProducedSpeciesInReaction);
	try {
		system.solveEquilibrium();
	} catch(const TooManyProducedSpeciesInReaction& e) {
		ASSERT_THAT(e.getChemicalSystem(), Ref(system));
		ASSERT_THAT(e.getInvalidReaction(), Ref(system.getReaction("O2")));
		CHEM_LOGV(5) << "Thrown InvalidReaction exception:" << std::endl << e.what();
	}
	EXPECT_THROW(system.solveEquilibrium(), TooManyProducedSpeciesInReaction);
	try {
		system.solveEquilibrium();
	} catch(const TooManyProducedSpeciesInReaction& e) {
		ASSERT_THAT(e.getChemicalSystem(), Ref(system));
		ASSERT_THAT(e.getInvalidReaction(), Ref(system.getReaction("O2")));
		CHEM_LOGV(5) << "Thrown InvalidReaction exception:" << std::endl << e.what();
	}
}

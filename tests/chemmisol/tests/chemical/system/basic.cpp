#include "chemmisol/tests/chemical/system/basic.h"
#include "chemmisol/tests/chemical/system/utils.h"
#include <random>
#include <regex>

TEST_F(BasicAqueousChemicalSystemTest, get_components) {
	auto& components = chemical_system.getComponents();

	for(std::size_t i = 0; i < components.size(); i++)
		ASSERT_THAT(components[i]->getIndex(), i);
}

TEST_F(BasicAqueousChemicalSystemTest, get_species) {
	auto& species = chemical_system.getSpecies();

	for(std::size_t i = 0; i < species.size(); i++)
		ASSERT_THAT(species[i]->getIndex(), i);
}

TEST_F(BasicAqueousChemicalSystemTest, get_reaction) {
	auto& reactions = chemical_system.getReactions();

	for(std::size_t i = 0; i < reactions.size(); i++)
		ASSERT_THAT(reactions[i]->getIndex(), i);
	ASSERT_THAT(chemical_system.getReaction("OH-"), Ref(*oh_reaction));
	ASSERT_THAT(chemical_system.getReaction("NaCl"), Ref(*nacl_reaction));
	ASSERT_THAT(chemical_system.getReaction("NaOH"), Ref(*naoh_reaction));
}

TEST_F(BasicAqueousChemicalSystemTest, get_component_reagent) {
	ASSERT_THAT(chemical_system.getComponentReagents(chemical_system.getReaction("OH-")),
				UnorderedElementsAre(
					AllOf(
						Field(&ComponentReagent::coefficient, -1),
						Field(
							&ComponentReagent::component,
							&chemical_system.getComponent("H+")
							)
						)
				));
	ASSERT_THAT(chemical_system.getComponentReagents(chemical_system.getReaction("NaCl")),
				UnorderedElementsAre(
					AllOf(
						Field(&ComponentReagent::coefficient, 1),
						Field(
							&ComponentReagent::component,
							&chemical_system.getComponent("Na+")
							)
						),
					AllOf(
						Field(&ComponentReagent::coefficient, 1),
						Field(
							&ComponentReagent::component,
							&chemical_system.getComponent("Cl-")
							)
						)
				));
	ASSERT_THAT(chemical_system.getComponentReagents(chemical_system.getReaction("NaOH")),
				UnorderedElementsAre(
					AllOf(
						Field(&ComponentReagent::coefficient, -1),
						Field(
							&ComponentReagent::component,
							&chemical_system.getComponent("H+")
							)
						),
					AllOf(
						Field(&ComponentReagent::coefficient, 1),
						Field(
							&ComponentReagent::component,
							&chemical_system.getComponent("Na+")
							)
						)
					));

}

TEST_F(BasicAqueousChemicalSystemTest, get_species_reagents) {
	ASSERT_THAT(chemical_system.getSpeciesReagent(chemical_system.getReaction("OH-")),
			AllOf(
				Field(&ChemicalSpeciesReagent::coefficient, -1),
				Field(&ChemicalSpeciesReagent::species,
					&chemical_system.getSpecies("OH-"))
				)
			);

	ASSERT_THAT(chemical_system.getSpeciesReagent(chemical_system.getReaction("NaCl")),
			AllOf(
				Field(&ChemicalSpeciesReagent::coefficient, -1),
				Field(&ChemicalSpeciesReagent::species,
					&chemical_system.getSpecies("NaCl"))
				)
			);

	ASSERT_THAT(chemical_system.getSpeciesReagent(chemical_system.getReaction("NaOH")),
			AllOf(
				Field(&ChemicalSpeciesReagent::coefficient, -1),
				Field(&ChemicalSpeciesReagent::species,
					&chemical_system.getSpecies("NaOH"))
				)
			);

}

TEST_F(BasicAqueousChemicalSystemTest, aqueous_species) {
	ASSERT_THAT(chemical_system.getComponents(), UnorderedElementsAre(
				Pointee(AllOf(
						Property(&ChemicalComponent::getSpecies,
							Pointee(Property(&ChemicalSpecies::getName, "H+"))
						),
						Property(&ChemicalComponent::getTotalQuantity,
							std::pow(10, -7) * AqueousSpecies::V
						)
						)),
				Pointee(AllOf(
						Property(&ChemicalComponent::getSpecies,
							Pointee(Property(&ChemicalSpecies::getName, "Na+"))
						),
						Property(&ChemicalComponent::getTotalQuantity, 0.1*mol/l)
						)),
				Pointee(AllOf(
						Property(&ChemicalComponent::getSpecies,
						Pointee(Property(&ChemicalSpecies::getName, "Cl-"))
						),
						Property(&ChemicalComponent::getTotalQuantity, 0.1*mol/l)
						)),
				Pointee(AllOf(
						Property(&ChemicalComponent::getSpecies,
						Pointee(Property(&ChemicalSpecies::getName, "H2O"))
						)
						))
				));
	// Test fixed components handling
	for(auto& component : chemical_system.getComponents()) {
		if(component->isFixed())
			ASSERT_THAT(
					component.get()->getSpecies(),
					WhenDynamicCastTo<FixedAqueousSpecies*>(Not(IsNull()))
					);
		else
			ASSERT_THAT(
					component.get()->getSpecies(),
					WhenDynamicCastTo<AqueousSpecies*>(Not(IsNull()))
					);
	}

	ASSERT_THAT(chemical_system.getSpecies(), UnorderedElementsAre(
				Pointee(Property(&ChemicalSpecies::getName, "H2O")),
				Pointee(Property(&ChemicalSpecies::getName, "H+")),
				Pointee(Property(&ChemicalSpecies::getName, "OH-")),
				Pointee(Property(&ChemicalSpecies::getName, "Na+")),
				Pointee(Property(&ChemicalSpecies::getName, "Cl-")),
				Pointee(Property(&ChemicalSpecies::getName, "NaCl")),
				Pointee(Property(&ChemicalSpecies::getName, "NaOH"))
				));
	// Test all species
	for(auto& species : chemical_system.getSpecies()) {
		if(species->getName() == "H2O") {
			ASSERT_THAT(species.get(), WhenDynamicCastTo<Solvent*>(Not(IsNull())));
		} else {
			ASSERT_THAT(
					species.get(),
					// An aqueous species might be Fixed or not
					AnyOf(
						WhenDynamicCastTo<AqueousSpecies*>(Not(IsNull())),
						WhenDynamicCastTo<FixedAqueousSpecies*>(Not(IsNull()))
						)
					);
		}
	}
}

TEST_F(BasicAqueousChemicalSystemTest, copy_constructor) {
	ChemicalSystem other_system(chemical_system);

	ASSERT_THAT(
			chemical_system.getComponents(),
			SizeIs(other_system.getComponents().size())
			);
	for(std::size_t i = 0; i < chemical_system.getComponents().size(); i++) {
		const ChemicalComponent& component = *chemical_system.getComponents()[i];
		const ChemicalComponent& other_component = *other_system.getComponents()[i];
		ASSERT_THAT(component,
				AllOf(
					Property(&ChemicalComponent::getIndex, other_component.getIndex()),
					Property(&ChemicalComponent::isFixed, other_component.isFixed()),
					Property(&ChemicalComponent::getSpecies, Pointee(
							Property(
								&ChemicalSpecies::getIndex,
								other_component.getSpecies()->getIndex()
								)
							))
					)
				);
	}

	ASSERT_THAT(
			chemical_system.getSpecies(),
			SizeIs(other_system.getSpecies().size())
			);
	for(std::size_t i = 0; i < chemical_system.getSpecies().size(); i++) {
		const ChemicalSpecies& species = *chemical_system.getSpecies()[i];
		const ChemicalSpecies& other_species = *other_system.getSpecies()[i];
		ASSERT_THAT(species,
				AllOf(
					Property(&ChemicalSpecies::getIndex, other_species.getIndex()),
					Property(&ChemicalSpecies::getName, other_species.getName()),
					Property(&ChemicalSpecies::getPhase, other_species.getPhase()),
					Property(&ChemicalSpecies::concentration, other_species.concentration())
					)
				);
	}

	ASSERT_THAT(
			chemical_system.getReactions(),
			SizeIs(other_system.getReactions().size())
			);
	for(std::size_t i = 0; i < chemical_system.getReactions().size(); i++) {
		const Reaction& reaction = *chemical_system.getReactions()[i];
		const Reaction& other_reaction = *other_system.getReactions()[i];
		ASSERT_THAT(reaction,
				AllOf(
					Property(&Reaction::getIndex, other_reaction.getIndex()),
					Property(&Reaction::getName, other_reaction.getName()),
					Property(&Reaction::getReagents, ElementsAreArray(
							other_reaction.getReagents())
						)
					)
				);
	}

	chemical_system.solveEquilibrium();
	other_system.solveEquilibrium();

	for(std::size_t i = 0; i < chemical_system.getSpecies().size(); i++) {
		ASSERT_FLOAT_EQ(
				chemical_system.getSpecies()[i]->concentration(),
				other_system.getSpecies()[i]->concentration()
				);
	}
}

TEST_F(BasicAqueousChemicalSystemTest, initial_aqueous_concentrations) {
	// User specified components
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("Na+").concentration(), 0.1*mol/l);
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("Cl-").concentration(), 0.1*mol/l);
	
	// Automatically added species
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("NaCl").concentration(), 0);
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("NaOH").concentration(), 0);
}

TEST_F(BasicAqueousChemicalSystemTest, activity) {
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("Na+").activity(), 0.1);
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("H+").activity(), std::pow(10, -7));
}

TEST_F(BasicAqueousChemicalSystemTest, H2O_activity) {
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("H2O").activity(), 1.0);
}

TEST_F(BasicAqueousChemicalSystemTest, mass_conservation_law) {
	std::vector<double> total_components(chemical_system.getComponents().size());
	chemical_system.massConservationLaw(total_components);

	double total_h = chemical_system.getComponent("H+").getTotalQuantity();
	double total_na = chemical_system.getComponent("Na+").getTotalQuantity();
	double total_cl = chemical_system.getComponent("Cl-").getTotalQuantity();

	const ChemicalSpecies& ho = chemical_system.getSpecies("OH-");
	const ChemicalSpecies& h = chemical_system.getSpecies("H+");
	const ChemicalSpecies& na = chemical_system.getSpecies("Na+");
	const ChemicalSpecies& cl = chemical_system.getSpecies("Cl-");
	const ChemicalSpecies& nacl = chemical_system.getSpecies("NaCl");
	const ChemicalSpecies& naoh = chemical_system.getSpecies("NaOH");

	ASSERT_FLOAT_EQ(
			total_components[chemical_system.getComponent("H+").getIndex()],
			h.quantity() - ho.quantity() - naoh.quantity() - total_h
			);
	ASSERT_FLOAT_EQ(
			total_components[chemical_system.getComponent("Na+").getIndex()],
			na.quantity() + nacl.quantity() + naoh.quantity() - total_na
			);
	ASSERT_FLOAT_EQ(
			total_components[chemical_system.getComponent("Cl-").getIndex()],
			cl.quantity() + nacl.quantity() - total_cl);

	std::minstd_rand rd;
	// Makes the system evolve randomly for N iterations, and check that the
	// total quantity of each component is always preserved according to the
	// massConservationLaw() implementation.
	for(std::size_t i = 0; i < 100; i++) {
		{
			// OH- random extent
			std::uniform_real_distribution<double> random_extent(
					-h.quantity(),
					ho.quantity()
					);
			chemical_system.proceed(
					chemical_system.getReaction("OH-"), random_extent(rd));
		}
		{
			// NaCl random extent
			std::uniform_real_distribution<double> random_extent(
					-std::min(na.quantity(), cl.quantity()),
					nacl.quantity()
					);
			chemical_system.proceed(
					chemical_system.getReaction("NaCl"), random_extent(rd));
		}
		{
			// NaOH random extent
			std::uniform_real_distribution<double> random_extent(
					-na.quantity(),
					naoh.quantity()
					);
			chemical_system.proceed(
					chemical_system.getReaction("NaOH"), random_extent(rd));
		}

		chemical_system.massConservationLaw(total_components);

		ASSERT_FLOAT_EQ(
				total_components[chemical_system.getComponent("H+").getIndex()],
				h.quantity() - ho.quantity() - naoh.quantity() - total_h
				);
		ASSERT_FLOAT_EQ(
				total_components[chemical_system.getComponent("Na+").getIndex()],
				na.quantity() + nacl.quantity() + naoh.quantity() - total_na
				);
		ASSERT_FLOAT_EQ(
				total_components[chemical_system.getComponent("Cl-").getIndex()],
				cl.quantity() + nacl.quantity() - total_cl);
	}
}

TEST_F(BasicMineralChemicalSystemTest, aqueous_and_mineral_species) {
	const std::regex surface_complex_regex("=S.*");

	ASSERT_THAT(chemical_system.getComponents(), UnorderedElementsAre(
				Pointee(Property(&ChemicalComponent::getSpecies,
						Pointee(Property(&ChemicalSpecies::getName, "H+"))
						)),
				Pointee(Property(&ChemicalComponent::getSpecies,
						Pointee(Property(&ChemicalSpecies::getName, "H2O"))
						)),
				Pointee(Property(&ChemicalComponent::getSpecies,
						Pointee(Property(&ChemicalSpecies::getName, "=SOH"))
						))
				));
	// Test fixed components handling
	for(auto& component : chemical_system.getComponents()) {
		if(component->isFixed()) {
			if (std::regex_match(component->getSpecies()->getName(), surface_complex_regex)) {
				ASSERT_THAT(
						component->getSpecies(),
						WhenDynamicCastTo<MineralSpecies*>(Not(IsNull())));
			} else {
				ASSERT_THAT(
						component.get()->getSpecies(),
						WhenDynamicCastTo<FixedAqueousSpecies*>(Not(IsNull()))
						);
			}
		} else {
			if (std::regex_match(component->getSpecies()->getName(), surface_complex_regex)) {
				ASSERT_THAT(
						component.get()->getSpecies(),
						WhenDynamicCastTo<MineralSpecies*>(Not(IsNull()))
						);
			} else {
				ASSERT_THAT(
						component.get()->getSpecies(),
						WhenDynamicCastTo<AqueousSpecies*>(Not(IsNull()))
						);
			}
		}
	}

	ASSERT_THAT(chemical_system.getSpecies(), UnorderedElementsAre(
				Pointee(Property(&ChemicalSpecies::getName, "H2O")),
				Pointee(Property(&ChemicalSpecies::getName, "H+")),
				Pointee(Property(&ChemicalSpecies::getName, "OH-")),
				Pointee(Property(&ChemicalSpecies::getName, "=SOH")),
				Pointee(Property(&ChemicalSpecies::getName, "=SOH2"))
				));
	for(auto& species : chemical_system.getSpecies()) {
		if(species->getName() == "H2O") {
			ASSERT_THAT(species.get(), WhenDynamicCastTo<Solvent*>(Not(IsNull())));
		} else if (std::regex_match(species->getName(), surface_complex_regex)) {
			ASSERT_THAT(species.get(), AnyOf(
						WhenDynamicCastTo<MineralSpecies*>(Not(IsNull())),
						WhenDynamicCastTo<FixedMineralSpecies*>(Not(IsNull()))
						));
		} else {
			ASSERT_THAT(species.get(), AnyOf(
					WhenDynamicCastTo<AqueousSpecies*>(Not(IsNull())),
					WhenDynamicCastTo<FixedAqueousSpecies*>(Not(IsNull()))
					));
		}
	}
}

TEST_F(BasicMineralChemicalSystemTest, initial_mineral_concentrations) {
	ASSERT_FLOAT_EQ(
			chemical_system.sitesQuantity(), 
			2.5 * g/l * 24.2 * m2/g * 0.8 * entities/nm2
			);
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("=SOH2").concentration(), 0);

	ASSERT_FLOAT_EQ(chemical_system.getSpecies("=SOH").concentration(), 1.0);
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("=SOH").activity(), 1.0);
	ASSERT_FLOAT_EQ(
			chemical_system.getSpecies("=SOH").quantity(),
			chemical_system.sitesQuantity()	
			);
}

TEST_F(BasicMineralChemicalSystemTest, proceed_mineral_reaction) {
	const double init_quantity = chemical_system.getSpecies("=SOH").quantity();

	chemical_system.proceed(chemical_system.getReaction("=SOH2"), -1e-6);

	ASSERT_FLOAT_EQ(
			chemical_system.getSpecies("=SOH").quantity(),
			init_quantity - 1e-6
			);
	ASSERT_FLOAT_EQ(
			chemical_system.getSpecies("=SOH").concentration(),
			(init_quantity - 1e-6) / init_quantity
			);
	ASSERT_FLOAT_EQ(
			chemical_system.getSpecies("=SOH").activity(),
			chemical_system.getSpecies("=SOH").concentration()
			);

	ASSERT_FLOAT_EQ(
			chemical_system.getSpecies("=SOH2").quantity(),
			1e-6
			);
	ASSERT_FLOAT_EQ(
			chemical_system.getSpecies("=SOH2").concentration(),
			1e-6 / init_quantity
			);
	ASSERT_FLOAT_EQ(
			chemical_system.getSpecies("=SOH2").activity(),
			chemical_system.getSpecies("=SOH2").concentration()
			);
}

class BadMineralChemicalSystemTest : public Test {
	protected:
		ChemicalSystem chemical_system;

		void SetUp() override {
			chemical_system.addReaction("OH-", -13.997, {
					{"OH-", -1},
					{"H+", -1},
					{"H2O", 1}
					});
			chemical_system.addSolvent("H2O");
			chemical_system.fixPH(7);

			chemical_system.addReaction("=SOH2", 3.46, {
					{"=SOH2", MINERAL, -1},
					{"=SOH", MINERAL, 1}, // AQUEOUS specified here only for test
										  // purpose, since it should be the default
										  // phase
					{"H+", AQUEOUS, 1}
					});
		}
};

TEST_F(BadMineralChemicalSystemTest, add_component) {
	EXPECT_THROW(
			chemical_system.addComponent("=SOH", MINERAL, 1.0),
			InvalidMineralSpeciesWithUndefinedSitesCount
			);
	try {
			chemical_system.addComponent("=SOH", MINERAL, 1.0);
	} catch(const InvalidMineralSpeciesWithUndefinedSitesCount& e) {
		CHEM_LOGV(6) << "Exception: " << e.what();
		ASSERT_THAT(e.getChemicalSystem(), Ref(chemical_system));
		ASSERT_EQ(e.getName(), "=SOH");
		ASSERT_EQ(e.getPhase(), MINERAL);
	}
}

TEST_F(BadMineralChemicalSystemTest, fix_component) {
	EXPECT_THROW(
			chemical_system.fixComponent("=SOH", MINERAL, 1.0),
			InvalidMineralSpeciesWithUndefinedSitesCount
			);
	try {
			chemical_system.fixComponent("=SOH", MINERAL, 1.0);
	} catch(const InvalidMineralSpeciesWithUndefinedSitesCount& e) {
		CHEM_LOGV(6) << "Exception: " << e.what();
		ASSERT_THAT(e.getChemicalSystem(), Ref(chemical_system));
		ASSERT_EQ(e.getName(), "=SOH");
		ASSERT_EQ(e.getPhase(), MINERAL);
	}
}


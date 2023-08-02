#include "gmock/gmock.h"
#include "chemmisol/chemical/system.h"
#include <cstddef>
#include <regex>
#include <random>

namespace chemmisol {
	std::ostream& operator<<(std::ostream& o, const ChemicalSpecies& s) {
		o << "{i:" << s.getIndex() << ", name:" << s.getName() << ", phase:" << s.getPhase() << "}";
		return o;
	}
}

using namespace testing;
using namespace chemmisol;

class ChemicalSystemTest : public Test {
	protected:
	ChemicalSystem chemical_system;

	void SetUp() override {
		chemical_system.addReaction("OH-", -13.997, {
				{"OH-", -1},
				{"H+", -1},
				{"H2O", 1}
				});
		chemical_system.addReaction("NaCl", -0.3, {
				{"NaCl", -1},
				{"Na+", AQUEOUS, 1}, // AQUEOUS specified here only for test
									 // purpose, since it should be the default
									 // phase
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
		// Automatically adds the H+ component
		chemical_system.fixPH(7);

		chemical_system.setUp();
	}

	std::array<double, 7> analytical_concentrations(const std::array<double, 3>& x) const {
		const double x_0 = x[chemical_system.getReaction("OH-").getIndex()];
		const double x_1 = x[chemical_system.getReaction("NaCl").getIndex()];
		const double x_2 = x[chemical_system.getReaction("NaOH").getIndex()];

		double h = chemical_system.getSpecies("H+").concentration();
		double ho = chemical_system.getSpecies("OH-").concentration();
		double nacl = chemical_system.getSpecies("NaCl").concentration();
		double na = chemical_system.getSpecies("Na+").concentration();
		double cl = chemical_system.getSpecies("Cl-").concentration();
		double naoh = chemical_system.getSpecies("NaOH").concentration();
		const double V = AqueousSpecies::V;

		//h += (-x_0-x_2)/V;
		h += 0; // Fixed
		ho += -x_0/V;
		nacl += -x_1/V;
		na += (x_1+x_2)/V;
		cl += x_1/V;
		naoh += -x_2/V;

		std::array<double, 7> concentrations;
		concentrations[chemical_system.getSpecies("H+").getIndex()] = h;
		concentrations[chemical_system.getSpecies("OH-").getIndex()] = ho;
		concentrations[chemical_system.getSpecies("NaCl").getIndex()] = nacl;
		concentrations[chemical_system.getSpecies("Na+").getIndex()] = na;
		concentrations[chemical_system.getSpecies("Cl-").getIndex()] = cl;
		concentrations[chemical_system.getSpecies("NaOH").getIndex()] = naoh;

		return concentrations;
	}
	std::array<double, 3> analytical_f(const std::array<double, 3>& x) const {
		auto concentrations = analytical_concentrations(x);
		
		const double h = concentrations[chemical_system.getSpecies("H+").getIndex()];
		const double ho = concentrations[chemical_system.getSpecies("OH-").getIndex()];
		const double nacl = concentrations[chemical_system.getSpecies("NaCl").getIndex()];
		const double na = concentrations[chemical_system.getSpecies("Na+").getIndex()];
		const double cl = concentrations[chemical_system.getSpecies("Cl-").getIndex()];
		const double naoh = concentrations[chemical_system.getSpecies("NaOH").getIndex()];

		std::array<double, 3> f_x;
		f_x[chemical_system.getReaction("OH-").getIndex()]
			= -std::log10(h/(1*mol/l)) - std::log10(ho/(1*mol/l)) - 13.997;
		f_x[chemical_system.getReaction("NaCl").getIndex()]
			= - std::log10(nacl) + std::log10(na) + std::log10(cl) - 0.3;
		f_x[chemical_system.getReaction("NaOH").getIndex()]
			= - std::log10(naoh) - std::log10(h) + std::log10(na) - 13.897;
		return f_x;
	}

	std::array<std::array<double, 3>, 3> analytical_df(const std::array<double, 3>& x) const {
		std::array<std::array<double, 3>, 3> df_x;
		auto concentrations = analytical_concentrations(x);
		
		const double ho = concentrations[chemical_system.getSpecies("OH-").getIndex()];
		const double nacl = concentrations[chemical_system.getSpecies("NaCl").getIndex()];
		const double na = concentrations[chemical_system.getSpecies("Na+").getIndex()];
		const double cl = concentrations[chemical_system.getSpecies("Cl-").getIndex()];
		const double naoh = concentrations[chemical_system.getSpecies("NaOH").getIndex()];

		// df HO-
		{
			auto& df_h = df_x[chemical_system.getReaction("OH-").getIndex()];
			df_h[chemical_system.getReaction("OH-").getIndex()] = 1/ho*1/ln10;
			df_h[chemical_system.getReaction("NaCl").getIndex()] = 0;
			df_h[chemical_system.getReaction("NaOH").getIndex()] = 0;
		}

		// df NaCl
		{
			auto& df_nacl = df_x[chemical_system.getReaction("NaCl").getIndex()];
			df_nacl[chemical_system.getReaction("OH-").getIndex()] = 0;
			df_nacl[chemical_system.getReaction("NaCl").getIndex()] = (1/nacl + 1/na + 1/cl) * 1/ln10;
			df_nacl[chemical_system.getReaction("NaOH").getIndex()] = 1/na * 1/ln10;
		}

		// df NaOH
		{
			auto& df_naoh = df_x[chemical_system.getReaction("NaOH").getIndex()];
			df_naoh[chemical_system.getReaction("OH-").getIndex()] = 0;
			df_naoh[chemical_system.getReaction("NaCl").getIndex()] = 1/na * 1/ln10;
			df_naoh[chemical_system.getReaction("NaOH").getIndex()] = (1/naoh + 1/na) * 1/ln10;
		}
		return df_x;
	}

	void checkEquilibrium() {
		{
			// Check OH-
			double H2O = chemical_system.getSpecies("H2O").activity();
			double H = chemical_system.getSpecies("H+").activity();
			double OH = chemical_system.getSpecies("OH-").activity();
			ASSERT_FLOAT_EQ(
					std::log10((H*OH)/H2O),
					chemical_system.getReaction("OH-").getLogK()
					);
		}
		{
			// Check NaCl
			double NaCl = chemical_system.getSpecies("NaCl").activity();
			double Na = chemical_system.getSpecies("Na+").activity();
			double Cl = chemical_system.getSpecies("Cl-").activity();
			ASSERT_FLOAT_EQ(
					std::log10(NaCl/(Na*Cl)),
					chemical_system.getReaction("NaCl").getLogK()
					);
		}
		{
			// Check NaOH
			double NaOH = chemical_system.getSpecies("NaOH").activity();
			double Na = chemical_system.getSpecies("Na+").activity();
			double H = chemical_system.getSpecies("H+").activity();
			ASSERT_FLOAT_EQ(
					std::log10((NaOH*H)/Na),
					chemical_system.getReaction("NaOH").getLogK()
					);
		}
		{
			// Check pH
			double H = chemical_system.getSpecies("H+").activity();
			ASSERT_FLOAT_EQ(-std::log10(H), 7);
		}
	}
};

TEST_F(ChemicalSystemTest, get_components) {
	auto& components = chemical_system.getComponents();

	for(std::size_t i = 0; i < components.size(); i++)
		ASSERT_THAT(components[i]->getIndex(), i);
}

TEST_F(ChemicalSystemTest, get_species) {
	auto& species = chemical_system.getSpecies();

	for(std::size_t i = 0; i < species.size(); i++)
		ASSERT_THAT(species[i]->getIndex(), i);
}

TEST_F(ChemicalSystemTest, aqueous_species) {
	ASSERT_THAT(chemical_system.getComponents(), UnorderedElementsAre(
				Pointee(AllOf(
						Property(&Component::getSpecies,
							Pointee(Property(&ChemicalSpecies::getName, "H+"))
						),
						Property(&Component::getTotalQuantity,
							std::pow(10, -7) * AqueousSpecies::V
						)
						)),
				Pointee(AllOf(
						Property(&Component::getSpecies,
							Pointee(Property(&ChemicalSpecies::getName, "Na+"))
						),
						Property(&Component::getTotalQuantity, 0.1*mol/l)
						)),
				Pointee(AllOf(
						Property(&Component::getSpecies,
						Pointee(Property(&ChemicalSpecies::getName, "Cl-"))
						),
						Property(&Component::getTotalQuantity, 0.1*mol/l)
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

TEST_F(ChemicalSystemTest, initial_aqueous_concentrations) {
	// User specified components
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("Na+").concentration(), 0.1*mol/l);
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("Cl-").concentration(), 0.1*mol/l);
	
	// Automatically added species
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("NaCl").concentration(), 0);
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("NaOH").concentration(), 0);
}

TEST_F(ChemicalSystemTest, activity) {
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("Na+").activity(), 0.1);
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("H+").activity(), std::pow(10, -7));
}

TEST_F(ChemicalSystemTest, H2O_activity) {
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("H2O").activity(), 1.0);
}

TEST_F(ChemicalSystemTest, mass_conservation_law) {
	std::vector<double> total_components(chemical_system.getComponents().size());
	chemical_system.massConservationLaw(total_components);

	double total_h = chemical_system.getComponent("H+").getTotalQuantity();
	double total_na = chemical_system.getComponent("Na+").getTotalQuantity();
	double total_cl = chemical_system.getComponent("Cl-").getTotalQuantity();

	ASSERT_FLOAT_EQ(
			total_components[chemical_system.getComponent("H+").getIndex()], total_h);
	ASSERT_FLOAT_EQ(
			total_components[chemical_system.getComponent("Na+").getIndex()], total_na);
	ASSERT_FLOAT_EQ(
			total_components[chemical_system.getComponent("Cl-").getIndex()], total_cl);

	const ChemicalSpecies& ho = chemical_system.getSpecies("OH-");
	const ChemicalSpecies& h = chemical_system.getSpecies("H+");
	const ChemicalSpecies& na = chemical_system.getSpecies("Na+");
	const ChemicalSpecies& cl = chemical_system.getSpecies("Cl-");
	const ChemicalSpecies& nacl = chemical_system.getSpecies("NaCl");
	const ChemicalSpecies& naoh = chemical_system.getSpecies("NaOH");

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

		ASSERT_FLOAT_EQ(
				total_components[chemical_system.getComponent("H+").getIndex()], total_h);
		ASSERT_FLOAT_EQ(
				total_components[chemical_system.getComponent("Na+").getIndex()], total_na);
		ASSERT_FLOAT_EQ(
				total_components[chemical_system.getComponent("Cl-").getIndex()], total_cl);
	}
}

TEST_F(ChemicalSystemTest, guess_extent) {
	auto guesses = chemical_system.guessInitialExtents();

	// Checks that the initial guess allows to solve the equilibrium
	chemical_system.solveEquilibrium();
	checkEquilibrium();
}

TEST_F(ChemicalSystemTest, guess_extent_similar_logK) {
	// Adding a fake reaction that consumes Na+ with a logK similar to the NaOH
	// reaction
	chemical_system.addReaction("NaOH2", -14, {
			{"NaOH2", -1},
			{"H+", -2},
			{"Na+", 1},
			{"H2O", 1}
			});
	chemical_system.setUp();

	auto guesses = chemical_system.guessInitialExtents();

	// Checks that the initial guess allows to solve the equilibrium
	chemical_system.solveEquilibrium();
	checkEquilibrium();
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_matrix) {
	auto reaction_matrix = chemical_system.getReactionMatrix();

	ASSERT_THAT(reaction_matrix, SizeIs(3));
	{
		// Check OH-
		const std::vector<double>& reaction
			= reaction_matrix[chemical_system.getReaction("OH-").getIndex()];
		ASSERT_THAT(reaction[chemical_system.getSpecies("OH-").getIndex()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getSpecies("H+").getIndex()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getSpecies("H2O").getIndex()], Eq(1));
	}
	{
		// Check NaCl
		const std::vector<double>& reaction
			= reaction_matrix[chemical_system.getReaction("NaCl").getIndex()];
		ASSERT_THAT(reaction[chemical_system.getSpecies("NaCl").getIndex()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getSpecies("Na+").getIndex()], Eq(1));
		ASSERT_THAT(reaction[chemical_system.getSpecies("Cl-").getIndex()], Eq(1));
	}
	{
		// Check NaOH
		const std::vector<double>& reaction
			= reaction_matrix[chemical_system.getReaction("NaOH").getIndex()];
		ASSERT_THAT(reaction[chemical_system.getSpecies("NaOH").getIndex()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getSpecies("H+").getIndex()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getSpecies("Na+").getIndex()], Eq(1));
		ASSERT_THAT(reaction[chemical_system.getSpecies("H2O").getIndex()], Eq(1));
	}
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_concentrations) {
	using namespace solver;
	F f(chemical_system);

	std::array<double, 3> x;
	x[chemical_system.getReaction("OH-").getIndex()] = 0;
	x[chemical_system.getReaction("NaCl").getIndex()] = -1e-16;
	x[chemical_system.getReaction("NaOH").getIndex()] = -1e-16;
	auto C = f.concentrations({x.begin(), x.end()});
	auto analytical_C = analytical_concentrations(x);
	for (auto& reaction : chemical_system.getReactions()) {
		ASSERT_FLOAT_EQ(C[reaction->getIndex()], analytical_C[reaction->getIndex()]);
	}
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_f) {
	using namespace solver;
	F f(chemical_system);

	std::array<double, 3> x;
	x[chemical_system.getReaction("OH-").getIndex()] = -1e-7;
	x[chemical_system.getReaction("NaCl").getIndex()] = -1e-16;
	x[chemical_system.getReaction("NaOH").getIndex()] = -1e-16;
	auto F = f.f({x.begin(), x.end()});
	auto analytical_F = analytical_f(x);
	for (auto& reaction : chemical_system.getReactions()) {
		ASSERT_FLOAT_EQ(F[reaction->getIndex()], analytical_F[reaction->getIndex()]);
	}
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_df) {
	using namespace solver;
	F f(chemical_system);

	std::array<double, 3> x;
	x[chemical_system.getReaction("OH-").getIndex()] = -1e-7;
	x[chemical_system.getReaction("NaCl").getIndex()] = -1e-16;
	x[chemical_system.getReaction("NaOH").getIndex()] = -1e-16;
	auto dF = f.df({x.begin(), x.end()});
	auto analytical_dF = analytical_df(x);
	for (auto& reaction : chemical_system.getReactions()) {
		for (auto& extent : chemical_system.getReactions()) {
			ASSERT_FLOAT_EQ(
					dF[reaction->getIndex()][extent->getIndex()],
					analytical_dF[reaction->getIndex()][extent->getIndex()]);
		}
	}
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction) {
	auto reaction_matrix = chemical_system.getReactionMatrix();

	/*
	 *chemical_system
	 *    .setInitialGuessExtent("NaCl", -std::pow(10, -2.3));
	 *chemical_system
	 *    .setInitialGuessExtent("NaOH", -std::pow(10, -7.897));
	 */

	chemical_system.setMaxIteration(10);
	chemical_system.solveEquilibrium();

	checkEquilibrium();

	CHEM_LOG(INFO) << "Basic NaCl reaction equilibrium test:";
	for(auto& species : chemical_system.getSpecies())
		CHEM_LOG(INFO) << "  " << species->getName() << ": " << species->concentration()/(1*mol/l) << " mol/l";
}

TEST_F(ChemicalSystemTest, test_pH) {
	auto reaction_matrix = chemical_system.getReactionMatrix();
	chemical_system.setMaxIteration(200);
	chemical_system.solveEquilibrium();

	// Check pH
	ASSERT_FLOAT_EQ(chemical_system.getPH(), 7);
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_quotient) {
	auto reaction_matrix = chemical_system.getReactionMatrix();
	chemical_system.setMaxIteration(200);
	chemical_system.solveEquilibrium();

	{
		// Check OH-
		double H2O = chemical_system.getSpecies("H2O").activity();
		double H = chemical_system.getSpecies("H+").activity();
		double OH = chemical_system.getSpecies("OH-").activity();
		ASSERT_FLOAT_EQ(
				chemical_system.reactionQuotient("OH-"),
				(H*OH)/H2O
				);
	}
	{
		// Check NaCl
		double NaCl = chemical_system.getSpecies("NaCl").activity();
		double Na = chemical_system.getSpecies("Na+").activity();
		double Cl = chemical_system.getSpecies("Cl-").activity();
		ASSERT_FLOAT_EQ(
				chemical_system.reactionQuotient("NaCl"),
				NaCl/(Na*Cl)
				);
	}
	{
		// Check NaOH
		double NaOH = chemical_system.getSpecies("NaOH").activity();
		double Na = chemical_system.getSpecies("Na+").activity();
		double H = chemical_system.getSpecies("H+").activity();
		ASSERT_FLOAT_EQ(
				chemical_system.reactionQuotient("NaOH"),
				(NaOH*H)/Na
				);
	}
}

TEST_F(ChemicalSystemTest, copy_constructor) {
	ChemicalSystem other_system(chemical_system);

	ASSERT_THAT(
			chemical_system.getComponents(),
			SizeIs(other_system.getComponents().size())
			);
	for(std::size_t i = 0; i < chemical_system.getComponents().size(); i++) {
		const Component& component = *chemical_system.getComponents()[i];
		const Component& other_component = *other_system.getComponents()[i];
		ASSERT_THAT(component,
				AllOf(
					Property(&Component::getIndex, other_component.getIndex()),
					Property(&Component::isFixed, other_component.isFixed()),
					Property(&Component::getSpecies, Pointee(
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

class AdsorptionTest : public Test {
	protected:
	ChemicalSystem chemical_system {
			2.5 * g/l,
			24.2 * m2/g,
			0.8 * entities/nm2,
			"=SOH"
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
		
		chemical_system.fixPH(7);

		chemical_system.setUp();
	}
};

TEST_F(AdsorptionTest, aqueous_and_mineral_species) {
	std::regex surface_complex_regex("=S.*");
	ASSERT_THAT(chemical_system.getComponents(), UnorderedElementsAre(
				Pointee(Property(&Component::getSpecies,
						Pointee(Property(&ChemicalSpecies::getName, "H+"))
						)),
				Pointee(Property(&Component::getSpecies,
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

TEST_F(AdsorptionTest, initial_mineral_concentrations) {
	// Automatically added components
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("=SOH").concentration(), 1.0);
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("=SOH2").concentration(), 0);
}


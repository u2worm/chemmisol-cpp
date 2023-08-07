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

void print_df(
		const ChemicalSystem& system,
		const solver::F& f,
		std::vector<std::vector<double>>& df,
		std::size_t unfixed_components) {

	std::ostringstream s;
	s << "  " << std::setw(6) << "";
	for(const auto& species : system.getSpecies()) {
		if(species->getName() != "H2O"
				&& species->getName() != "H+") {
			s << std::setw(16) << species->getName();
		}
	}
	CHEM_LOGV(5) << s.str();
	for(const auto& component : system.getComponents()) {
		s.str("");
		if(component->getSpecies()->getName() != "H2O"
				&& component->getSpecies()->getName() != "H+") {
			s << "  " << std::setw(6) << component->getSpecies()->getName();
			for(const auto& item : df[f.componentsIndexes()[component->getIndex()]])
				s << std::setw(16) << std::scientific << item;
			CHEM_LOGV(5) << s.str();
		}
	}
	for(const auto& reaction : system.getReactions()) {
		s.str("");
		s << "  " << std::setw(6) << reaction->getName();
			for(const auto& item : df[unfixed_components + reaction->getIndex()])
				s << std::setw(16) << std::scientific << item;
		CHEM_LOGV(5) << s.str();
	}
}

/**
 * Tests basic ChemicalSystem features and solving of a simple reaction system
 * with only unitary stoichiometric coefficients.
 */
class BasicChemicalSystemTest : public Test {
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
		chemical_system.addSolvent("H2O");

		// Automatically adds the H+ component
		chemical_system.fixPH(7);

		chemical_system.setUp();
	}

	std::vector<double> analytical_f(
			const std::vector<std::size_t>& components_indexes,
			const std::vector<std::size_t>& species_indexes,
			const std::vector<double>& activities) const {
		// Fixed component
		const double h = chemical_system.getSpecies("H+").activity();
		// Dynamic components
		const double ho = activities[species_indexes[chemical_system.getSpecies("OH-").getIndex()]];
		const double nacl = activities[species_indexes[chemical_system.getSpecies("NaCl").getIndex()]];
		const double na = activities[species_indexes[chemical_system.getSpecies("Na+").getIndex()]];
		const double cl = activities[species_indexes[chemical_system.getSpecies("Cl-").getIndex()]];
		const double naoh = activities[species_indexes[chemical_system.getSpecies("NaOH").getIndex()]];

		// Count of mass conservation law equations
		const std::size_t N = 2;

		std::vector<double> f_x(
				N
				// Count of reactions
				+ chemical_system.getReactions().size());
		// Mass conservation law
		f_x[components_indexes[chemical_system.getComponent("Na+").getIndex()]]
			= (na + nacl + naoh) - chemical_system.getComponent("Na+").getTotalQuantity();
		f_x[components_indexes[chemical_system.getComponent("Cl-").getIndex()]]
			= (cl + nacl) - chemical_system.getComponent("Cl-").getTotalQuantity();
		// H+ and H2O are fixed

		// Equilibriums
		f_x[N+chemical_system.getReaction("OH-").getIndex()]
			= h * ho - chemical_system.getReaction("OH-").getK() * 1.0;
		f_x[N+chemical_system.getReaction("NaCl").getIndex()]
			= nacl - chemical_system.getReaction("NaCl").getK() * na * cl;
		f_x[N+chemical_system.getReaction("NaOH").getIndex()]
			= naoh * h - chemical_system.getReaction("NaOH").getK() * na;
		return f_x;
	}

	std::vector<std::vector<double>> analytical_df(
			const std::vector<std::size_t>& components_indexes,
			const std::vector<std::size_t>& species_indexes,
			const std::vector<double>& activities) const {
		std::vector<std::vector<double>> df_x(
				chemical_system.getComponents().size()
				+ chemical_system.getReactions().size());
		// Fixed component
		const double h = chemical_system.getSpecies("H+").activity();
		// Dynamic components
		const double na = activities[species_indexes[chemical_system.getSpecies("Na+").getIndex()]];
		const double cl = activities[species_indexes[chemical_system.getSpecies("Cl-").getIndex()]];

		// Count of fixed species (H+ and H2O)
		const std::size_t N = 2;
		// Count of mass conservation law equations
		const std::size_t M = chemical_system.getComponents().size()-N;
		// Mass conservation law
		{
			{
				// Na+
				auto& df = df_x[components_indexes[chemical_system.getComponent("Na+").getIndex()]];
				df.resize(chemical_system.getSpecies().size()-N);
				df[species_indexes[chemical_system.getSpecies("Na+").getIndex()]] = 1.0;
				df[species_indexes[chemical_system.getSpecies("NaCl").getIndex()]] = 1.0;
				df[species_indexes[chemical_system.getSpecies("NaOH").getIndex()]] = 1.0;
				// df = 0.0 for other species
			}
			{
				// Cl-
				auto& df = df_x[chemical_system.getComponent("Cl-").getIndex()];
				df.resize(chemical_system.getSpecies().size()-N);
				df[species_indexes[chemical_system.getSpecies("Cl-").getIndex()]] = 1.0;
				df[species_indexes[chemical_system.getSpecies("NaCl").getIndex()]] = 1.0;
				// df = 0.0 for other species
			}
		}

		// df HO-
		{
			auto& df_h = df_x[M+chemical_system.getReaction("OH-").getIndex()];
			df_h.resize(chemical_system.getSpecies().size()-N);
			df_h[species_indexes[chemical_system.getSpecies("OH-").getIndex()]] = h;
		}

		// df NaCl
		{
			double K = chemical_system.getReaction("NaCl").getK();
			auto& df_nacl = df_x[M+chemical_system.getReaction("NaCl").getIndex()];
			df_nacl.resize(chemical_system.getSpecies().size()-N);
			df_nacl[species_indexes[chemical_system.getSpecies("NaCl").getIndex()]] = 1;
			df_nacl[species_indexes[chemical_system.getSpecies("Na+").getIndex()]] = -K * cl;
			df_nacl[species_indexes[chemical_system.getSpecies("Cl-").getIndex()]] = -K * na;
		}

		// df NaOH
		{
			double K = chemical_system.getReaction("NaOH").getK();
			auto& df_naoh = df_x[M+chemical_system.getReaction("NaOH").getIndex()];
			df_naoh.resize(chemical_system.getSpecies().size()-N);
			df_naoh[species_indexes[chemical_system.getSpecies("NaOH").getIndex()]] = h;
			df_naoh[species_indexes[chemical_system.getSpecies("Na+").getIndex()]] = -K;
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

TEST_F(BasicChemicalSystemTest, get_components) {
	auto& components = chemical_system.getComponents();

	for(std::size_t i = 0; i < components.size(); i++)
		ASSERT_THAT(components[i]->getIndex(), i);
}

TEST_F(BasicChemicalSystemTest, get_species) {
	auto& species = chemical_system.getSpecies();

	for(std::size_t i = 0; i < species.size(); i++)
		ASSERT_THAT(species[i]->getIndex(), i);
}

TEST_F(BasicChemicalSystemTest, get_component_reagent) {
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

TEST_F(BasicChemicalSystemTest, get_species_reagents) {
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

TEST_F(BasicChemicalSystemTest, aqueous_species) {
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
						)),
				Pointee(AllOf(
						Property(&Component::getSpecies,
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

TEST_F(BasicChemicalSystemTest, copy_constructor) {
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

TEST_F(BasicChemicalSystemTest, initial_aqueous_concentrations) {
	// User specified components
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("Na+").concentration(), 0.1*mol/l);
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("Cl-").concentration(), 0.1*mol/l);
	
	// Automatically added species
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("NaCl").concentration(), 0);
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("NaOH").concentration(), 0);
}

TEST_F(BasicChemicalSystemTest, activity) {
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("Na+").activity(), 0.1);
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("H+").activity(), std::pow(10, -7));
}

TEST_F(BasicChemicalSystemTest, H2O_activity) {
	ASSERT_FLOAT_EQ(chemical_system.getSpecies("H2O").activity(), 1.0);
}

TEST_F(BasicChemicalSystemTest, mass_conservation_law) {
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

TEST_F(BasicChemicalSystemTest, basic_NaCl_reaction_matrix) {
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

TEST_F(BasicChemicalSystemTest, basic_NaCl_reaction_f) {
	using namespace solver;
	F f(chemical_system);

	std::vector<double> activities(chemical_system.getSpecies().size()-2);
	{
		ChemicalSystem fake_system(chemical_system);
		// Example extents
		fake_system.proceed(chemical_system.getReaction("OH-"), -1e-7);
		fake_system.proceed(chemical_system.getReaction("NaCl"), -1e-3);
		fake_system.proceed(chemical_system.getReaction("NaOH"), -1e-4);
		for(const auto& species : fake_system.getSpecies()) {
			if(species->getName() != "H2O" && species->getName() != "H+") {
				activities[f.speciesIndexes()[species->getIndex()]]
					= species->activity();
			}
		}
	}

	auto F = f.f(activities);
	auto analytical_F = analytical_f(
			f.componentsIndexes(),
			f.speciesIndexes(),
			activities);
	for (auto& component : chemical_system.getComponents()) {
		if(component->getSpecies()->getName() != "H2O"
				&& component->getSpecies()->getName() != "H+") {
			CHEM_LOGV(5) << "Checking mass conservation of component: "
				<< component->getSpecies()->getName();
			ASSERT_FLOAT_EQ(
					F[f.speciesIndexes()[component->getSpecies()->getIndex()]],
					analytical_F[f.speciesIndexes()[component->getSpecies()->getIndex()]]);
		}
	}
	for (auto& reaction : chemical_system.getReactions()) {
		CHEM_LOGV(5) << "Checking distance to equilibrium of reaction: "
			<< reaction->getName() << " (logK=" << reaction->getLogK() << ")";
		ASSERT_FLOAT_EQ(F[2+reaction->getIndex()], analytical_F[2+reaction->getIndex()]);
	}
}

TEST_F(BasicChemicalSystemTest, basic_NaCl_reaction_df) {
	using namespace solver;
	F f(chemical_system);

	std::vector<double> activities(chemical_system.getSpecies().size()-2);
	{
		ChemicalSystem fake_system(chemical_system);
		// Example extents
		fake_system.proceed(chemical_system.getReaction("OH-"), -1e-7);
		fake_system.proceed(chemical_system.getReaction("NaCl"), -1e-3);
		fake_system.proceed(chemical_system.getReaction("NaOH"), -1e-4);
		for(const auto& species : fake_system.getSpecies()) {
			if(species->getName() != "H2O" && species->getName() != "H+") {
				activities[f.speciesIndexes()[species->getIndex()]]
					= species->activity();
			}
		}
	}
	
	for(const auto& species : chemical_system.getSpecies()) {
		if(species->getName() != "H2O" && species->getName() != "H+") {
			activities[f.speciesIndexes()[species->getIndex()]] = species->activity();
		}
	}

	auto dF = f.df(activities);
	auto analytical_dF = analytical_df(
			f.componentsIndexes(),
			f.speciesIndexes(),
			activities);

	CHEM_LOGV(5) << "Computed jacobian:";
	print_df(chemical_system, f, dF, 2);

	CHEM_LOGV(5) << "Expected jacobian:";
	print_df(chemical_system, f, analytical_dF, 2);

	for(const auto& component : chemical_system.getComponents()) {
		if(component->getSpecies()->getName() != "H2O"
				&& component->getSpecies()->getName() != "H+") {
			CHEM_LOGV(6) << "F: mass_conservation_law of "
				<< component->getSpecies()->getName();
			for(const auto& species : chemical_system.getSpecies()) {
				if(species->getName() != "H2O"
						&& species->getName() != "H+") {
					CHEM_LOGV(6) << "  dX: " << species->getName();
					ASSERT_FLOAT_EQ(
							dF
							[f.componentsIndexes()[component->getIndex()]]
							[f.speciesIndexes()[species->getIndex()]],
							analytical_dF
							[f.componentsIndexes()[component->getIndex()]]
							[f.speciesIndexes()[species->getIndex()]]
							);
				}
			}
		}
	}
	for (auto& reaction : chemical_system.getReactions()) {
		CHEM_LOGV(6) << "F: equilibrium of reaction "
			<< reaction->getName();
		for(const auto& species : chemical_system.getSpecies()) {
			if(species->getName() != "H2O" && species->getName() != "H+") {
				CHEM_LOGV(6) << "  dX: " << species->getName();
				ASSERT_FLOAT_EQ(
						dF
						[2 + reaction->getIndex()]
						[f.speciesIndexes()[species->getIndex()]],
						analytical_dF
						[2 + reaction->getIndex()]
						[f.speciesIndexes()[species->getIndex()]]);
			}
		}
	}
}

TEST_F(BasicChemicalSystemTest, basic_NaCl_reaction) {
	auto reaction_matrix = chemical_system.getReactionMatrix();

	chemical_system.setMaxIteration(10);
	chemical_system.solveEquilibrium();

	checkEquilibrium();

	CHEM_LOG(INFO) << "Basic NaCl reaction equilibrium test:";
	for(auto& species : chemical_system.getSpecies())
		CHEM_LOG(INFO) << "  " << species->getName() << ": " << species->concentration()/(1*mol/l) << " mol/l";
}

TEST_F(BasicChemicalSystemTest, test_pH) {
	auto reaction_matrix = chemical_system.getReactionMatrix();
	chemical_system.setMaxIteration(200);
	chemical_system.solveEquilibrium();

	// Check pH
	ASSERT_FLOAT_EQ(chemical_system.getPH(), 7);
}

TEST_F(BasicChemicalSystemTest, basic_NaCl_reaction_quotient) {
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

/**
 * Tests ChemicalSystem equilibrium solving features with a reaction system than
 * BasicChemicalSystemTest, with non unitary stoichiometric coefficients.
 */
class ComplexChemicalSystemTest : public Test {
	protected:
	ChemicalSystem chemical_system;

	void SetUp() override {
		chemical_system.addReaction("OH-", -13.997, {
				{"OH-", -1},
				{"H+", -1},
				{"H2O", 1}
				});

		chemical_system.addReaction("H3PO4", -10.5, {
				{"H3PO4", -1},
				{"PO4-3", 1},
				{"H+", 3}
				});

		chemical_system.addComponent("PO4-3", 6.7e-6*mol/l);
		chemical_system.addSolvent("H2O");

		// Automatically adds the H+ component
		chemical_system.fixPH(7);

		chemical_system.setUp();
	}

	std::vector<double> analytical_f(
			const std::vector<std::size_t>& components_indexes,
			const std::vector<std::size_t>& species_indexes,
			const std::vector<double>& activities) const {
		// Fixed component
		const double h = chemical_system.getSpecies("H+").activity();
		// Dynamic components
		const double ho = activities[species_indexes[chemical_system.getSpecies("OH-").getIndex()]];
		const double h3po4 = activities[species_indexes[chemical_system.getSpecies("H3PO4").getIndex()]];
		const double po4 = activities[species_indexes[chemical_system.getSpecies("PO4-3").getIndex()]];

		// Count of mass conservation law equations
		const std::size_t N = 1;

		std::vector<double> f_x(
				N
				// Count of reactions
				+ chemical_system.getReactions().size());
		// Mass conservation law
		f_x[components_indexes[chemical_system.getComponent("PO4-3").getIndex()]]
			= (po4 + h3po4) - chemical_system.getComponent("PO4-3").getTotalQuantity();
		// H+ and H2O are fixed

		// Equilibriums
		f_x[N+chemical_system.getReaction("OH-").getIndex()]
			= h * ho - chemical_system.getReaction("OH-").getK() * 1.0;
		f_x[N+chemical_system.getReaction("H3PO4").getIndex()]
			= h3po4 - chemical_system.getReaction("H3PO4").getK() * po4 * std::pow(h, 3);
		return f_x;
	}

	std::vector<std::vector<double>> analytical_df(
			const std::vector<std::size_t>& components_indexes,
			const std::vector<std::size_t>& species_indexes,
			const std::vector<double>& activities) const {
		std::vector<std::vector<double>> df_x(
				chemical_system.getComponents().size()
				+ chemical_system.getReactions().size());
		// Fixed component
		const double h = chemical_system.getSpecies("H+").activity();
		// Dynamic components
		const double po4 = activities[species_indexes[chemical_system.getSpecies("PO4-3").getIndex()]];

		// Count of fixed species (H+ and H2O)
		const std::size_t N = 2;
		// Count of mass conservation law equations
		const std::size_t M = chemical_system.getComponents().size()-N;
		// Mass conservation law
		{
				// PO4+
				auto& df = df_x[components_indexes[chemical_system.getComponent("PO4-3").getIndex()]];
				df.resize(chemical_system.getSpecies().size()-N);
				df[species_indexes[chemical_system.getSpecies("PO4-3").getIndex()]] = 1.0;
				df[species_indexes[chemical_system.getSpecies("H3PO4").getIndex()]] = 1.0;
				// df = 0.0 for other species
		}

		// df HO-
		{
			auto& df_h = df_x[M+chemical_system.getReaction("OH-").getIndex()];
			df_h.resize(chemical_system.getSpecies().size()-N);
			df_h[species_indexes[chemical_system.getSpecies("OH-").getIndex()]] = h;
		}

		// df H3PO4
		{
			double K = chemical_system.getReaction("H3PO4").getK();
			auto& df = df_x[M+chemical_system.getReaction("H3PO4").getIndex()];
			df.resize(chemical_system.getSpecies().size()-N);
			df[species_indexes[chemical_system.getSpecies("H3PO4").getIndex()]] = 1;
			df[species_indexes[chemical_system.getSpecies("PO4-3").getIndex()]] = -K * std::pow(h, 3);
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
			// Check NaOH
			double h3po4 = chemical_system.getSpecies("H3PO4").activity();
			double po4 = chemical_system.getSpecies("PO4-3").activity();
			double h = chemical_system.getSpecies("H+").activity();
			ASSERT_FLOAT_EQ(
					std::log10((h3po4)/(std::pow(h, 3) * po4)),
					chemical_system.getReaction("H3PO4").getLogK()
					);
		}
		{
			// Check pH
			double H = chemical_system.getSpecies("H+").activity();
			ASSERT_FLOAT_EQ(-std::log10(H), 7);
		}
	}
};

TEST_F(ComplexChemicalSystemTest, complex_H3PO4_reaction_f) {
	using namespace solver;
	F f(chemical_system);

	std::vector<double> activities(chemical_system.getSpecies().size()-2);
	{
		ChemicalSystem fake_system(chemical_system);
		// Example extents
		fake_system.proceed(chemical_system.getReaction("OH-"), -1e-7);
		fake_system.proceed(chemical_system.getReaction("H3PO4"), -1e-6);
		for(const auto& species : fake_system.getSpecies()) {
			if(species->getName() != "H2O" && species->getName() != "H+") {
				activities[f.speciesIndexes()[species->getIndex()]]
					= species->activity();
			}
		}
	}

	auto F = f.f(activities);
	auto analytical_F = analytical_f(
			f.componentsIndexes(),
			f.speciesIndexes(),
			activities);
	CHEM_LOGV(6) << "Computed F: " << F;
	CHEM_LOGV(6) << "Expected F: " << analytical_F;
	for (auto& component : chemical_system.getComponents()) {
		if(component->getSpecies()->getName() != "H2O"
				&& component->getSpecies()->getName() != "H+") {
			CHEM_LOGV(5) << "Checking mass conservation of component: "
				<< component->getSpecies()->getName();
			ASSERT_FLOAT_EQ(
					F[f.speciesIndexes()[component->getSpecies()->getIndex()]],
					analytical_F[f.speciesIndexes()[component->getSpecies()->getIndex()]]);
		}
	}
	for (auto& reaction : chemical_system.getReactions()) {
		CHEM_LOGV(5) << "Checking distance to equilibrium of reaction: "
			<< reaction->getName() << " (logK=" << reaction->getLogK() << ")";
		ASSERT_FLOAT_EQ(F[1+reaction->getIndex()], analytical_F[1+reaction->getIndex()]);
	}
}

TEST_F(ComplexChemicalSystemTest, complex_H3PO4_reaction_df) {
	using namespace solver;
	F f(chemical_system);

	std::vector<double> activities(chemical_system.getSpecies().size()-2);
	{
		ChemicalSystem fake_system(chemical_system);
		// Example extents
		fake_system.proceed(chemical_system.getReaction("OH-"), -1e-7);
		fake_system.proceed(chemical_system.getReaction("H3PO4"), -1e-6);
		for(const auto& species : fake_system.getSpecies()) {
			if(species->getName() != "H2O" && species->getName() != "H+") {
				activities[f.speciesIndexes()[species->getIndex()]]
					= species->activity();
			}
		}
	}

	for(const auto& species : chemical_system.getSpecies()) {
		if(species->getName() != "H2O" && species->getName() != "H+") {
			activities[f.speciesIndexes()[species->getIndex()]] = species->activity();
		}
	}

	auto dF = f.df(activities);
	auto analytical_dF = analytical_df(
			f.componentsIndexes(),
			f.speciesIndexes(),
			activities);

	CHEM_LOGV(5) << "Computed jacobian:";
	print_df(chemical_system, f, dF, 1);

	CHEM_LOGV(5) << "Expected jacobian:";
	print_df(chemical_system, f, analytical_dF, 1);

	for(const auto& component : chemical_system.getComponents()) {
		if(component->getSpecies()->getName() != "H2O"
				&& component->getSpecies()->getName() != "H+") {
			CHEM_LOGV(6) << "F: mass_conservation_law of "
				<< component->getSpecies()->getName();
			for(const auto& species : chemical_system.getSpecies()) {
				if(species->getName() != "H2O"
						&& species->getName() != "H+") {
					CHEM_LOGV(6) << "  dX: " << species->getName();
					ASSERT_FLOAT_EQ(
							dF
							[f.componentsIndexes()[component->getIndex()]]
							[f.speciesIndexes()[species->getIndex()]],
							analytical_dF
							[f.componentsIndexes()[component->getIndex()]]
							[f.speciesIndexes()[species->getIndex()]]
							);
				}
			}
		}
	}
	for (auto& reaction : chemical_system.getReactions()) {
		CHEM_LOGV(6) << "F: equilibrium of reaction "
			<< reaction->getName();
		for(const auto& species : chemical_system.getSpecies()) {
			if(species->getName() != "H2O" && species->getName() != "H+") {
				CHEM_LOGV(6) << "  dX: " << species->getName();
				ASSERT_FLOAT_EQ(
						dF
						[1 + reaction->getIndex()]
						[f.speciesIndexes()[species->getIndex()]],
						analytical_dF
						[1 + reaction->getIndex()]
						[f.speciesIndexes()[species->getIndex()]]);
			}
		}
	}
}

TEST_F(ComplexChemicalSystemTest, complex_H3PO4_reaction) {
	auto reaction_matrix = chemical_system.getReactionMatrix();

	chemical_system.setMaxIteration(10);
	chemical_system.solveEquilibrium();

	checkEquilibrium();

	CHEM_LOG(INFO) << "Complex H3PO4 reaction equilibrium test:";
	for(auto& species : chemical_system.getSpecies())
		CHEM_LOG(INFO) << "  " << species->getName() << ": " << species->concentration()/(1*mol/l) << " mol/l";
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


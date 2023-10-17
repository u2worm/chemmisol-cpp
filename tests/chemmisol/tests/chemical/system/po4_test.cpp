#include "gmock/gmock.h"
#include "chemmisol/tests/chemical/system/utils.h"

using namespace testing;
using namespace chemmisol;

/**
 * Tests ChemicalSystem equilibrium solving features with a reaction system
 * slightly more complex than NaClChemicalSystemTest, with non unitary
 * stoichiometric coefficients. However the only species with a non unitary
 * coefficient (3H+ in the H3PO4 reaction) is fixed so this coefficient only has
 * a real influence on F(X) but not on dF(X) or on the mass conservation law.
 */
class PO4ChemicalSystemTest : public Test {
	protected:
	ChemicalSystem chemical_system;
	solver::AbsoluteNewton absolute_newton;
	static constexpr double default_ph = 7;

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
			const std::vector<double>&) const {
		std::vector<std::vector<double>> df_x(
				chemical_system.getComponents().size()
				+ chemical_system.getReactions().size());
		// Fixed component
		const double h = chemical_system.getSpecies("H+").activity();

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
		// Mass conservation
		// Note: no mass conservation for H+ since H+ is fixed
		{
			double PO4 = chemical_system.getSpecies("PO4-3").quantity();
			double H3PO4 = chemical_system.getSpecies("H3PO4").quantity();
			ASSERT_FLOAT_EQ(
					PO4 + H3PO4,
					chemical_system.getComponent("PO4-3").getTotalQuantity()
					);
		}

		// Reaction quotients
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
			// Check H3PO4
			double h3po4 = chemical_system.getSpecies("H3PO4").activity();
			double po4 = chemical_system.getSpecies("PO4-3").activity();
			double h = chemical_system.getSpecies("H+").activity();
			ASSERT_FLOAT_EQ(
					std::log10((h3po4)/(std::pow(h, 3) * po4)),
					chemical_system.getReaction("H3PO4").getLogK()
					);
		}
	}

	void checkPh(double ph) {
		// Check pH
		double H = chemical_system.getSpecies("H+").activity();
		ASSERT_FLOAT_EQ(-std::log10(H), ph);
	}
};

TEST_F(PO4ChemicalSystemTest, complex_H3PO4_reaction_f) {
	using namespace solver;
	// Real domain F
	solver::ReducedChemicalSystem<solver::X> reduced_system(chemical_system);
	F<solver::X, solver::M> f(reduced_system, chemical_system);

	std::vector<double> activities(chemical_system.getSpecies().size()-2);
	{
		ChemicalSystem fake_system(chemical_system);
		// Example extents
		fake_system.proceed(fake_system.getReaction("OH-"), -1e-7);
		fake_system.proceed(fake_system.getReaction("H3PO4"), -1e-6);
		for(const auto& species : fake_system.getSpecies()) {
			if(species->getName() != "H2O" && species->getName() != "H+") {
				activities[reduced_system.speciesIndexes()[species->getIndex()]]
					= species->activity();
			}
		}
	}

	auto F = f.f(activities);
	auto analytical_F = analytical_f(
			reduced_system.componentsIndexes(),
			reduced_system.speciesIndexes(),
			activities);
	CHEM_LOGV(6) << "Computed F: " << F;
	CHEM_LOGV(6) << "Expected F: " << analytical_F;
	for (auto& component : chemical_system.getComponents()) {
		if(component->getSpecies()->getName() != "H2O"
				&& component->getSpecies()->getName() != "H+") {
			CHEM_LOGV(5) << "Checking mass conservation of component: "
				<< component->getSpecies()->getName();
			ASSERT_FLOAT_EQ(
					F[reduced_system.speciesIndexes()[component->getSpecies()->getIndex()]],
					analytical_F[reduced_system.speciesIndexes()[component->getSpecies()->getIndex()]]);
		}
	}
	for (auto& reaction : chemical_system.getReactions()) {
		CHEM_LOGV(5) << "Checking distance to equilibrium of reaction: "
			<< reaction->getName() << " (logK=" << reaction->getLogK() << ")";
		ASSERT_FLOAT_EQ(F[1+reaction->getIndex()], analytical_F[1+reaction->getIndex()]);
	}
}

TEST_F(PO4ChemicalSystemTest, complex_H3PO4_reaction_df) {
	using namespace solver;
	// Real domain F
	solver::ReducedChemicalSystem<solver::X> reduced_system(chemical_system);
	F<solver::X, solver::M> f(reduced_system, chemical_system);

	std::vector<double> activities(chemical_system.getSpecies().size()-2);
	{
		ChemicalSystem fake_system(chemical_system);
		// Example extents
		fake_system.proceed(fake_system.getReaction("OH-"), -1e-7);
		fake_system.proceed(fake_system.getReaction("H3PO4"), -1e-6);
		for(const auto& species : fake_system.getSpecies()) {
			if(species->getName() != "H2O" && species->getName() != "H+") {
				activities[reduced_system.speciesIndexes()[species->getIndex()]]
					= species->activity();
			}
		}
	}

	for(const auto& species : chemical_system.getSpecies()) {
		if(species->getName() != "H2O" && species->getName() != "H+") {
			activities[reduced_system.speciesIndexes()[species->getIndex()]] = species->activity();
		}
	}

	auto dF = f.df(activities);
	auto analytical_dF = analytical_df(
			reduced_system.componentsIndexes(),
			reduced_system.speciesIndexes(),
			activities);

	CHEM_LOGV(5) << "Computed jacobian:";
	print_df(chemical_system, reduced_system, dF, 1);

	CHEM_LOGV(5) << "Expected jacobian:";
	print_df(chemical_system, reduced_system, analytical_dF, 1);

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
							[reduced_system.componentsIndexes()[component->getIndex()]]
							[reduced_system.speciesIndexes()[species->getIndex()]],
							analytical_dF
							[reduced_system.componentsIndexes()[component->getIndex()]]
							[reduced_system.speciesIndexes()[species->getIndex()]]
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
						[reduced_system.speciesIndexes()[species->getIndex()]],
						analytical_dF
						[1 + reaction->getIndex()]
						[reduced_system.speciesIndexes()[species->getIndex()]]);
			}
		}
	}
}

TEST_F(PO4ChemicalSystemTest, solve_equilibrium_absolute_newton) {
	chemical_system.setMaxIteration(10);
	chemical_system.solveEquilibrium(absolute_newton);

	checkEquilibrium();
	checkPh(default_ph);
}

class PO4ChemicalSystemHomotopyTest :
	public PO4ChemicalSystemTest, public WithParamInterface<int> {
		protected:
			// seeds are defined in utils.h
			solver::HomotopyContinuation<std::minstd_rand> homotopy {
					std::minstd_rand {seeds[GetParam()]}, 100, 100
					};
};


TEST_P(PO4ChemicalSystemHomotopyTest, solve_equilibrium_homotopy) {
	chemical_system.solveEquilibrium(homotopy);
	checkEquilibrium();
	checkPh(default_ph);
}
INSTANTIATE_TEST_SUITE_P(HomotopyTest, PO4ChemicalSystemHomotopyTest, Range(0, 10));

TEST_F(PO4ChemicalSystemTest, solve_equilibrium_brutal_ph_homotopy) {
	solver::HomotopyContinuation<std::minstd_rand> homotopy(
			std::minstd_rand {}, 100, 100
			);

	chemical_system.fixPH(10);
	chemical_system.solveEquilibrium(homotopy);
	checkEquilibrium();
	checkPh(10);

	chemical_system.fixPH(2);
	chemical_system.solveEquilibrium(homotopy);
	checkEquilibrium();
	checkPh(2);
}

TEST_F(PO4ChemicalSystemTest, set_total_concentration) {
	chemical_system.setTotalConcentration(
			chemical_system.getComponent("PO4-3"),
			1.5e-6 * mol/l);

	chemical_system.solveEquilibrium(absolute_newton);

	ASSERT_FLOAT_EQ(
			chemical_system.getComponent("PO4-3").getTotalQuantity(),
			1.5e-6 * mol/l * AqueousSpecies::V);

	checkEquilibrium();
	checkPh(default_ph);
}


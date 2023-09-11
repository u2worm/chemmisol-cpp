#include "gmock/gmock.h"
#include "chemmisol/tests/chemical/system/utils.h"

using namespace testing;
using namespace chemmisol;

/**
 * Tests ChemicalSystem equilibrium solving features with a reaction system
 * slightly more complex than PO4ChemicalSystemTest, with non unitary
 * stoichiometric coefficients. This time, the species with a non unitary
 * coefficient (3Cl- in the AgCl3 reaction) is not fixed so it has an influence
 * on the mass conservation law, F(X) and dF(X).
 */
class AgClChemicalSystemTest : public Test {
	protected:
	ChemicalSystem chemical_system;

	void SetUp() override {
		chemical_system.addReaction("OH-", -13.997, {
				{"OH-", -1},
				{"H+", -1},
				{"H2O", 1}
				});

		chemical_system.addReaction("AgCl", 3.31, {
				{"AgCl", -1},
				{"Ag+", 1},
				{"Cl-", 1}
				});
		chemical_system.addReaction("AgCl2-", 5.25, {
				{"AgCl2-", -1},
				{"Ag+", 1},
				{"Cl-", 2}
				});
		chemical_system.addReaction("AgCl3-2", 5.2, {
				{"AgCl3-2", -1},
				{"Ag+", 1},
				{"Cl-", 3}
				});

		chemical_system.addComponent("Ag+", 0.1*mol/l);
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
		const double cl = activities[species_indexes[chemical_system.getSpecies("Cl-").getIndex()]];
		const double ag = activities[species_indexes[chemical_system.getSpecies("Ag+").getIndex()]];
		const double agcl = activities[species_indexes[chemical_system.getSpecies("AgCl").getIndex()]];
		const double agcl2 = activities[species_indexes[chemical_system.getSpecies("AgCl2-").getIndex()]];
		const double agcl3 = activities[species_indexes[chemical_system.getSpecies("AgCl3-2").getIndex()]];

		// Count of mass conservation law equations
		const std::size_t N = 2;

		std::vector<double> f_x(
				N
				// Count of reactions
				+ chemical_system.getReactions().size());
		// Mass conservation law
		f_x[components_indexes[chemical_system.getComponent("Ag+").getIndex()]]
			= (ag + agcl + agcl2 + agcl3) - chemical_system.getComponent("Ag+").getTotalQuantity();
		f_x[components_indexes[chemical_system.getComponent("Cl-").getIndex()]]
			= (cl + agcl + 2*agcl2 + 3*agcl3) - chemical_system.getComponent("Cl-").getTotalQuantity();
		// H+ and H2O are fixed

		// Equilibriums
		f_x[N+chemical_system.getReaction("OH-").getIndex()]
			= h * ho - chemical_system.getReaction("OH-").getK() * 1.0;
		f_x[N+chemical_system.getReaction("AgCl").getIndex()]
			= agcl - chemical_system.getReaction("AgCl").getK() * ag * cl;
		f_x[N+chemical_system.getReaction("AgCl2-").getIndex()]
			= agcl2 - chemical_system.getReaction("AgCl2-").getK() * ag * std::pow(cl, 2);
		f_x[N+chemical_system.getReaction("AgCl3-2").getIndex()]
			= agcl3 - chemical_system.getReaction("AgCl3-2").getK() * ag * std::pow(cl, 3);
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
		const double cl = activities[species_indexes[chemical_system.getSpecies("Cl-").getIndex()]];
		const double ag = activities[species_indexes[chemical_system.getSpecies("Ag+").getIndex()]];


		// Count of fixed species (H+ and H2O)
		const std::size_t N = 2;
		// Count of mass conservation law equations
		const std::size_t M = chemical_system.getComponents().size()-N;
		// Mass conservation law
		{
			{
				// Ag+
				auto& df = df_x[components_indexes[chemical_system.getComponent("Ag+").getIndex()]];
				df.resize(chemical_system.getSpecies().size()-N);
				df[species_indexes[chemical_system.getSpecies("Ag+").getIndex()]] = 1.0;
				df[species_indexes[chemical_system.getSpecies("AgCl").getIndex()]] = 1.0;
				df[species_indexes[chemical_system.getSpecies("AgCl2-").getIndex()]] = 1.0;
				df[species_indexes[chemical_system.getSpecies("AgCl3-2").getIndex()]] = 1.0;
				// df = 0.0 for other species
			}
			{
				// Cl-
				auto& df = df_x[components_indexes[chemical_system.getComponent("Cl-").getIndex()]];
				df.resize(chemical_system.getSpecies().size()-N);
				df[species_indexes[chemical_system.getSpecies("Cl-").getIndex()]] = 1.0;
				df[species_indexes[chemical_system.getSpecies("AgCl").getIndex()]] = 1.0;
				df[species_indexes[chemical_system.getSpecies("AgCl2-").getIndex()]] = 2.0;
				df[species_indexes[chemical_system.getSpecies("AgCl3-2").getIndex()]] = 3.0;
				// df = 0.0 for other species
			}
		}

		// df HO-
		{
			auto& df_h = df_x[M+chemical_system.getReaction("OH-").getIndex()];
			df_h.resize(chemical_system.getSpecies().size()-N);
			df_h[species_indexes[chemical_system.getSpecies("OH-").getIndex()]] = h;
		}

		// df AgCl
		{
			double K = chemical_system.getReaction("AgCl").getK();
			auto& df = df_x[M+chemical_system.getReaction("AgCl").getIndex()];
			df.resize(chemical_system.getSpecies().size()-N);
			df[species_indexes[chemical_system.getSpecies("AgCl").getIndex()]] = 1;
			df[species_indexes[chemical_system.getSpecies("Ag+").getIndex()]] = -K * cl;
			df[species_indexes[chemical_system.getSpecies("Cl-").getIndex()]] = -K * ag;
		}
		// df AgCl2
		{
			double K = chemical_system.getReaction("AgCl2-").getK();
			auto& df = df_x[M+chemical_system.getReaction("AgCl2-").getIndex()];
			df.resize(chemical_system.getSpecies().size()-N);
			df[species_indexes[chemical_system.getSpecies("AgCl2-").getIndex()]] = 1;
			df[species_indexes[chemical_system.getSpecies("Ag+").getIndex()]] = -K * std::pow(cl, 2);
			df[species_indexes[chemical_system.getSpecies("Cl-").getIndex()]] = -K * 2 * ag * cl;
		}
		// df AgCl3
		{
			double K = chemical_system.getReaction("AgCl3-2").getK();
			auto& df = df_x[M+chemical_system.getReaction("AgCl3-2").getIndex()];
			df.resize(chemical_system.getSpecies().size()-N);
			df[species_indexes[chemical_system.getSpecies("AgCl3-2").getIndex()]] = 1;
			df[species_indexes[chemical_system.getSpecies("Ag+").getIndex()]] = -K * std::pow(cl, 3);
			df[species_indexes[chemical_system.getSpecies("Cl-").getIndex()]] = -K * 3 * ag * std::pow(cl, 2);
		}
		return df_x;
	}

	void checkEquilibrium() {
		// Mass conservation
		// Note: no mass conservation for H+ since H+ is fixed
		{
			// Ag+
			double Ag = chemical_system.getSpecies("Ag+").quantity();
			double AgCl = chemical_system.getSpecies("AgCl").quantity();
			double AgCl2 = chemical_system.getSpecies("AgCl2-").quantity();
			double AgCl3 = chemical_system.getSpecies("AgCl3-2").quantity();
			ASSERT_FLOAT_EQ(
					Ag + AgCl + AgCl2 + AgCl3,
					chemical_system.getComponent("Ag+").getTotalQuantity()
					);
		}
		{
			// Cl-
			double Cl = chemical_system.getSpecies("Cl-").quantity();
			double AgCl = chemical_system.getSpecies("AgCl").quantity();
			double AgCl2 = chemical_system.getSpecies("AgCl2-").quantity();
			double AgCl3 = chemical_system.getSpecies("AgCl3-2").quantity();
			ASSERT_FLOAT_EQ(
					Cl + AgCl + 2*AgCl2 + 3*AgCl3,
					chemical_system.getComponent("Cl-").getTotalQuantity()
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
			// Check AgCl
			double agcl = chemical_system.getSpecies("AgCl").activity();
			double ag = chemical_system.getSpecies("Ag+").activity();
			double cl = chemical_system.getSpecies("Cl-").activity();
			ASSERT_FLOAT_EQ(
					std::log10(agcl/(ag * cl)),
					chemical_system.getReaction("AgCl").getLogK()
					);
		}
		{
			// Check AgCl2
			double agcl2 = chemical_system.getSpecies("AgCl2-").activity();
			double ag = chemical_system.getSpecies("Ag+").activity();
			double cl = chemical_system.getSpecies("Cl-").activity();
			ASSERT_FLOAT_EQ(
					std::log10(agcl2/(ag * std::pow(cl, 2))),
					chemical_system.getReaction("AgCl2-").getLogK()
					);
		}
		{
			// Check AgCl3
			double agcl3 = chemical_system.getSpecies("AgCl3-2").activity();
			double ag = chemical_system.getSpecies("Ag+").activity();
			double cl = chemical_system.getSpecies("Cl-").activity();
			ASSERT_FLOAT_EQ(
					std::log10(agcl3/(ag * std::pow(cl, 3))),
					chemical_system.getReaction("AgCl3-2").getLogK()
					);
		}
		{
			// Check pH
			double H = chemical_system.getSpecies("H+").activity();
			ASSERT_FLOAT_EQ(-std::log10(H), 7);
		}
	}
};

TEST_F(AgClChemicalSystemTest, complex_AgCl_reaction_f) {
	using namespace solver;
	F f(chemical_system);

	std::vector<double> activities(chemical_system.getSpecies().size()-2);
	{
		ChemicalSystem fake_system(chemical_system);
		// Example extents
		// TODO: add randomness
		fake_system.proceed(fake_system.getReaction("OH-"), -1e-7);
		fake_system.proceed(fake_system.getReaction("AgCl"), -1e-2);
		fake_system.proceed(fake_system.getReaction("AgCl2-"), -1e-2);
		fake_system.proceed(fake_system.getReaction("AgCl3-2"), -1e-2);
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

TEST_F(AgClChemicalSystemTest, complex_AgCl_reaction_df) {
	using namespace solver;
	F f(chemical_system);

	std::vector<double> activities(chemical_system.getSpecies().size()-2);
	{
		ChemicalSystem fake_system(chemical_system);
		// Example extents
		// TODO: add randomness
		fake_system.proceed(fake_system.getReaction("OH-"), -1e-7);
		fake_system.proceed(fake_system.getReaction("AgCl"), -1e-2);
		fake_system.proceed(fake_system.getReaction("AgCl2-"), -1e-2);
		fake_system.proceed(fake_system.getReaction("AgCl3-2"), -1e-2);
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

TEST_F(AgClChemicalSystemTest, solve_equilibrium) {
	chemical_system.setMaxIteration(10);
	chemical_system.solveEquilibrium();

	checkEquilibrium();
}

TEST_F(AgClChemicalSystemTest, set_total_concentration) {
	chemical_system.setTotalConcentration(
			chemical_system.getComponent("Ag+"),
			0.002 * mol/l);
	chemical_system.setTotalConcentration(
			chemical_system.getComponent("Cl-"),
			0.8 * mol/l);

	chemical_system.solveEquilibrium();

	ASSERT_FLOAT_EQ(
			chemical_system.getComponent("Ag+").getTotalQuantity(),
			0.002 * mol/l * AqueousSpecies::V);
	ASSERT_FLOAT_EQ(
			chemical_system.getComponent("Cl-").getTotalQuantity(),
			0.8 * mol/l * AqueousSpecies::V);

	checkEquilibrium();
}


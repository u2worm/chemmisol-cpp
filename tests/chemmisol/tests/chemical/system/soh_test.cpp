#include "gmock/gmock.h"
#include "chemmisol/tests/chemical/system/basic.h"
#include "chemmisol/tests/chemical/system/utils.h"

using namespace testing;
using namespace chemmisol;

/*
 * Tests ChemicalSystem equilibrium solving features with a basic reaction
 * system that involves mineral species.
 */
class SohTest : public BasicMineralChemicalSystemTest {
	protected:
		solver::AbsoluteNewton absolute_newton;

		std::vector<double> analytical_f(
				const std::vector<std::size_t>& components_indexes,
				const std::vector<std::size_t>& species_indexes,
				const std::vector<double>& activities) const {
			// Fixed component
			const double h = chemical_system.getSpecies("H+").activity();
			// Dynamic components
			const double oh = 
				activities[species_indexes[chemical_system.getSpecies("OH-").getIndex()]];
			const double n_soh = chemical_system.getSpecies("=SOH").quantity(
					activities[species_indexes[chemical_system.getSpecies("=SOH").getIndex()]]
					);
			const double soh =
				activities[species_indexes[chemical_system.getSpecies("=SOH").getIndex()]];
			const double n_soh2 = chemical_system.getSpecies("=SOH2").quantity(
					activities[species_indexes[chemical_system.getSpecies("=SOH2").getIndex()]]
					);
			const double soh2 =
				activities[species_indexes[chemical_system.getSpecies("=SOH2").getIndex()]];

			// Count of mass conservation law equations
			const std::size_t N = 1;

			std::vector<double> f_x(
					N
					// Count of reactions
					+ chemical_system.getReactions().size());
			// Mass conservation law
			f_x[components_indexes[chemical_system.getComponent("=SOH").getIndex()]]
				= (n_soh + n_soh2) - chemical_system.getComponent("=SOH").getTotalQuantity();
			// H+ and H2O are fixed

			// Equilibriums
			f_x[N+chemical_system.getReaction("OH-").getIndex()]
				= h * oh - chemical_system.getReaction("OH-").getK() * 1.0;
			f_x[N+chemical_system.getReaction("=SOH2").getIndex()]
				= soh2 - chemical_system.getReaction("=SOH2").getK() * soh * h;
			return f_x;
		}

		std::vector<std::vector<double>> analytical_df(
				const std::vector<std::size_t>& components_indexes,
				const std::vector<std::size_t>& species_indexes,
				const std::vector<double>& /*activities*/,
				double sites_quantity) const {
			std::vector<std::vector<double>> df_x(
					// Mass conservation of SOH
					1
					// Equilibrium of 2 reactions
					+ 2
					);
			// Fixed component
			const double h = chemical_system.getSpecies("H+").activity();


			// Count of fixed species (H+ and H2O)
			const std::size_t N = 2;
			// Count of mass conservation law equations
			const std::size_t M = 1;
			// Mass conservation law
			{
				{
					// =SOH
					auto& df = df_x[components_indexes[chemical_system.getComponent("=SOH").getIndex()]];
					df.resize(chemical_system.getSpecies().size()-N);
					df[species_indexes[chemical_system.getSpecies("=SOH").getIndex()]] = sites_quantity;
					df[species_indexes[chemical_system.getSpecies("=SOH2").getIndex()]] = sites_quantity;
					// df = 0.0 for other species
				}
			}

			// df HO-
			{
				auto& df_h = df_x[M+chemical_system.getReaction("OH-").getIndex()];
				df_h.resize(chemical_system.getSpecies().size()-N);
				df_h[species_indexes[chemical_system.getSpecies("OH-").getIndex()]] = h;
			}

			// df SOH2
			{
				double K = chemical_system.getReaction("=SOH2").getK();
				auto& df = df_x[M+chemical_system.getReaction("=SOH2").getIndex()];
				df.resize(chemical_system.getSpecies().size()-N);
				df[species_indexes[chemical_system.getSpecies("=SOH2").getIndex()]] = 1;
				df[species_indexes[chemical_system.getSpecies("=SOH").getIndex()]] = -K * h;
			}
			return df_x;
		}

		void checkEquilibrium() {
			// Check mass conservation
			// Note: no mass conservation for H+ since H+ is fixed
			{
				// =SOH
				double soh = chemical_system.getSpecies("=SOH").quantity();
				double soh2 = chemical_system.getSpecies("=SOH2").quantity();
				ASSERT_FLOAT_EQ(
						soh + soh2,
						chemical_system.getComponent("=SOH").getTotalQuantity()
						);
				// Bonus check: checks the sum of molar fractions is equal to 1
				ASSERT_FLOAT_EQ(
						chemical_system.getSpecies("=SOH").concentration() +
						chemical_system.getSpecies("=SOH2").concentration(),
						1.0
						);
			}

			// Reaction quotients
			{
				// Check OH-
				double h2o = chemical_system.getSpecies("H2O").activity();
				double h = chemical_system.getSpecies("H+").activity();
				double oh = chemical_system.getSpecies("OH-").activity();
				ASSERT_FLOAT_EQ(
						std::log10((h*oh)/h2o),
						chemical_system.getReaction("OH-").getLogK()
						);
			}
			{
				// Check =SOH2
				double soh = chemical_system.getSpecies("=SOH").activity();
				double soh2 = chemical_system.getSpecies("=SOH2").activity();
				double h = chemical_system.getSpecies("H+").activity();
				ASSERT_FLOAT_EQ(
						std::log10(soh2/(soh*h)),
						chemical_system.getReaction("=SOH2").getLogK()
						);
			}
			{
				// Check pH
				double H = chemical_system.getSpecies("H+").activity();
				ASSERT_FLOAT_EQ(-std::log10(H), 7);
			}
		}
};

TEST_F(SohTest, reaction_f) {
	using namespace solver;
	// Real domain F
	solver::ReducedChemicalSystem<solver::X> reduced_system(chemical_system);
	F<solver::X, solver::M> f(reduced_system, chemical_system);

	std::vector<double> activities(chemical_system.getSpecies().size()
			// Two fixed species, H2O and H+, that are removed from the
			// activities considered by the solver
			-2);
	{
		ChemicalSystem fake_system(chemical_system);
		// Example extents
		// TODO: add randomness
		fake_system.proceed(fake_system.getReaction("OH-"), -1e-7);
		fake_system.proceed(fake_system.getReaction("=SOH2"), -1e-6);
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

TEST_F(SohTest, reaction_df) {
	using namespace solver;
	// Real domain F
	solver::ReducedChemicalSystem<solver::X> reduced_system(chemical_system);
	F<solver::X, solver::M> f(reduced_system, chemical_system);

	std::vector<double> activities(chemical_system.getSpecies().size()
			// Two fixed species, H2O and H+, that are removed from the
			// activities considered by the solver
			-2);
	{
		ChemicalSystem fake_system(chemical_system);
		// Example extents
		// TODO: add randomness
		fake_system.proceed(fake_system.getReaction("OH-"), -1e-7);
		fake_system.proceed(fake_system.getReaction("=SOH2"), -1e-6);
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
			activities,
			chemical_system.sitesQuantity()
			);

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

TEST_F(SohTest, solve_equilibrium_absolute_newton) {
	chemical_system.setMaxIteration(10);
	chemical_system.solveEquilibrium(absolute_newton);

	checkEquilibrium();
}

TEST_F(SohTest, solve_equilibrium_homotopy) {
	solver::HomotopyContinuation<std::minstd_rand> homotopy(
			std::minstd_rand {}, 10, 10
			);
	chemical_system.solveEquilibrium(homotopy);

	checkEquilibrium();
}

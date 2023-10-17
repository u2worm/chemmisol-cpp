#include "chemmisol/tests/chemical/system/basic.h"
#include "chemmisol/tests/chemical/system/utils.h"

/*
 * Test equilibrium solving features for a system with unitary stoichiometric
 * coefficients only.
 */
class NaClChemicalSystemTest : public BasicAqueousChemicalSystemTest {
	protected:
		solver::AbsoluteNewton absolute_newton;
	public:
		static constexpr double default_ph = 7;
		void SetUp() override {
			BasicAqueousChemicalSystemTest::SetUp();

			chemical_system.fixPH(default_ph);
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
			// Check mass conservation
			// Note: no mass conservation for H+ since H+ is fixed
			{
				// Na
				double Na = chemical_system.getSpecies("Na+").quantity();
				double NaCl = chemical_system.getSpecies("NaCl").quantity();
				double NaOH = chemical_system.getSpecies("NaOH").quantity();
				ASSERT_FLOAT_EQ(
						Na + NaCl + NaOH,
						chemical_system.getComponent("Na+").getTotalQuantity()
						);
			}
			{
				// Cl
				double Cl = chemical_system.getSpecies("Cl-").quantity();
				double NaCl = chemical_system.getSpecies("NaCl").quantity();
				ASSERT_FLOAT_EQ(
						Cl + NaCl,
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
		}

		void checkPh(double ph) {
			// Check pH
			double H = chemical_system.getSpecies("H+").activity();
			ASSERT_FLOAT_EQ(-std::log10(H), ph);
		}
};

TEST_F(NaClChemicalSystemTest, NaCl_reaction_matrix) {
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

TEST_F(NaClChemicalSystemTest, NaCl_reaction_f) {
	using namespace solver;
	// Real domain F
	solver::ReducedChemicalSystem<solver::X> reduced_system(chemical_system);
	F<solver::X, solver::M> f(reduced_system, chemical_system);

	std::vector<double> activities(chemical_system.getSpecies().size()-2);
	{
		ChemicalSystem fake_system(chemical_system);
		// Example extents
		fake_system.proceed(fake_system.getReaction("OH-"), -1e-7);
		fake_system.proceed(fake_system.getReaction("NaCl"), -1e-3);
		fake_system.proceed(fake_system.getReaction("NaOH"), -1e-4);
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
	for (auto& component : chemical_system.getComponents()) {
		if(component->getSpecies()->getName() != "H2O"
				&& component->getSpecies()->getName() != "H+") {
			CHEM_LOGV(5) << "Checking mass conservation of component: "
				<< component->getSpecies()->getName();
			ASSERT_FLOAT_EQ(
					F[reduced_system.componentsIndexes()[component->getIndex()]],
					analytical_F[reduced_system.componentsIndexes()[component->getIndex()]]);
		}
	}
	for (auto& reaction : chemical_system.getReactions()) {
		CHEM_LOGV(5) << "Checking distance to equilibrium of reaction: "
			<< reaction->getName() << " (logK=" << reaction->getLogK() << ")";
		ASSERT_FLOAT_EQ(F[2+reaction->getIndex()], analytical_F[2+reaction->getIndex()]);
	}
}

TEST_F(NaClChemicalSystemTest, NaCl_reaction_df) {
	using namespace solver;
	// Real domain F
	solver::ReducedChemicalSystem<solver::X> reduced_system(chemical_system);
	F<solver::X, solver::M> f(reduced_system, chemical_system);

	std::vector<double> activities(chemical_system.getSpecies().size()-2);
	{
		ChemicalSystem fake_system(chemical_system);
		// Example extents
		fake_system.proceed(fake_system.getReaction("OH-"), -1e-7);
		fake_system.proceed(fake_system.getReaction("NaCl"), -1e-3);
		fake_system.proceed(fake_system.getReaction("NaOH"), -1e-4);
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
	print_df(chemical_system, reduced_system, dF, 2);

	CHEM_LOGV(5) << "Expected jacobian:";
	print_df(chemical_system, reduced_system, analytical_dF, 2);

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
						[2 + reaction->getIndex()]
						[reduced_system.speciesIndexes()[species->getIndex()]],
						analytical_dF
						[2 + reaction->getIndex()]
						[reduced_system.speciesIndexes()[species->getIndex()]]);
			}
		}
	}
}

TEST_F(NaClChemicalSystemTest, solve_equilibrium_absolute_newton) {
	chemical_system.setMaxIteration(10);
	chemical_system.solveEquilibrium(absolute_newton);

	checkEquilibrium();
	checkPh(default_ph);
}

TEST_F(NaClChemicalSystemTest, homotopy_G0) {
	solver::ReducedChemicalSystem<solver::CX> reduced_system(chemical_system);
	std::minstd_rand rd;
	solver::G g(reduced_system, rd);

	// Degrees of the F system:
	// - deg(F0) = 1 (Na mass conservation law)
	// - deg(F1) = 1 (Cl mass conservation law)
	// - def(F2) = 1 ([OH][H+] - [H2O] * K = 0 with [H+] and [H2O] fixed)
	// - def(F3) = 2 ([NaCl] - [Na][Cl] * K = 0)
	// - def(F4) = 1 ([NaOH][H+] - [Na][H2O] * K = 0 with [H+] and [H2O] fixed)
	// So the complete solution set of G(x) = 0 is of size 1*1*1*2*1=2.
	std::vector<solver::CX> G0 = g.initValues();
	ASSERT_THAT(G0, SizeIs(2));
	std::size_t i = 0;
	for(const auto& item : G0) {
		CHEM_LOGV(6) << "G0, " << i << " = " << item;
		++i;
		solver::CX gx = g.g(item);
		for(const auto& v : gx) {
			ASSERT_NEAR(v.imag(), 0.0, 1e-12);
			ASSERT_NEAR(v.real(), 0.0, 1e-12);
		}
	}
}

TEST_F(NaClChemicalSystemTest, solve_equilibrium_homotopy) {
	solver::HomotopyContinuation<std::minstd_rand> homotopy(
			std::minstd_rand {}, 1000, 100
			);

	// Run the solver several times to ensure the algorithm converges from any
	// random a/b selected to build G.
	for(std::size_t i = 0; i < 10; i++) {
		chemical_system.solveEquilibrium(homotopy);
		checkEquilibrium();
		checkPh(default_ph);
	}
}

TEST_F(NaClChemicalSystemTest, solve_equilibrium_brutal_ph_homotopy) {
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

TEST_F(NaClChemicalSystemTest, set_total_quantity) {
	chemical_system.setTotalQuantity(
			chemical_system.getComponent("Na+"),
			0.16*mol/l*AqueousSpecies::V);
	chemical_system.setTotalQuantity(
			chemical_system.getComponent("Cl-"),
			0.02*mol/l*AqueousSpecies::V);

	chemical_system.solveEquilibrium(absolute_newton);

	ASSERT_FLOAT_EQ(
			chemical_system.getComponent("Na+").getTotalQuantity(),
			0.16*mol/l*AqueousSpecies::V);
	ASSERT_FLOAT_EQ(
			chemical_system.getComponent("Cl-").getTotalQuantity(),
			0.02*mol/l*AqueousSpecies::V);

	checkEquilibrium();
}

TEST_F(NaClChemicalSystemTest, set_total_concentration) {
	chemical_system.setTotalConcentration(
			chemical_system.getComponent("Na+"),
			0.16*mol/l);
	chemical_system.setTotalConcentration(
			chemical_system.getComponent("Cl-"),
			0.02*mol/l);

	chemical_system.solveEquilibrium(absolute_newton);

	ASSERT_FLOAT_EQ(
			chemical_system.getComponent("Na+").getTotalQuantity(),
			0.16*mol/l*AqueousSpecies::V);
	ASSERT_FLOAT_EQ(
			chemical_system.getComponent("Cl-").getTotalQuantity(),
			0.02*mol/l*AqueousSpecies::V);

	checkEquilibrium();
}

TEST_F(NaClChemicalSystemTest, test_pH) {
	auto reaction_matrix = chemical_system.getReactionMatrix();
	chemical_system.setMaxIteration(200);
	chemical_system.solveEquilibrium(absolute_newton);

	// Check pH
	ASSERT_FLOAT_EQ(chemical_system.getPH(), 7);
}

TEST_F(NaClChemicalSystemTest, NaCl_reaction_quotient) {
	auto reaction_matrix = chemical_system.getReactionMatrix();
	chemical_system.setMaxIteration(200);
	chemical_system.solveEquilibrium(absolute_newton);

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


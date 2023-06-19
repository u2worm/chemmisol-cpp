#include "gmock/gmock.h"
#include <gtest/gtest.h>
#include <limits>
#include "chemical.h"

using namespace testing;
using namespace mineral;

class ChemicalSystemTest : public Test {
	protected:
	ChemicalSystem chemical_system;

	void SetUp() override {
		chemical_system.addReaction("OH-", 13.997, {
				{"OH-", -1},
				{"H+", -1},
				{"H2O", 1}
				});
		chemical_system.addReaction("NaCl", 0.3, {
				{"NaCl", -1},
				{"Na+", 1},
				{"Cl-", 1}
				});
		chemical_system.addReaction("NaOH", 13.897, {
				{"NaOH", -1},
				{"H+", -1},
				{"Na+", 1},
				{"H2O", 1}
				});

		chemical_system.addComponent("Na+", 0.1*mol/l);
		chemical_system.addComponent("Cl-", 0.1*mol/l);
		chemical_system.fixPH(7);

		chemical_system.initReactionMatrix();
	}

	std::array<double, 7> analytical_concentrations(const std::array<double, 4>& x) const {
		const double x_0 = x[chemical_system.getReaction("OH-").getId()];
		const double x_1 = x[chemical_system.getReaction("NaCl").getId()];
		const double x_2 = x[chemical_system.getReaction("NaOH").getId()];
		const double x_3 = x[chemical_system.getReaction("pH").getId()];

		double h = chemical_system.getComponent("H+").concentration();
		double ho = chemical_system.getComponent("OH-").concentration();
		double nacl = chemical_system.getComponent("NaCl").concentration();
		double na = chemical_system.getComponent("Na+").concentration();
		double cl = chemical_system.getComponent("Cl-").concentration();
		double naoh = chemical_system.getComponent("NaOH").concentration();
		const double V = Component::V;

		h += (-x_0-x_2-x_3)/V;
		ho += -x_0/V;
		nacl += -x_1/V;
		na += (x_1+x_2)/V;
		cl += x_1/V;
		naoh += -x_2/V;

		std::array<double, 7> concentrations;
		concentrations[chemical_system.getComponent("H+").getId()] = h;
		concentrations[chemical_system.getComponent("OH-").getId()] = ho;
		concentrations[chemical_system.getComponent("NaCl").getId()] = nacl;
		concentrations[chemical_system.getComponent("Na+").getId()] = na;
		concentrations[chemical_system.getComponent("Cl-").getId()] = cl;
		concentrations[chemical_system.getComponent("NaOH").getId()] = naoh;

		return concentrations;
	}
	std::array<double, 4> analytical_f(const std::array<double, 4>& x) const {
		auto concentrations = analytical_concentrations(x);
		
		const double h = concentrations[chemical_system.getComponent("H+").getId()];
		const double ho = concentrations[chemical_system.getComponent("OH-").getId()];
		const double nacl = concentrations[chemical_system.getComponent("NaCl").getId()];
		const double na = concentrations[chemical_system.getComponent("Na+").getId()];
		const double cl = concentrations[chemical_system.getComponent("Cl-").getId()];
		const double naoh = concentrations[chemical_system.getComponent("NaOH").getId()];

		std::array<double, 4> f_x;
		f_x[chemical_system.getReaction("OH-").getId()]
			= -std::log10(h/(1*mol/l)) - std::log10(ho/(1*mol/l)) - 13.997;
		f_x[chemical_system.getReaction("NaCl").getId()]
			= - std::log10(nacl) + std::log10(na) + std::log10(cl) - 0.3;
		f_x[chemical_system.getReaction("NaOH").getId()]
			= - std::log10(naoh) - std::log10(h) + std::log10(na) - 13.897;
		f_x[chemical_system.getReaction("pH").getId()]
			= -std::log10(h) - 7;
		return f_x;
	}

	std::array<std::array<double, 4>, 4> analytical_df(const std::array<double, 4>& x) const {
		std::array<std::array<double, 4>, 4> df_x;
		auto concentrations = analytical_concentrations(x);
		
		const double h = concentrations[chemical_system.getComponent("H+").getId()];
		const double ho = concentrations[chemical_system.getComponent("OH-").getId()];
		const double nacl = concentrations[chemical_system.getComponent("NaCl").getId()];
		const double na = concentrations[chemical_system.getComponent("Na+").getId()];
		const double cl = concentrations[chemical_system.getComponent("Cl-").getId()];
		const double naoh = concentrations[chemical_system.getComponent("NaOH").getId()];

		// df HO-
		{
			auto& df_h = df_x[chemical_system.getReaction("OH-").getId()];
			df_h[chemical_system.getReaction("OH-").getId()] = (1/h + 1/ho)*1/ln10;
			df_h[chemical_system.getReaction("NaCl").getId()] = 0;
			df_h[chemical_system.getReaction("NaOH").getId()] = 1/h * 1/ln10;
			df_h[chemical_system.getReaction("pH").getId()] = 1/h * 1/ln10;
		}

		// df NaCl
		{
			auto& df_nacl = df_x[chemical_system.getReaction("NaCl").getId()];
			df_nacl[chemical_system.getReaction("OH-").getId()] = 0;
			df_nacl[chemical_system.getReaction("NaCl").getId()] = (1/nacl + 1/na + 1/cl) * 1/ln10;
			df_nacl[chemical_system.getReaction("NaOH").getId()] = 1/na * 1/ln10;
			df_nacl[chemical_system.getReaction("pH").getId()] = 0;
		}

		// df NaOH
		{
			auto& df_naoh = df_x[chemical_system.getReaction("NaOH").getId()];
			df_naoh[chemical_system.getReaction("OH-").getId()] = 1/h * 1/ln10;
			df_naoh[chemical_system.getReaction("NaCl").getId()] = 1/na * 1/ln10;
			df_naoh[chemical_system.getReaction("NaOH").getId()] = (1/naoh + 1/na + 1/h) * 1/ln10;
			df_naoh[chemical_system.getReaction("pH").getId()] = 1/h * 1/ln10;
		}

		// df pH
		{
			auto& df_naoh = df_x[chemical_system.getReaction("pH").getId()];
			df_naoh[chemical_system.getReaction("OH-").getId()] = 1/h * 1/ln10;
			df_naoh[chemical_system.getReaction("NaCl").getId()] = 0;
			df_naoh[chemical_system.getReaction("NaOH").getId()] = 1/h * 1/ln10;
			df_naoh[chemical_system.getReaction("pH").getId()] = 1/h * 1/ln10;
		}
		return df_x;
	}
};

TEST_F(ChemicalSystemTest, activity) {
	ASSERT_FLOAT_EQ(chemical_system.getComponent("Na+").activity(), 0.1);
	ASSERT_FLOAT_EQ(chemical_system.getComponent("OH-").activity(), std::pow(10, -7));
	ASSERT_FLOAT_EQ(chemical_system.getComponent("H+").activity(), std::pow(10, -7));
}

TEST_F(ChemicalSystemTest, H2O_activity) {
	ASSERT_FLOAT_EQ(chemical_system.getComponent("H2O").activity(), 1.0);
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_matrix) {
	auto reaction_matrix = chemical_system.getReactionMatrix();

	ASSERT_THAT(reaction_matrix, SizeIs(3+1));
	{
		// Check OH-
		const std::vector<double>& reaction
			= reaction_matrix[chemical_system.getReaction("OH-").getId()];
		ASSERT_THAT(reaction[chemical_system.getComponent("OH-").getId()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getComponent("H+").getId()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getComponent("H2O").getId()], Eq(1));
	}
	{
		// Check NaCl
		const std::vector<double>& reaction
			= reaction_matrix[chemical_system.getReaction("NaCl").getId()];
		ASSERT_THAT(reaction[chemical_system.getComponent("NaCl").getId()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getComponent("Na+").getId()], Eq(1));
		ASSERT_THAT(reaction[chemical_system.getComponent("Cl-").getId()], Eq(1));
	}
	{
		// Check NaOH
		const std::vector<double>& reaction
			= reaction_matrix[chemical_system.getReaction("NaOH").getId()];
		ASSERT_THAT(reaction[chemical_system.getComponent("NaOH").getId()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getComponent("H+").getId()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getComponent("Na+").getId()], Eq(1));
		ASSERT_THAT(reaction[chemical_system.getComponent("H2O").getId()], Eq(1));
	}
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_concentrations) {
	using namespace solver;
	Solver::F f(chemical_system);

	std::array<double, 4> x;
	x[chemical_system.getReaction("OH-").getId()] = 0;
	x[chemical_system.getReaction("NaCl").getId()] = -1e-16;
	x[chemical_system.getReaction("NaOH").getId()] = -1e-16;
	x[chemical_system.getReaction("pH").getId()] = 0;
	auto C = f.concentrations({x.begin(), x.end()});
	auto analytical_C = analytical_concentrations(x);
	for (auto& reaction : chemical_system.getReactions()) {
		ASSERT_FLOAT_EQ(C[reaction->getId()], analytical_C[reaction->getId()]);
	}
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_f) {
	using namespace solver;
	Solver::F f(chemical_system);

	std::array<double, 4> x;
	x[chemical_system.getReaction("OH-").getId()] = 0;
	x[chemical_system.getReaction("NaCl").getId()] = -1e-16;
	x[chemical_system.getReaction("NaOH").getId()] = -1e-16;
	x[chemical_system.getReaction("pH").getId()] = 0;
	auto F = f.f({x.begin(), x.end()});
	auto analytical_F = analytical_f(x);
	for (auto& reaction : chemical_system.getReactions()) {
		ASSERT_FLOAT_EQ(F[reaction->getId()], analytical_F[reaction->getId()]);
	}
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_df) {
	using namespace solver;
	Solver::F f(chemical_system);

	std::array<double, 4> x;
	x[chemical_system.getReaction("OH-").getId()] = 0;
	x[chemical_system.getReaction("NaCl").getId()] = -1e-16;
	x[chemical_system.getReaction("NaOH").getId()] = -1e-16;
	x[chemical_system.getReaction("pH").getId()] = 0;
	auto dF = f.df({x.begin(), x.end()});
	auto analytical_dF = analytical_df(x);
	for (auto& reaction : chemical_system.getReactions()) {
		for (auto& extent : chemical_system.getReactions()) {
			ASSERT_FLOAT_EQ(
					dF[reaction->getId()][extent->getId()],
					analytical_dF[reaction->getId()][extent->getId()]);
		}
	}
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction) {
	auto reaction_matrix = chemical_system.getReactionMatrix();
	chemical_system.setMaxIteration(200);
	chemical_system.solveEquilibrium();

	{
		// Check OH-
		double H2O = chemical_system.getComponent("H2O").activity();
		double H = chemical_system.getComponent("H+").activity();
		double OH = chemical_system.getComponent("OH-").activity();
		ASSERT_FLOAT_EQ(
				std::log10(H2O/(H*OH)),
				chemical_system.getReaction("OH-").getLogK()
				);
	}
	{
		// Check NaCl
		double NaCl = chemical_system.getComponent("NaCl").activity();
		double Na = chemical_system.getComponent("Na+").activity();
		double Cl = chemical_system.getComponent("Cl-").activity();
		ASSERT_FLOAT_EQ(
				std::log10((Na*Cl)/NaCl),
				chemical_system.getReaction("NaCl").getLogK()
				);
	}
	{
		// Check NaOH
		double NaOH = chemical_system.getComponent("NaOH").activity();
		double Na = chemical_system.getComponent("Na+").activity();
		double H = chemical_system.getComponent("H+").activity();
		ASSERT_FLOAT_EQ(
				std::log10(Na/(NaOH*H)),
				chemical_system.getReaction("NaOH").getLogK()
				);
	}
	{
		// Check pH
		double H = chemical_system.getComponent("H+").activity();
		ASSERT_FLOAT_EQ(-std::log10(H), 7);
	}
	for(auto& component : chemical_system.getComponents())
		std::cout << component->getName() << ": " << component->concentration()/(1*mol/l) << " mol/l" << std::endl;
}

/*
 *TEST(ChemicalSystem, defaultSoil) {
 *    ChemicalSystem::defaultSoil();
 *}
 */

/*
 *class DerivativeTest : public testing::TestWithParam<int> {
 *};
 *
 *TEST_P(DerivativeTest, dx) {
 *    ChemicalSystem problem = ChemicalSystem::defaultEquilibrium();
 *    pHSolver::logPHSolver::F f(0, problem);
 *
 *    int dxn = GetParam();
 *    double delta = std::numeric_limits<double>::epsilon();
 *    for(int n = 0; n < 10; n++) {
 *        pHSolver::X x0 {{0, 0, 0, 0}};
 *        x0[dxn] = n*delta;
 *        pHSolver::X dx {{0, 0, 0, 0}};
 *        dx[dxn] = delta;
 *        auto x1 = x0 + dx;
 *        auto f_x0 = f.f(x0);
 *        auto f_x1 = f.f(x1);
 *
 *        auto df = f.df(x0);
 *        std::cout << "dx" << dxn+1 << std::endl;
 *        std::cout << "delta: " << delta << ", x0: " << x0[dxn] << ", x1:" << x1[dxn] << std::endl;
 *        for(std::size_t i = 0; i < 4; i++) {
 *            std::cout << "f" << i+1 << "(x1): " << f_x1[i] << ", f" << i+1 << "(x0): " << f_x0[i] << std::endl;
 *            std::cout << "f" << i+1 << ": " << (f_x1[i]-f_x0[i])/delta << " / "
 *                << df[i][dxn] << std::endl;
 *            if(df[i][dxn] != 0)
 *                ASSERT_NEAR((f_x1[i]-f_x0[i])/delta / df[i][dxn], 1, 1e-2);
 *            else
 *                ASSERT_NEAR((f_x1[i]-f_x0[i])/delta, 0, 1e-2);
 *
 *        }
 *        //ASSERT_NEAR((f_x1[0]-f_x0[0])/deltaX1 / f.df(x0)[0][0], 1, 1e-2);
 *        //ASSERT_NEAR((f_x1[1]-f_x0[1])/deltaX1 / f.df(x0)[1][0], 1, 1e-2);
 *        //ASSERT_NEAR((f_x1[2]-f_x0[2])/deltaX1 / f.df(x0)[2][0], 1, 1e-2);
 *        //ASSERT_NEAR((f_x1[3]-f_x0[3])/deltaX1 / f.df(x0)[3][0], 1, 1e-2);
 *    }
 *}
 *
 *INSTANTIATE_TEST_SUITE_P(Derivative,
 *                         DerivativeTest,
 *                         testing::Values(0, 1, 2, 3));
 */


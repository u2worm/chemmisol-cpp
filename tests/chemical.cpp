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
		chemical_system.addReaction("OH-", -13.997, {
				{"OH-", -1},
				{"H+", -1},
				{"H2O", 1}
				});
		chemical_system.addReaction("NaCl", -0.3, {
				{"NaCl", -1},
				{"Na+", 1},
				{"Cl-", 1}
				});
		chemical_system.addReaction("NaOH", -13.897, {
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
		const double x_0 = x[chemical_system.getReaction("OH-").getIndex()];
		const double x_1 = x[chemical_system.getReaction("NaCl").getIndex()];
		const double x_2 = x[chemical_system.getReaction("NaOH").getIndex()];
		const double x_3 = x[chemical_system.getReaction("pH").getIndex()];

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
		concentrations[chemical_system.getComponent("H+").getIndex()] = h;
		concentrations[chemical_system.getComponent("OH-").getIndex()] = ho;
		concentrations[chemical_system.getComponent("NaCl").getIndex()] = nacl;
		concentrations[chemical_system.getComponent("Na+").getIndex()] = na;
		concentrations[chemical_system.getComponent("Cl-").getIndex()] = cl;
		concentrations[chemical_system.getComponent("NaOH").getIndex()] = naoh;

		return concentrations;
	}
	std::array<double, 4> analytical_f(const std::array<double, 4>& x) const {
		auto concentrations = analytical_concentrations(x);
		
		const double h = concentrations[chemical_system.getComponent("H+").getIndex()];
		const double ho = concentrations[chemical_system.getComponent("OH-").getIndex()];
		const double nacl = concentrations[chemical_system.getComponent("NaCl").getIndex()];
		const double na = concentrations[chemical_system.getComponent("Na+").getIndex()];
		const double cl = concentrations[chemical_system.getComponent("Cl-").getIndex()];
		const double naoh = concentrations[chemical_system.getComponent("NaOH").getIndex()];

		std::array<double, 4> f_x;
		f_x[chemical_system.getReaction("OH-").getIndex()]
			= -std::log10(h/(1*mol/l)) - std::log10(ho/(1*mol/l)) - 13.997;
		f_x[chemical_system.getReaction("NaCl").getIndex()]
			= - std::log10(nacl) + std::log10(na) + std::log10(cl) - 0.3;
		f_x[chemical_system.getReaction("NaOH").getIndex()]
			= - std::log10(naoh) - std::log10(h) + std::log10(na) - 13.897;
		f_x[chemical_system.getReaction("pH").getIndex()]
			= -std::log10(h) - 7;
		return f_x;
	}

	std::array<std::array<double, 4>, 4> analytical_df(const std::array<double, 4>& x) const {
		std::array<std::array<double, 4>, 4> df_x;
		auto concentrations = analytical_concentrations(x);
		
		const double h = concentrations[chemical_system.getComponent("H+").getIndex()];
		const double ho = concentrations[chemical_system.getComponent("OH-").getIndex()];
		const double nacl = concentrations[chemical_system.getComponent("NaCl").getIndex()];
		const double na = concentrations[chemical_system.getComponent("Na+").getIndex()];
		const double cl = concentrations[chemical_system.getComponent("Cl-").getIndex()];
		const double naoh = concentrations[chemical_system.getComponent("NaOH").getIndex()];

		// df HO-
		{
			auto& df_h = df_x[chemical_system.getReaction("OH-").getIndex()];
			df_h[chemical_system.getReaction("OH-").getIndex()] = (1/h + 1/ho)*1/ln10;
			df_h[chemical_system.getReaction("NaCl").getIndex()] = 0;
			df_h[chemical_system.getReaction("NaOH").getIndex()] = 1/h * 1/ln10;
			df_h[chemical_system.getReaction("pH").getIndex()] = 1/h * 1/ln10;
		}

		// df NaCl
		{
			auto& df_nacl = df_x[chemical_system.getReaction("NaCl").getIndex()];
			df_nacl[chemical_system.getReaction("OH-").getIndex()] = 0;
			df_nacl[chemical_system.getReaction("NaCl").getIndex()] = (1/nacl + 1/na + 1/cl) * 1/ln10;
			df_nacl[chemical_system.getReaction("NaOH").getIndex()] = 1/na * 1/ln10;
			df_nacl[chemical_system.getReaction("pH").getIndex()] = 0;
		}

		// df NaOH
		{
			auto& df_naoh = df_x[chemical_system.getReaction("NaOH").getIndex()];
			df_naoh[chemical_system.getReaction("OH-").getIndex()] = 1/h * 1/ln10;
			df_naoh[chemical_system.getReaction("NaCl").getIndex()] = 1/na * 1/ln10;
			df_naoh[chemical_system.getReaction("NaOH").getIndex()] = (1/naoh + 1/na + 1/h) * 1/ln10;
			df_naoh[chemical_system.getReaction("pH").getIndex()] = 1/h * 1/ln10;
		}

		// df pH
		{
			auto& df_naoh = df_x[chemical_system.getReaction("pH").getIndex()];
			df_naoh[chemical_system.getReaction("OH-").getIndex()] = 1/h * 1/ln10;
			df_naoh[chemical_system.getReaction("NaCl").getIndex()] = 0;
			df_naoh[chemical_system.getReaction("NaOH").getIndex()] = 1/h * 1/ln10;
			df_naoh[chemical_system.getReaction("pH").getIndex()] = 1/h * 1/ln10;
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
			= reaction_matrix[chemical_system.getReaction("OH-").getIndex()];
		ASSERT_THAT(reaction[chemical_system.getComponent("OH-").getIndex()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getComponent("H+").getIndex()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getComponent("H2O").getIndex()], Eq(1));
	}
	{
		// Check NaCl
		const std::vector<double>& reaction
			= reaction_matrix[chemical_system.getReaction("NaCl").getIndex()];
		ASSERT_THAT(reaction[chemical_system.getComponent("NaCl").getIndex()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getComponent("Na+").getIndex()], Eq(1));
		ASSERT_THAT(reaction[chemical_system.getComponent("Cl-").getIndex()], Eq(1));
	}
	{
		// Check NaOH
		const std::vector<double>& reaction
			= reaction_matrix[chemical_system.getReaction("NaOH").getIndex()];
		ASSERT_THAT(reaction[chemical_system.getComponent("NaOH").getIndex()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getComponent("H+").getIndex()], Eq(-1));
		ASSERT_THAT(reaction[chemical_system.getComponent("Na+").getIndex()], Eq(1));
		ASSERT_THAT(reaction[chemical_system.getComponent("H2O").getIndex()], Eq(1));
	}
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_concentrations) {
	using namespace solver;
	F f(chemical_system);

	std::array<double, 4> x;
	x[chemical_system.getReaction("OH-").getIndex()] = 0;
	x[chemical_system.getReaction("NaCl").getIndex()] = -1e-16;
	x[chemical_system.getReaction("NaOH").getIndex()] = -1e-16;
	x[chemical_system.getReaction("pH").getIndex()] = 0;
	auto C = f.concentrations({x.begin(), x.end()});
	auto analytical_C = analytical_concentrations(x);
	for (auto& reaction : chemical_system.getReactions()) {
		ASSERT_FLOAT_EQ(C[reaction->getIndex()], analytical_C[reaction->getIndex()]);
	}
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_f) {
	using namespace solver;
	F f(chemical_system);

	std::array<double, 4> x;
	x[chemical_system.getReaction("OH-").getIndex()] = 0;
	x[chemical_system.getReaction("NaCl").getIndex()] = -1e-16;
	x[chemical_system.getReaction("NaOH").getIndex()] = -1e-16;
	x[chemical_system.getReaction("pH").getIndex()] = 0;
	auto F = f.f({x.begin(), x.end()});
	auto analytical_F = analytical_f(x);
	for (auto& reaction : chemical_system.getReactions()) {
		ASSERT_FLOAT_EQ(F[reaction->getIndex()], analytical_F[reaction->getIndex()]);
	}
}

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_df) {
	using namespace solver;
	F f(chemical_system);

	std::array<double, 4> x;
	x[chemical_system.getReaction("OH-").getIndex()] = 0;
	x[chemical_system.getReaction("NaCl").getIndex()] = -1e-16;
	x[chemical_system.getReaction("NaOH").getIndex()] = -1e-16;
	x[chemical_system.getReaction("pH").getIndex()] = 0;
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
	chemical_system.setMaxIteration(200);
	chemical_system.solveEquilibrium();

	{
		// Check OH-
		double H2O = chemical_system.getComponent("H2O").activity();
		double H = chemical_system.getComponent("H+").activity();
		double OH = chemical_system.getComponent("OH-").activity();
		ASSERT_FLOAT_EQ(
				std::log10((H*OH)/H2O),
				chemical_system.getReaction("OH-").getLogK()
				);
	}
	{
		// Check NaCl
		double NaCl = chemical_system.getComponent("NaCl").activity();
		double Na = chemical_system.getComponent("Na+").activity();
		double Cl = chemical_system.getComponent("Cl-").activity();
		ASSERT_FLOAT_EQ(
				std::log10(NaCl/(Na*Cl)),
				chemical_system.getReaction("NaCl").getLogK()
				);
	}
	{
		// Check NaOH
		double NaOH = chemical_system.getComponent("NaOH").activity();
		double Na = chemical_system.getComponent("Na+").activity();
		double H = chemical_system.getComponent("H+").activity();
		ASSERT_FLOAT_EQ(
				std::log10((NaOH*H)/Na),
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

TEST_F(ChemicalSystemTest, basic_NaCl_reaction_quotient) {
	auto reaction_matrix = chemical_system.getReactionMatrix();
	chemical_system.setMaxIteration(200);
	chemical_system.solveEquilibrium();

	{
		// Check OH-
		double H2O = chemical_system.getComponent("H2O").activity();
		double H = chemical_system.getComponent("H+").activity();
		double OH = chemical_system.getComponent("OH-").activity();
		ASSERT_FLOAT_EQ(
				chemical_system.reactionQuotient("OH-"),
				(H*OH)/H2O
				);
	}
	{
		// Check NaCl
		double NaCl = chemical_system.getComponent("NaCl").activity();
		double Na = chemical_system.getComponent("Na+").activity();
		double Cl = chemical_system.getComponent("Cl-").activity();
		ASSERT_FLOAT_EQ(
				chemical_system.reactionQuotient("NaCl"),
				NaCl/(Na*Cl)
				);
	}
	{
		// Check NaOH
		double NaOH = chemical_system.getComponent("NaOH").activity();
		double Na = chemical_system.getComponent("Na+").activity();
		double H = chemical_system.getComponent("H+").activity();
		ASSERT_FLOAT_EQ(
				chemical_system.reactionQuotient("NaOH"),
				(NaOH*H)/Na
				);
	}
	{
		// Check pH
		double H = chemical_system.getComponent("H+").activity();
		ASSERT_FLOAT_EQ(H, chemical_system.reactionQuotient("pH"));
	}
}


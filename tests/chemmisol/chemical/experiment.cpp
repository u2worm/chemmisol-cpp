#include "gmock/gmock.h"
#include "chemmisol/chemical/experiment.h"

using namespace testing;
using namespace chemmisol;

class TestExperiment : public Test {
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
				{"Na+", 1}, // AQUEOUS specified here only for test
									 // purpose, since it should be the default
									 // phase
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
	}
};


TEST_F(TestExperiment, test_ph_experiment) {
	unsigned long iterations = 100;
	std::vector<double> ph_values(iterations);
	for(std::size_t i = 0; i < iterations; i++) {
		ph_values[i] = 3.0 + i * (10.0-3.0)/iterations;
	}
	Experiment exp(
			"test.csv",
			chemical_system,
			[&ph_values] (unsigned long time) {
				return ph_values[time];
			}
			);

	auto results = exp.run(iterations);

	for(std::size_t i = 0; i < iterations; i++) {
		ASSERT_FLOAT_EQ(
				-std::log10(results["H+"][i]),
				ph_values[i]
				);
	}
}

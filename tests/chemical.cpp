#include "gmock/gmock.h"
#include <limits>
#include "chemical.h"

using namespace testing;
using namespace mineral;


/*
 *TEST(ChemicalSystem, defaultSoil) {
 *    ChemicalSystem::defaultSoil();
 *}
 */

class DerivativeTest : public testing::TestWithParam<int> {
};

TEST_P(DerivativeTest, dx) {
	ChemicalSystem problem = ChemicalSystem::defaultEquilibrium();
	pHSolver::F f(0, problem);

	int dxn = GetParam();
	double delta = std::numeric_limits<double>::epsilon();
	for(int n = 0; n < 10; n++) {
		pHSolver::X x0 {{0, 0, 0, 0}};
		x0[dxn] = n*delta;
		pHSolver::X dx {{0, 0, 0, 0}};
		dx[dxn] = delta;
		auto x1 = x0 + dx;
		auto f_x0 = f.f(x0);
		auto f_x1 = f.f(x1);

		auto df = f.df(x0);
		std::cout << "dx" << dxn+1 << std::endl;
		std::cout << "delta: " << delta << ", x0: " << x0[dxn] << ", x1:" << x1[dxn] << std::endl;
		for(std::size_t i = 0; i < 4; i++) {
			std::cout << "f" << i+1 << "(x1): " << f_x1[i] << ", f" << i+1 << "(x0): " << f_x0[i] << std::endl;
			std::cout << "f" << i+1 << ": " << (f_x1[i]-f_x0[i])/delta << " / "
				<< df[i][dxn] << std::endl;
			if(df[i][dxn] != 0)
				ASSERT_NEAR((f_x1[i]-f_x0[i])/delta / df[i][dxn], 1, 1e-2);
			else
				ASSERT_NEAR((f_x1[i]-f_x0[i])/delta, 0, 1e-2);

		}
		//ASSERT_NEAR((f_x1[0]-f_x0[0])/deltaX1 / f.df(x0)[0][0], 1, 1e-2);
		//ASSERT_NEAR((f_x1[1]-f_x0[1])/deltaX1 / f.df(x0)[1][0], 1, 1e-2);
		//ASSERT_NEAR((f_x1[2]-f_x0[2])/deltaX1 / f.df(x0)[2][0], 1, 1e-2);
		//ASSERT_NEAR((f_x1[3]-f_x0[3])/deltaX1 / f.df(x0)[3][0], 1, 1e-2);
	}
}

INSTANTIATE_TEST_SUITE_P(Derivative,
                         DerivativeTest,
                         testing::Values(0, 1, 2, 3));


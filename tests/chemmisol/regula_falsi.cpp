#include <gmock/gmock.h>
#include "chemmisol/regula_falsi.h"

TEST(RegulaFalsi, solve_iter) {
	double x = chemmisol::RegulaFalsi<double>(0, 10, [] (const double& x) {
			return std::pow(x, 3.7) - 7;
			}).solve_iter(10000);

	ASSERT_NEAR(x, std::pow(7, 1/3.7), 1e-3);
}

TEST(RegulaFalsi, solve_eps) {
	double x = chemmisol::RegulaFalsi<double>(0, 10, [] (const double& x) {
			return std::pow(x, 3.7) - 7;
			}).solve_eps(1e-3);

	ASSERT_NEAR(x, std::pow(7, 1/3.7), 1e-3);
}

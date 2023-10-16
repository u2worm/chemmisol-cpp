#include <gmock/gmock.h>
#include "chemmisol/math/regula_falsi.h"

TEST(RegulaFalsi, solve_iter) {
	double x = chemmisol::RegulaFalsi<double>().solve_iter(
			0, 10, [] (const double& x) {
			return std::pow(x, 3.7) - 7;
			},
			10000);

	ASSERT_NEAR(x, std::pow(7, 1/3.7), 1e-3);
}

TEST(RegulaFalsi, solve_eps) {
	double x = chemmisol::RegulaFalsi<double>().solve_eps(
			0, 10, [] (const double& x) {
			return std::pow(x, 3.7) - 7;
			},
			1e-3);

	ASSERT_NEAR(x, std::pow(7, 1/3.7), 1e-3);
}

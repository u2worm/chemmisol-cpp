#include "gtest/gtest.h"
#include "chemmisol/newton.h"

TEST(NewtonTest, square_3) {
	double x = chemmisol::Newton<double, double>(
			12.7f,
			[] (double x) -> double {return std::pow(x, 3) - 7;},
			[] (double x) -> double {return 3*std::pow(x, 2);}
		  ).solve_eps(1e-15);
	ASSERT_NEAR(x, 1.91f, 0.01f);
}

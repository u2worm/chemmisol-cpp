#include "gtest/gtest.h"
#include "chemmisol/math/newton.h"

TEST(NewtonTest, square_3) {
	auto result = chemmisol::Newton<double, double>().solve_eps(
			  12.7f,
			  [] (double x) -> double {return std::pow(x, 3) - 7;},
			  [] (double x) -> double {return 3*std::pow(x, 2);},
			  1e-15
			  );
	ASSERT_TRUE(result.isFinite());
	ASSERT_NEAR(result.x, 1.91f, 0.01f);
}

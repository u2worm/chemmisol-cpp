#include "gtest/gtest.h"
#include "chemmisol/math/linear.h"
#include "chemmisol/math/gauss.h"
#include <array>

using namespace testing;
using namespace chemmisol;

TEST(GaussTest, solve_array) {
	M<float, 3> m({{
			{{2, 1, -1}},
			{{-3, -1, 2}},
			{{-2, 1, 2}}
			}});
	X<float, 3> y = {{
		8,
		-11,
		-3
	}};

	auto x = gauss::solve(m, y);

	ASSERT_EQ(x[0], 2);
	ASSERT_EQ(x[1], 3);
	ASSERT_EQ(x[2], -1);
}

TEST(GaussTest, solve_vec) {
	VecM<float> m({
		{2, 1, -1},
		{-3, -1, 2},
		{-2, 1, 2}
	});
	VecX<float> y = {{
		8,
		-11,
		-3
	}};

	auto x = gauss::solve(m, y);

	ASSERT_EQ(x[0], 2);
	ASSERT_EQ(x[1], 3);
	ASSERT_EQ(x[2], -1);
}

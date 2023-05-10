#include "gtest/gtest.h"
#include "gauss.h"
#include "linear.h"
#include <array>

using namespace testing;
using namespace mineral;

TEST(GaussTest, solve) {
	M<float, 3> m({
		{2, 1, -1},
		{-3, -1, 2},
		{-2, 1, 2}
	});
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

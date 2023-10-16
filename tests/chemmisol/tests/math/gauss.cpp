#include "gmock/gmock.h"
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

	ASSERT_THAT(x, ElementsAre(2, 3, -1));
}

TEST(GaussTest, solve_array_view) {
	float x {};
	M<float, 4, 5> m({{
		{x, x, 2, 1, -1},
		{x, x, -3, -1, 2},
		{x, x, -2, 1, 2},
		{x, x, x, x, x}
	}});
	X<float, 4> y = {{
		x,
		8,
		-11,
		-3
	}};

	auto _x = gauss::solve(
			mview(m, 0, 3, 2, 5),
			xview(y, 1, 4)
			);

	ASSERT_THAT(_x, ElementsAre(_, 2, 3, -1));
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

	ASSERT_THAT(x, ElementsAre(2, 3, -1));
}

TEST(GaussTest, solve_vec_view) {
	float x {};
	VecM<float> m({
		{x, x, 2, 1, -1},
		{x, x, -3, -1, 2},
		{x, x, -2, 1, 2},
		{x, x, x, x, x}
	});
	VecX<float> y = {{
		x,
		8,
		-11,
		-3
	}};

	auto _x = gauss::solve(
			mview(m, 0, 3, 2, 5),
			xview(y, 1, 4)
			);

	ASSERT_THAT(_x, ElementsAre(_, 2, 3, -1));
}

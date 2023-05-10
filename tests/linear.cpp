#include "gmock/gmock.h"
#include "linear.h"

using namespace mineral;
using namespace testing;

TEST(LinearTest, X) {
	X<int, 4> x({3, 7, 12, 14});
	ASSERT_EQ(x[0], 3);
	ASSERT_EQ(x[1], 7);
	ASSERT_EQ(x[2], 12);
	ASSERT_EQ(x[3], 14);
}

TEST(LinearTest, M) {
	M<int, 4> x({
			{3, 7, 12, 14},
			{1, 4, 22, -6},
			{8, 0, 1, 16},
			{2, 9, 2, -5},
			});
	ASSERT_EQ(x[0][0], 3);
	ASSERT_EQ(x[0][1], 7);
	ASSERT_EQ(x[0][2], 12);
	ASSERT_EQ(x[0][3], 14);

	ASSERT_EQ(x[1][0], 1);
	ASSERT_EQ(x[1][1], 4);
	ASSERT_EQ(x[1][2], 22);
	ASSERT_EQ(x[1][3], -6);

	ASSERT_EQ(x[2][0], 8);
	ASSERT_EQ(x[2][1], 0);
	ASSERT_EQ(x[2][2], 1);
	ASSERT_EQ(x[2][3], 16);

	ASSERT_EQ(x[3][0], 2);
	ASSERT_EQ(x[3][1], 9);
	ASSERT_EQ(x[3][2], 2);
	ASSERT_EQ(x[3][3], -5);
}

TEST(LinearTest, product) {
	M<double, 4> m1({
		{1, 2, 3, 4},
		{5, 3, 7, 9},
		{3, 0, 2, 7},
		{4, 4, 1, 12}
	});

	X<double, 4> x = {{
		6,
		3,
		1.5,
		1
	}};

	X<double, 4> x1 = m1*x;
	ASSERT_THAT(x1, ElementsAre(
				Eq(1*6 + 2*3 + 3*1.5 + 4*1),
				Eq(5*6 + 3*3 + 7*1.5 + 1*9),
				Eq(3*6 + 0*3 + 2*1.5 + 7*1),
				Eq(4*6 + 4*3 + 1*1.5 + 12*1)
				));
}


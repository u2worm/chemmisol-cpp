#include "gmock/gmock.h"
#include "chemmisol/math/linear.h"

using namespace testing;
using namespace chemmisol;

TEST(LinearTest, X_array) {
	// Checks that the redefined constructor does the job
	X<int, 4> x({3, 7, 12, 14});
	ASSERT_THAT(x, ElementsAre(3, 7, 12, 14));
}

TEST(LinearTest, M_array) {
	// Checks that the redefined constructor does the job
	M<int, 4> x({{
			{{3, 7, 12, 14}},
			{{1, 4, 22, -6}},
			{{8, 0, 1, 16}},
			{{2, 9, 2, -5}},
			}});

	ASSERT_THAT(x, ElementsAre(
				ElementsAre(3, 7, 12, 14),
				ElementsAre(1, 4, 22, -6),
				ElementsAre(8, 0, 1, 16),
				ElementsAre(2, 9, 2, -5)
				));
}

TEST(LinearTest, minus_array) {
	// Checks that the redefined constructor does the job
	X<int, 4> x({3, 7, 12, 14});
	ASSERT_THAT(-x, ElementsAre(-3, -7, -12, -14));
}

TEST(LinearTest, sum_array) {
	// Checks that the redefined constructor does the job
	X<int, 4> x1({3, 7, 12, 14});
	X<int, 4> x2({8, 27, -11, 9});
	ASSERT_THAT(x1+x2, ElementsAre(3+8, 7+27, 12-11, 14+9));
}

TEST(LinearTest, product_array) {
	M<double, 4> m1({{
			{{1, 2, 3, 4}},
			{{5, 3, 7, 9}},
			{{3, 0, 2, 7}},
			{{4, 4, 1, 12}}
	}});

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

TEST(LinearTest, scalar_product_array) {
	M<double, 4> m({{
			{{1, 2, 3, 4}},
			{{5, 3, 7, 9}},
			{{3, 0, 2, 7}},
			{{4, 4, 1, 12}}
	}});

	X<double, 4> x = {{
		6,
		3,
		1.5,
		1
	}};

	X<double, 4> x1 = 7.2*x;
	ASSERT_THAT(x1, ElementsAre(
				Eq(7.2*x[0]),
				Eq(7.2*x[1]),
				Eq(7.2*x[2]),
				Eq(7.2*x[3])
				));

	M<double, 4> m1 = 4.3*m;
	ASSERT_THAT(m1, ElementsAre(
				ElementsAre(4.3*1, 4.3*2, 4.3*3, 4.3*4),
				ElementsAre(4.3*5, 4.3*3, 4.3*7, 4.3*9),
				ElementsAre(4.3*3, 4.3*0, 4.3*2, 4.3*7),
				ElementsAre(4.3*4, 4.3*4, 4.3*1, 4.3*12)
				));
}

TEST(LinearTest, augment_array) {
	M<double, 3> m1({{
			{{1, 2, 3}},
			{{5, 3, 7}},
			{{3, 0, 2}}
			}});

	X<double, 3> x = {{
		4,
		9,
		7
	}};

	auto m2 = augment(m1, x);
	ASSERT_THAT(m2, ElementsAre(
				ElementsAre(1, 2, 3, 4),
				ElementsAre(5, 3, 7, 9),
				ElementsAre(3, 0, 2, 7)
				));
}

TEST(LinearTest, minus_array_view) {
	// Checks that the redefined constructor does the job
	X<int, 4> x({3, 7, 12, 14});
	X<int, 4> min_x = -xview(x, 1, 3);
	ASSERT_THAT(-x, ElementsAre(_, -7, -12, _));
}

TEST(LinearTest, sum_array_view) {
	// Checks that the redefined constructor does the job
	X<int, 4> x1({3, 7, 12, 14});
	X<int, 4> x2({8, 27, -11, 9});
	X<int, 4> sum = xview(x1, 2, 4) + xview(x2, 2, 4);
	ASSERT_THAT(sum, ElementsAre(_, _, 12-11, 14+9));
}

TEST(LinearTest, product_array_view) {
	M<double, 4> m1({{
			{{1, 2, 3, 4}},
			{{5, 3, 7, 9}},
			{{3, 0, 2, 7}},
			{{4, 4, 1, 12}}
	}});

	X<double, 4> x = {{
		6,
		3,
		1.5,
		1
	}};

	X<double, 4> x1 = mview(m1,0, 0, 2, 3) * xview(x, 0, 3);
	ASSERT_THAT(x1, ElementsAre(
				Eq(1*6 + 2*3 + 3*1.5),
				Eq(5*6 + 3*3 + 7*1.5),
				Eq(3*6 + 0*3 + 2*1.5),
				_
				));
}

TEST(LinearTest, scalar_product_array_view) {
	M<double, 4> m({{
			{{1, 2, 3, 4}},
			{{5, 3, 7, 9}},
			{{3, 0, 2, 7}},
			{{4, 4, 1, 12}}
	}});

	X<double, 4> x = {{
		6,
		3,
		1.5,
		1
	}};

	X<double, 4> x1 = 7.2*xview(x, 0, 3);
	ASSERT_THAT(x1, ElementsAre(
				Eq(7.2*x[0]),
				Eq(7.2*x[1]),
				Eq(7.2*x[2]),
				_
				));

	M<double, 4> m1 = 4.3*mview(m, 2, 2, 4, 4);
	ASSERT_THAT(m1, ElementsAre(
				ElementsAre(_, _, _, _),
				ElementsAre(_, _, _, _),
				ElementsAre(_, _, 4.3*2, 4.3*7),
				ElementsAre(_, _, 4.3*1, 4.3*12)
				));
}

TEST(LinearTest, augment_array_view) {
	M<double, 3> m1({{
			{{1, 2, 3}},
			{{5, 3, 7}},
			{{3, 0, 2}}
			}});

	X<double, 3> x = {{
		4,
		9,
		7
	}};

	auto m2 = augment(mview(m1, 0, 1, 2, 3), xview(x, 0, 2));
	ASSERT_THAT(m2, ElementsAre(
				ElementsAre(_, 2, 3, 4),
				ElementsAre(_, 3, 7, 9),
				ElementsAre(_, _, _, _)
				));
}

TEST(LinearTest, minus_vec) {
	// Checks that the redefined constructor does the job
	VecX<int> x({3, 7, 12, 14});
	ASSERT_THAT(-x, ElementsAre(-3, -7, -12, -14));
}

TEST(LinearTest, sum_vec) {
	// Checks that the redefined constructor does the job
	VecX<int> x1({3, 7, 12, 14});
	VecX<int> x2({8, 27, -11, 9});
	ASSERT_THAT(x1+x2, ElementsAre(3+8, 7+27, 12-11, 14+9));
}

TEST(LinearTest, product_vec) {
	VecM<double> m1({
		{1, 2, 3, 4},
		{5, 3, 7, 9},
		{3, 0, 2, 7},
		{4, 4, 1, 12}
	});

	VecX<double> x = {{
		6,
		3,
		1.5,
		1
	}};

	VecX<double> x1 = m1*x;
	ASSERT_THAT(x1, ElementsAre(
				Eq(1*6 + 2*3 + 3*1.5 + 4*1),
				Eq(5*6 + 3*3 + 7*1.5 + 1*9),
				Eq(3*6 + 0*3 + 2*1.5 + 7*1),
				Eq(4*6 + 4*3 + 1*1.5 + 12*1)
				));
}

TEST(LinearTest, scalar_product_vec) {
	VecM<double> m({{
			{{1, 2, 3, 4}},
			{{5, 3, 7, 9}},
			{{3, 0, 2, 7}},
			{{4, 4, 1, 12}}
	}});

	VecX<double> x = {{
		6,
		3,
		1.5,
		1
	}};

	VecX<double> x1 = 7.2*x;
	ASSERT_THAT(x1, ElementsAre(
				Eq(7.2*x[0]),
				Eq(7.2*x[1]),
				Eq(7.2*x[2]),
				Eq(7.2*x[3])
				));

	VecM<double> m1 = 4.3*m;
	ASSERT_THAT(m1, ElementsAre(
				ElementsAre(4.3*1, 4.3*2, 4.3*3, 4.3*4),
				ElementsAre(4.3*5, 4.3*3, 4.3*7, 4.3*9),
				ElementsAre(4.3*3, 4.3*0, 4.3*2, 4.3*7),
				ElementsAre(4.3*4, 4.3*4, 4.3*1, 4.3*12)
				));
}
TEST(LinearTest, augment_vec) {
	VecM<double> m1({
		{1, 2, 3},
		{5, 3, 7},
		{3, 0, 2}
	});

	VecX<double> x = {{
		4,
		9,
		7
	}};

	auto m2 = augment(m1, x);
	ASSERT_THAT(m2, ElementsAre(
				ElementsAre(1, 2, 3, 4),
				ElementsAre(5, 3, 7, 9),
				ElementsAre(3, 0, 2, 7)
				));
}

TEST(LinearTest, minus_vec_view) {
	// Checks that the redefined constructor does the job
	VecX<int> x({3, 7, 12, 14});
	VecX<int> min_x = -xview(x, 1, 3);
	ASSERT_THAT(-x, ElementsAre(_, -7, -12, _));
}

TEST(LinearTest, sum_vec_view) {
	// Checks that the redefined constructor does the job
	VecX<int> x1({3, 7, 12, 14});
	VecX<int> x2({8, 27, -11, 9});
	VecX<int> sum = xview(x1, 2, 4) + xview(x2, 2, 4);
	ASSERT_THAT(sum, ElementsAre(_, _, 12-11, 14+9));
}

TEST(LinearTest, product_vec_view) {
	VecM<double> m1({{
			{{1, 2, 3, 4}},
			{{5, 3, 7, 9}},
			{{3, 0, 2, 7}},
			{{4, 4, 1, 12}}
	}});

	VecX<double> x = {{
		6,
		3,
		1.5,
		1
	}};

	VecX<double> x1 = mview(m1,0, 0, 2, 3) * xview(x, 0, 3);
	ASSERT_THAT(x1, ElementsAre(
				Eq(1*6 + 2*3 + 3*1.5),
				Eq(5*6 + 3*3 + 7*1.5),
				Eq(3*6 + 0*3 + 2*1.5),
				_
				));
}

TEST(LinearTest, scalar_product_vec_view) {
	VecM<double> m({{
			{{1, 2, 3, 4}},
			{{5, 3, 7, 9}},
			{{3, 0, 2, 7}},
			{{4, 4, 1, 12}}
	}});

	VecX<double> x = {{
		6,
		3,
		1.5,
		1
	}};

	VecX<double> x1 = 7.2*xview(x, 0, 3);
	ASSERT_THAT(x1, ElementsAre(
				Eq(7.2*x[0]),
				Eq(7.2*x[1]),
				Eq(7.2*x[2]),
				_
				));

	VecM<double> m1 = 4.3*mview(m, 2, 2, 4, 4);
	ASSERT_THAT(m1, ElementsAre(
				ElementsAre(_, _, _, _),
				ElementsAre(_, _, _, _),
				ElementsAre(_, _, 4.3*2, 4.3*7),
				ElementsAre(_, _, 4.3*1, 4.3*12)
				));
}

TEST(LinearTest, augment_vec_view) {
	VecM<double> m1({{
			{{1, 2, 3}},
			{{5, 3, 7}},
			{{3, 0, 2}}
			}});

	VecX<double> x = {{
		4,
		9,
		7
	}};

	auto m2 = augment(mview(m1, 0, 1, 2, 3), xview(x, 0, 2));
	ASSERT_THAT(m2, ElementsAre(
				ElementsAre(_, 2, 3, 4),
				ElementsAre(_, 3, 7, 9),
				ElementsAre(_, _, _, _)
				));
}


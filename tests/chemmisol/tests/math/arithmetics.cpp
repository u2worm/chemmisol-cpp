#include "gmock/gmock.h"
#include "chemmisol/math/arithmetics.h"

#include <random>

using namespace chemmisol;
using namespace testing;

TEST(Arithmetics, units_roots) {
	auto r = unit_roots<double>(14);

	ASSERT_THAT(r, SizeIs(14));
	auto sum = std::accumulate(r.begin(), r.end(), std::complex<double>(0, 0));
	ASSERT_NEAR(std::abs(sum), 0.0, 1e-14);

	for(auto _r : r) {
		auto pow = std::pow(_r, 14);
		ASSERT_NEAR(pow.real(), 1.0, 1e-14);
		ASSERT_NEAR(pow.imag(), 0.0, 1e-14);
	}
}

TEST(Arithmetics, complex_roots) {
	std::minstd_rand rg;
	std::uniform_real_distribution<double> rd(-10, 10);

	std::complex<double> c {rd(rg), rd(rg)};
	auto r = roots(c, 14);

	ASSERT_THAT(r, SizeIs(14));
	auto sum = std::accumulate(r.begin(), r.end(), std::complex<double>(0, 0));
	ASSERT_NEAR(std::abs(sum), 0.0, 1e-10);

	for(auto _r : r) {
		auto pow = std::pow(_r, 14);
		ASSERT_NEAR(pow.real(), c.real(), 1e-10);
		ASSERT_NEAR(pow.imag(), c.imag(), 1e-10);
	}
}

#include "gtest/gtest.h"
#include <complex>
#include "chemmisol/math/homotopy.h"
#include <random>
#include "chemmisol/math/arithmetics.h"

using namespace testing;
using namespace chemmisol;

std::minstd_rand rg;
std::uniform_real_distribution<double> rd {0.0, 1.0};

std::complex<double> random_complex() {
	return {rd(rg), rd(rg)};
}

class HomotopyTest : public Test {
	public:
		static const std::complex<double> a;
		static const std::complex<double> b;
		typedef std::complex<double> X;
		typedef std::complex<double> M;

		static X f(X x) {
			return 7.0*std::pow(x, 14) + 4.0*std::pow(x, 3) + 11.0 * x + 2.0;
		}

		static M df(X x) {
			return 14.0*7.0*std::pow(x, 13) + 3.0*4.0*std::pow(x, 2) + 11.0;
		}

		static X g(X x) {
			return a * std::pow(x, 14) - b;
		}

		static M dg(X x) {
			return 14.0*a*std::pow(x, 13);
		}

		static std::list<X> roots() {
			// Returns the roots of g
			return chemmisol::roots(b/a, 14);
		}
};
const std::complex<double> HomotopyTest::a = random_complex();
const std::complex<double> HomotopyTest::b = random_complex();

TEST_F(HomotopyTest, test) {
	Homotopy<X, M> homotopy(roots(), f, df, g, dg);

	std::list<std::complex<double>> x = homotopy.solve_iter(1000);

	for(auto v : x) {
		CHEM_LOGV(6) << "v: " << v << ", f(v): " << f(v);
		ASSERT_NEAR(std::abs(f(v)), 0.0, 1e-12);
	}
}

class Homotopy2Test : public Test {
	public:
		typedef std::array<std::complex<double>, 2> X;
		typedef std::array<std::array<std::complex<double>, 2>, 2> M;
		static const X a;
		static const X b;

		static X f(X x) {
			// Some random polynoms in two variables
			return {
				7.0*std::pow(x[0], 14)*std::pow(x[1], 3) + 4.0*std::pow(x[0], 3) + 0.2 * std::pow(x[1], 7) + 11.0 * std::pow(x[0], 2) * x[1] + 2.0,
				-18.0*std::pow(x[0], 2) + 28.0*std::pow(x[0], 3)*std::pow(x[1], 12) - 32.0 * x[0]  + 20.0 * x[1] + 4.0,
			};
		}

		static M df(X x) {
			return {{
				{{
					 // df0_dx0
					 7.0*13.0*std::pow(x[0], 13)*std::pow(x[1], 3) + 4.0*3.0*std::pow(x[0], 2) + 11.0*2.0*x[0]*x[1],
					 // df0_dx1
					 7.0*std::pow(x[0], 14)*3.0*std::pow(x[1], 2) + 7.0*0.2 * std::pow(x[1], 6) + 11.0 * std::pow(x[0], 2)
				}},
				{{
					 //df1_dx0
					 -18.0*2.0*x[0] + 28.0*3.0*std::pow(x[0], 2)*std::pow(x[1], 12) - 32.0,
					 //df1_dx1
					 28.0*std::pow(x[0], 3)*12.0*std::pow(x[1], 11) + 20.0
				}}
			}};
		}

		static X g(X x) {
			return {{
				a[0] * std::pow(x[0], 14+3) - b[0],
				a[1] * std::pow(x[1], 3+12) - b[1]
			}};
		}

		static M dg(X x) {
			return {{
				{{a[0] * (14.0+3.0) * std::pow(x[0], 14+3-1), 0}},
				{{0, a[1] * (3.0+12.0) * std::pow(x[1], 3+12-1)}}
			}};
		}

		static std::list<X> roots() {
			// Returns the roots of g
			auto g0_roots = chemmisol::roots(b[0]/a[0], 14+3);
			auto g1_roots = chemmisol::roots(b[1]/a[1], 3+12);
			auto it0 = g0_roots.begin();
			auto it1 = g1_roots.begin();
			std::list<X> roots;
			roots.push_back(X {*it0, *it1});
			it0++;
			roots.push_back(X {*it0, *it1});
			it0--;
			it1++;
			roots.push_back(X {*it0, *it1});
			it0++;
			it0++;
			roots.push_back(X {*it0, *it1});
			return roots;
		}
};

const Homotopy2Test::X Homotopy2Test::a {{ random_complex(), random_complex() }};
const Homotopy2Test::X Homotopy2Test::b {{ random_complex(), random_complex() }};

TEST_F(Homotopy2Test, roots) {
	auto r = roots();
	for(const auto& x : r) {
		auto gv = g(x);
		for(const auto& v : gv) {
			ASSERT_NEAR(v.real(), 0.0, 1e-12);
			ASSERT_NEAR(v.imag(), 0.0, 1e-12);
		}
	}
}

TEST_F(Homotopy2Test, test) {
	Homotopy<X, M> homotopy(roots(), f, df, g, dg);

	std::list<X> x = homotopy.solve_iter(1000);

	for(auto v : x) {
		CHEM_LOGV(6) << "v: " << v << ", f(v): " << f(v);
		X fv = f(v);
		for(const auto& v : fv) {
			ASSERT_NEAR(std::abs(v), 0.0, 1e-12);
		}
	}
}

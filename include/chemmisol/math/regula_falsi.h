#include <functional>
#include <stdexcept>
#include <iostream>
#include <math.h>
#include <string>

/**
 * @file chemmisol/math/regula_falsi.h
 *
 * Implementation of the [regula
 * falsi](https://en.wikipedia.org/wiki/Regula_falsi) method.
 */

namespace chemmisol {

	/**
	 * [Regula falsi](https://en.wikipedia.org/wiki/Regula_falsi) solver.
	 */
	template<typename T>
		class RegulaFalsi {
			private:
				T solve_iter(const std::function<T(const T&)>& f,
						T a, T b, T f_a, T f_b, std::size_t n) const;
				T solve_eps(const std::function<T(const T&)>& f,
						T a, T b, T f_a, T f_b, T eps) const;

			public:
				/**
				 * Runs the solver for n iterations and returns the found value
				 * of x such that `f(x)=0`.
				 *
				 * For the solver to work, f(a0) and f(b0) should be of opposite
				 * signs.
				 *
				 * @warning
				 * The algorithm is not guaranteed to converge.
				 *
				 * @param a0 Initial left bound of the bracketing interval.
				 * @param b0 Initial right bound of the bracketing interval.
				 * @param f Tunction to solve.
				 * @param n Count of iterations.
				 *
				 * @throws std::logic_error If f(a0) and f(b0) are not of
				 * opposite signs.
				 */
				T solve_iter(T a0, T b0, std::function<T(const T&)> f,
						std::size_t n) const;

				/**
				 * Runs the solver until `f(x) < eps` for n iterations and
				 * returns the found value of x such that `f(x) < eps`.
				 *
				 * For the solver to work, f(a0) and f(b0) should be of opposite
				 * signs.
				 *
				 * @warning
				 * The algorithm is not guaranteed to converge.
				 *
				 * @param a0 Initial left bound of the bracketing interval.
				 * @param b0 Initial right bound of the bracketing interval.
				 * @param f Tunction to solve.
				 * @param eps Required precision.
				 *
				 * @throws std::logic_error If f(a0) and f(b0) are not of
				 * opposite signs.
				 */
				T solve_eps(T a0, T b0, std::function<T(const T&)> f,
						T eps) const;
		};

	template<typename T>
		T RegulaFalsi<T>::solve_iter(
				T a0, T b0, std::function<T(const T&)> f,
				std::size_t n) const {
			T f_a = f(a0);
			T f_b = f(b0);
			if(f_a*f_b > 0)
				throw std::logic_error(
						" f(a_0) = " + std::to_string(f_a) +
						" and f(b_0) " + std::to_string(f_b) +
						" should be of opposite signs. (a_0=" +
						std::to_string(a0) + ", b_0=" + std::to_string(b0) + ")"
						);
			return solve_iter(f, a0, b0, f_a, f_b, n);
		}

	template<typename T>
		T RegulaFalsi<T>::solve_iter(
				const std::function<T(const T&)>& f,
				T a, T b, T f_a, T f_b, std::size_t n) const {
			if(f_a == f_b)
				return a;
			T c = (a * f_b - b * f_a) / (f_b - f_a);
			if(n==0)
				return c;
			T f_c = f(c);
			if(f_c == 0.0)
				return c;
			if(f_c * f_a < 0)
				return solve_iter(f, a, c, f_a, f_c, n-1);
			return solve_iter(f, c, b, f_c, f_b, n-1);
		}

	template<typename T>
		T RegulaFalsi<T>::solve_eps(
				T a0, T b0, std::function<T(const T&)> f,
				T eps) const {
			T f_a = f(a0);
			T f_b = f(b0);
			if(f_a*f_b > 0)
				throw std::logic_error(
						"f(a_0) = " + std::to_string(f_a) +
						" and f(b_0) " + std::to_string(f_b) +
						" should be of opposite signs. (a_0=" +
						std::to_string(a0) + ", b_0=" + std::to_string(b0) + ")"
						);
			return solve_eps(f, a0, b0, f_a, f_b, eps);
		}

	template<typename T>
		T RegulaFalsi<T>::solve_eps(
				const std::function<T(const T&)>& f,
				T a, T b, T f_a, T f_b, T eps) const {
			if(f_a == f_b)
				return a;
			T c = (a * f_b - b * f_a) / (f_b - f_a);
			T f_c = f(c);
			if(std::abs(f_c) < eps)
				return c;
			if(f_c * f_a < 0)
				return solve_eps(f, a, c, f_a, f_c, eps);
			return solve_eps(f, c, b, f_c, f_b, eps);
		}
}

#ifndef CHEMMISOL_NEWTON_H
#define CHEMMISOL_NEWTON_H

#include <functional>
#include <cmath>
#include "gauss.h"
#include "linear.h"

/**
 * @file chemmisol/math/newton.h
 *
 * Implementation of the [Newton-Raphson
 * method](https://en.wikipedia.org/wiki/Newton%27s_method).
 */

namespace chemmisol {
	
	/**
	 * Identity method, that returns x itself.
	 */
	template<typename T>
		T identity(const T& x) {
			return x;
		}

	/**
	 * Implementation of the [Newton-Raphson
	 * method](https://en.wikipedia.org/wiki/Newton%27s_method). The algorithm
	 * is generalized to solve systems of equations of [k variables and k
	 * functions](https://en.wikipedia.org/wiki/Newton%27s_method#k_variables,_k_functions).
	 *
	 * The algorithm recursively refines the x value using the following
	 * formulae:
	 *
	 *     x_n+1 = G(x_n - df(x_n)^-1 * f(x_n))
	 *
	 * where `f` denotes the function to solve so that `f(x)=0`, `df` denotes
	 * the [Jacobian
	 * matrix](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) of
	 * `f`, and `G` is a function that takes x as parameter and returns a vector
	 * of the same size.
	 *
	 * By default, G is the identity(), what allows to implement the classical
	 * Newton-Raphson method.
	 *
	 * When G is set to abs(), the so-called [AbsoluteNewton](@ref chemmisol::AbsoluteNewton) 
	 * method is implemented.
	 *
	 * The convergence is guaranteed only when the initial value is "close
	 * enough" to the desired solution.
	 *
	 * @tparam X Vector type of size n that defines the argument type and return
	 * value of the function `f`.
	 * @tparam M Matrix type of size n*n used to represent the Jacobian matrix
	 * of `f`, as returned by `df`.
	 * @tparam G Function applied to current x at the end of each iteration of
	 * the algorithm.
	 */
	template<typename X, typename M, X (&G)(const X&)=identity>
		class Newton {
			X x0;
			std::function<X(const X&)> f;
			std::function<M(const X&)> df;

			X solve_eps(const X& x, const X& f_x, float epsilon) const;
			X solve_iter(const X& x, const X& f_x, std::size_t n) const;

			public:
			/**
			 * Defines a Newton solver.
			 *
			 * @param x0 Initial x value, that should be "close enough" to the
			 * solution to guarantee convergence.
			 * @param f Function to solve so that `f(x)=0`.
			 * @param df Derivative of f, that returns the values of the
			 * Jacobian matrix of f for a specified x.
			 */
			Newton(
					const X& x0,
					std::function<X(const X&)> f,
					std::function<M(const X&)> df)
				: x0(x0), f(f), df(df) {
				}

			/**
			 * Runs the solver until `f(x) < eps` for n iterations and
			 * returns the found value of x such that `f(x) < eps`.
			 *
			 * @warning
			 * The algorithm is not guaranteed to converge if the initial x0 is
			 * not close enough to the solution.
			 *
			 * @param epsilon Required precision.
			 */
			X solve_eps(float epsilon) const;

			/**
			 * Runs the solver for n iterations and returns the found value
			 * of x such that `f(x)=0`.
			 *
			 * @warning
			 * The algorithm is not guaranteed to converge if the initial x0 is
			 * not close enough to the solution.
			 */
			X solve_iter(std::size_t n) const;
		};

	template<typename X, typename M, X (&G)(const X&)>
		X Newton<X, M, G>::solve_eps(float epsilon) const {
			return solve_eps(x0, f(x0), epsilon);
		}

	template<typename X, typename M, X (&G)(const X&)>
		X Newton<X, M, G>::solve_iter(std::size_t n) const {
			return solve_iter(x0, f(x0), n);
		}


	template<typename X, typename M, X (&G)(const X&)>
		X Newton<X, M, G>::solve_eps(const X& x, const X& f_x, float epsilon) const {
			CHEM_LOG(TRACE) << "[NEWTON] Current X: " << x;
			CHEM_LOG(TRACE) << "[NEWTON] Current F(X): " << f_x;
			CHEM_LOG(TRACE) << "[NEWTON] (e=" << epsilon << ") Epsilon: " << norm(f_x);
			if(norm(f_x) < epsilon)
				return x;
			X _x = gauss::solve(df(x), -f_x);
			X x1 = G(_x + x);
			return solve_eps(x1, f(x1), epsilon);
		}

	template<typename X, typename M, X (&G)(const X&)>
		X Newton<X, M, G>::solve_iter(const X& x, const X& f_x, std::size_t n) const {
			CHEM_LOG(TRACE) << "[NEWTON] Current X: " << x;
			CHEM_LOG(TRACE) << "[NEWTON] Current F(X): " << f_x;
			CHEM_LOG(TRACE) << "[NEWTON] (n=" << n << ") Epsilon: " << norm(f_x);
			if(n == 0)
				return x;
			X _x = gauss::solve(df(x), -f_x);
			X x1 = G(_x + x);
			return solve_iter(x1, f(x1), n-1);
		}

	/**
	 * Implements the absolute Newton-Raphson method, where the absolute value
	 * of `x_n+1` is taken as the result of each iteration.
	 *
	 * This variant is particularly suitable for chemical equilibrium solving,
	 * since it guarantees that x_n+1, that represents activity of species, is
	 * always positive.
	 *
	 * Moreover, the absolute Newton-Raphson method has proven to have a much
	 * larger [basin of
	 * attraction](https://en.wikipedia.org/wiki/Attractor#Basins_of_attraction)
	 * than the classical Newton-Raphson method [1].
	 *
	 * [1] K. Meintjes and A. P. Morgan, “A methodology for solving chemical
	 * equilibrium systems,” Applied Mathematics and Computation, vol. 22, no.
	 * 4, pp. 333–361, Jun. 1987, doi:
	 * [10.1016/0096-3003(87)90076-2](https://doi.org/10.1016/0096-3003(87)90076-2).

	 */
	template<typename X, typename M>
		using AbsoluteNewton = Newton<X, M, abs>;
}
#endif

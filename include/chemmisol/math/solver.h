#ifndef CHEMMISOL_MATH_SOLVER_H
#define CHEMMISOL_MATH_SOLVER_H

#include "linear.h"

namespace chemmisol {

	template<typename X>
		class SolverResult {
			public:
				X x;
				X f_x;
			private:
				bool is_finite;

			public:
				SolverResult() = default;
				SolverResult(const X& result, const X& f_x);

				bool isFinite() const {
					return is_finite;
				}
		};

	template<typename X>
		SolverResult<X>::SolverResult(const X& x, const X& f_x) :
			x(x), f_x(f_x), is_finite(norm(x)) {
			}

	template<typename X>
		MAKE_LOGGABLE(SolverResult<X>, result, os) {
			os << "(x=" << result.x << ", f(x)=" << result.f_x << ")";
			return os;
	}
}
#endif

#include <functional>
#include <cmath>
#include "gauss.h"
#include "../logging.h"

namespace chemmisol {
	template<typename T>
		double norm(const T& x) {
			return std::abs(x);
		}

	template<typename T>
		T abs(const T& x) {
			return std::abs(x);
		}

	template<typename T>
		T identity(const T& x) {
			return x;
		}

	template<typename X, typename M, X (&G)(const X&)=identity>
		class Newton {
			X x0;
			std::function<X(const X&)> f;
			std::function<M(const X&)> df;

			X solve_eps(const X& x, const X& f_x, float epsilon) const;
			X solve_iter(const X& x, const X& f_x, std::size_t n) const;

			public:
			Newton(
					const X& x0,
					std::function<X(const X&)> f,
					std::function<M(const X&)> df)
				: x0(x0), f(f), df(df) {
				}
			X solve_eps(float epsilon) const;
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

	template<typename X, typename M>
		using AbsoluteNewton = Newton<X, M, abs>;
}

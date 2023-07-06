#include <functional>
#include <cmath>
#include "gauss.h"

namespace chemmisol {
	template<typename T>
		double norm(const T& x) {
			return std::abs(x);
		}

	template<typename X, typename M>
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

	template<typename X, typename M>
		X Newton<X, M>::solve_eps(float epsilon) const {
			return solve_eps(x0, f(x0), epsilon);
		}

	template<typename X, typename M>
		X Newton<X, M>::solve_iter(std::size_t n) const {
			return solve_iter(x0, f(x0), n);
		}


	template<typename X, typename M>
		X Newton<X, M>::solve_eps(const X& x, const X& f_x, float epsilon) const {
			std::cout << "[NEWTON] Current X: " << x << std::endl;
			std::cout << "[NEWTON] Current F(X): " << f_x << std::endl;
			std::cout << "[NEWTON] Epsilon: " << norm(f_x) << std::endl;
			if(norm(f_x) < epsilon)
				return x;
			X _x = gauss::solve(df(x), -f_x);
			X x1 = _x + x;
			return solve_eps(x1, f(x1), epsilon);
		}

	template<typename X, typename M>
		X Newton<X, M>::solve_iter(const X& x, const X& f_x, std::size_t n) const {
			std::cout << "[NEWTON] Current X: " << x << std::endl;
			std::cout << "[NEWTON] Current F(X): " << f_x << std::endl;
			std::cout << "[NEWTON] Epsilon: " << norm(f_x) << std::endl;
			if(n == 0)
				return x;
			X _x = gauss::solve(df(x), -f_x);
			X x1 = _x + x;
			return solve_iter(x1, f(x1), n-1);
		}
}

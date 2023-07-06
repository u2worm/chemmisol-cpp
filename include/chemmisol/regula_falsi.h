#include <functional>
#include <stdexcept>
#include <iostream>
#include <math.h>
#include <string>

namespace chemmisol {

	template<typename F>
		class RegulaFalsi {
			private:
				F a_0;
				F b_0;

				std::function<F(const F&)> f;

				F solve_iter(F a, F b, F f_a, F f_b, std::size_t n) const;
				F solve_eps(F a, F b, F f_a, F f_b, F eps) const;

			public:
				RegulaFalsi(F a, F b, std::function<F(const F&)> f) :
					a_0(a), b_0(b), f(f) {
					}

				F solve_iter(std::size_t n) const;
				F solve_eps(F eps) const;
		};

	template<typename F>
		F RegulaFalsi<F>::solve_iter(std::size_t n) const {
			F f_a = f(a_0);
			F f_b = f(b_0);
			if(f_a*f_b > 0)
				throw std::logic_error(
						"f(a_0) = " + std::to_string(f_a) +
						" and f(b_0) " + std::to_string(f_b) +
						" should be of opposite signs. (a_0=" +
						std::to_string(a_0) + ", b_0=" + std::to_string(b_0) + ")"
						);
			return solve_iter(a_0, b_0, f_a, f_b, n);
		}

	template<typename F>
		F RegulaFalsi<F>::solve_iter(F a, F b, F f_a, F f_b, std::size_t n) const {
			//std::cout << "  a: " << a << std::endl;
			//std::cout << "  b: " << b << std::endl;
			//std::cout << "  fa: " << f_a << std::endl;
			//std::cout << "  fb: " << f_b << std::endl;
			if(f_a == f_b)
				return a;
			F c = (a * f_b - b * f_a) / (f_b - f_a);
			if(n==0)
				return c;
			F f_c = f(c);
			//std::cout << "  c: " << c << std::endl;
			//std::cout << "  fc: " << f_c << std::endl;
			if(f_c == 0.0)
				return c;
			if(f_c * f_a < 0)
				return solve_iter(a, c, f_a, f_c, n-1);
			return solve_iter(c, b, f_c, f_b, n-1);
		}

	template<typename F>
		F RegulaFalsi<F>::solve_eps(F eps) const {
			F f_a = f(a_0);
			F f_b = f(b_0);
			if(f_a*f_b > 0)
				throw std::logic_error(
						"f(a_0) = " + std::to_string(f_a) +
						" and f(b_0) " + std::to_string(f_b) +
						" should be of opposite signs. (a_0=" +
						std::to_string(a_0) + ", b_0=" + std::to_string(b_0) + ")"
						);
			return solve_eps(a_0, b_0, f_a, f_b, eps);
		}

	template<typename F>
		F RegulaFalsi<F>::solve_eps(F a, F b, F f_a, F f_b, F eps) const {
			if(f_a == f_b)
				return a;
			F c = (a * f_b - b * f_a) / (f_b - f_a);
			F f_c = f(c);
			if(std::abs(f_c) < eps)
				return c;
			if(f_c * f_a < 0)
				return solve_eps(a, c, f_a, f_c, eps);
			return solve_eps(c, b, f_c, f_b, eps);
		}
}

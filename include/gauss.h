#include <iostream>
#include <array>

namespace mineral { namespace gauss {
	template<typename X, typename A>
		X solve(const A& a);

	template<typename M, typename X>
		X solve(const M& m, const X& x) {
			return solve<X>(augment(m, x));
		}

	template<typename X, typename A>
		X solve(const A& m) {
			std::cout << "[GAUSS START]" << std::endl;
			A _m = m;
			std::size_t n = _m.size();
			std::cout << "n=" << n << std::endl;
			std::cout << "Step 0:" << std::endl;
			for(std::size_t i = 0; i < n; i++) {
				for(std::size_t j = 0; j < _m[i].size(); j++) {
					std::cout << _m[i][j] << ", ";
				}
				std::cout << std::endl;
			}

			for(std::size_t i = 0; i < n; i++) {
				for(std::size_t j = i+1; j < n; j++) {
					double ratio = _m[j][i]/_m[i][i];
					for(std::size_t k = 0; k < _m[i].size(); k++) {
						_m[j][k] = _m[j][k] - ratio*_m[i][k];
					}
				}
			}

			std::cout << "Step 1:" << std::endl;
			for(std::size_t i = 0; i < n; i++) {
				for(std::size_t j = 0; j < _m[i].size(); j++) {
					std::cout << _m[i][j] << ", ";
				}
				std::cout << std::endl;
			}
			X x;
			x[n-1] = _m[n-1][n]/_m[n-1][n-1];

			for(int i = n-2; i>=0; i--) {
				x[i] = _m[i][n];
				for(std::size_t j = i+1; j < n; j++) {
					x[i] = x[i] - _m[i][j]*x[j];
				}
				x[i] = x[i]/_m[i][i];
			}
			std::cout << "Step 2:" << std::endl;
			for(std::size_t i = 0; i < X::n; i++)
				std::cout << "x[" << i << "]=" << x[i] << std::endl;
			std::cout << "[GAUSS END]" << std::endl << std::endl;
			return x;
		}

	template<>
		double solve<double, double>(const double& f, const double& y);
}}

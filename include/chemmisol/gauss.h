#include <iostream>
#include <array>

namespace chemmisol { namespace gauss {
	template<typename X, typename A>
		X solve(const A& a);

	template<typename M, typename X>
		X solve(const M& m, const X& x) {
			std::cout << "[GAUSS START]" << std::endl;
			auto _m = augment(m, x);
			std::size_t n = _m.size();
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
			X _x = x;
			_x[n-1] = _m[n-1][n]/_m[n-1][n-1];

			for(int i = n-2; i>=0; i--) {
				_x[i] = _m[i][n];
				for(std::size_t j = i+1; j < n; j++) {
					_x[i] = _x[i] - _m[i][j]*_x[j];
				}
				_x[i] = _x[i]/_m[i][i];
			}
			std::cout << "Step 2:" << std::endl;
			for(std::size_t i = 0; i < _x.size(); i++)
				std::cout << "x[" << i << "]=" << _x[i] << std::endl;
			std::cout << "[GAUSS END]" << std::endl << std::endl;
			return _x;
		}

	template<>
		double solve<double, double>(const double& f, const double& y);
}}

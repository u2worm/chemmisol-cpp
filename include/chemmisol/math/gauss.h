#include <iostream>
#include <array>
#include "../logging.h"

namespace chemmisol { namespace gauss {
	template<typename X, typename A>
		X solve(const A& a);

	template<typename M, typename X>
		X solve(const M& m, const X& x) {
			CHEM_LOG(TRACE) << "[GAUSS START]";
			auto _m = augment(m, x);
			std::size_t n = _m.size();
			CHEM_LOG(TRACE) << "Step 0:";
			for(std::size_t i = 0; i < n; i++) {
				CHEM_LOG(TRACE) << _m[i];
			}

			for(std::size_t i = 0; i < n; i++) {
				for(std::size_t j = i+1; j < n; j++) {
					double ratio = _m[j][i]/_m[i][i];
					for(std::size_t k = 0; k < _m[i].size(); k++) {
						_m[j][k] = _m[j][k] - ratio*_m[i][k];
					}
				}
			}

			CHEM_LOG(TRACE) << "Step 1:";
			for(std::size_t i = 0; i < n; i++) {
				CHEM_LOG(TRACE) << _m[i];
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
			CHEM_LOG(TRACE) << "Step 2:";
			for(std::size_t i = 0; i < _x.size(); i++)
				CHEM_LOG(TRACE) << "x[" << i << "]=" << _x[i];
			CHEM_LOG(TRACE) << "[GAUSS END]";
			return _x;
		}

	template<>
		double solve<double, double>(const double& f, const double& y);
}}

#include <iostream>
#include <array>
#include "../logging.h"

/**
 * @file chemmisol/math/gauss.h
 *
 * Implementations of the [Gaussian elimination
 * algorithm](https://en.wikipedia.org/wiki/Gaussian_elimination).
 */

namespace chemmisol { namespace gauss {

	/**
	 * Solves m*x=y using the [Gaussian elimination
	 * algorithm](https://en.wikipedia.org/wiki/Gaussian_elimination).
	 *
	 * @param m Matrix of size n*n.
	 * @param y Matrix of size n.
	 * @return Vector x such that m*x=y.
	 *
	 * @tparam M Matrix type.
	 * @tparam X Vector type.
	 */
	template<typename M, typename X>
		X solve(const M& m, const X& y) {
			CHEM_LOG(TRACE) << "[GAUSS START]";
			auto _m = augment(m, y);
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
			
			X x = y;
			x[n-1] = _m[n-1][n]/_m[n-1][n-1];

			for(int i = n-2; i>=0; i--) {
				x[i] = _m[i][n];
				for(std::size_t j = i+1; j < n; j++) {
					x[i] = x[i] - _m[i][j]*x[j];
				}
				x[i] = x[i]/_m[i][i];
			}
			CHEM_LOG(TRACE) << "Step 2:";
			for(std::size_t i = 0; i < x.size(); i++)
				CHEM_LOG(TRACE) << "x[" << i << "]=" << x[i];
			CHEM_LOG(TRACE) << "[GAUSS END]";
			return x;
		}

	/**
	 * Specialisation of the Gaussian elimination for scalar value.
	 *
	 * The Gaussiant elimination is useless in this case, but this
	 * specialisation is necessary for the implementation of the Newton
	 * solver for scalar values.
	 *
	 * @param f Coefficient.
	 * @param y Value.
	 * @return Value x such that f*x=y.
	 */
	template<>
		double solve<double, double>(const double& f, const double& y);
}}

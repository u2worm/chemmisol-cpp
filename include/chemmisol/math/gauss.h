#ifndef CHEMMISOL_GAUSS_H
#define CHEMMISOL_GAUSS_H

#include <iostream>
#include <array>
#include <complex>
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
				CHEM_LOG(TRACE) << i << ": " << _m[i];
			}

			for(std::size_t i = 0; i < n; i++) {
				//CHEM_LOGV(9) << "Step 0." << i << "/" << n << ", current echelon:";
				//for(std::size_t k = 0; k < n; k++)
					//CHEM_LOGV(9) << k << ": " << _m[k];
				for(std::size_t j = i+1; j < n; j++) {
					auto m_j_i = _m[j][i];
					auto m_i_i = _m[i][i];
					//CHEM_LOGV(9) << "m[" << j << "] = m[" << j << "]-(" << _m[j][i] << "/" << _m[i][i] << ")*m[" << i << "]";
					// By definition, _m[j][i]-_m[i][i]/_m[i][i]*_m[j][i]=1.0
					_m[j][i] = 1.0;
					// Do not need to change coefficients from 0 to i+1.
					// They should all be 0 to build the echelon, but not used
					// in solving.
					for(std::size_t k = i+1; k < _m[i].size(); k++) {
						_m[j][k] = _m[j][k] - _m[i][k]/m_i_i*m_j_i;
					}
					//CHEM_LOGV(9) << j << ": " << _m[j];
				}
			}

			CHEM_LOG(TRACE) << "Step 1:";
			for(std::size_t i = 0; i < n; i++) {
				CHEM_LOG(TRACE) << i << ": " << _m[i];
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

	template<>
		std::complex<double> solve<std::complex<double>, std::complex<double>>(
				const std::complex<double>& f, const std::complex<double>& y);

}}
#endif

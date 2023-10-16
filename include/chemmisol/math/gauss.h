#ifndef CHEMMISOL_GAUSS_H
#define CHEMMISOL_GAUSS_H

#include <iostream>
#include <array>
#include <complex>
#include "linear.h"
#include "../logging.h"

/**
 * @file chemmisol/math/gauss.h
 *
 * Implementations of the [Gaussian elimination
 * algorithm](https://en.wikipedia.org/wiki/Gaussian_elimination).
 */

namespace chemmisol { namespace gauss {

	template<typename M>
		struct AView {
			MView<M> m_view;

			AView(const M& m,
					std::size_t a0, std::size_t a1, std::size_t b0, std::size_t b1)
				: m_view(mview(m, a0, a1, b0, b1)) {
				}
		};
	template<typename M>
		AView<M> aview(const M& m,
					std::size_t a0, std::size_t a1, std::size_t b0, std::size_t b1) {
			return {m, a0, a1, b0, b1};
		}

	template<typename M>
		inline MAKE_LOGGABLE(AView<M>, a_view, os) {
			for(std::size_t i = a_view.m_view.a0; i < a_view.m_view.a1; i++) {
				os << std::endl << std::setw((int) std::log10(a_view.m_view.a1)+1)
					<< i << ": " <<
					xview(a_view.m_view.m[i], a_view.m_view.b0, a_view.m_view.b1)
					<< xview(a_view.m_view.m[i], a_view.m_view.m[i].size()-1, a_view.m_view.m[i].size());
			}
			return os;
		}


	template<typename M, typename X>
	struct Gauss {
		static X solve(const M& m, const X& y);
		static X solve(const MView<M>& m, const XView<X>& y);
	};

	
	template<typename M, typename X>
		X Gauss<M, X>::solve(const M& m, const X& y) {
			return solve(mview(m), xview(y));
		}

	template<typename M, typename X>
		X Gauss<M, X>::solve(const MView<M>& m_view, const XView<X>& y_view) {
			CHEM_LOG(TRACE) << "[GAUSS START]";
			auto _m = augment(m_view, y_view);
			CHEM_LOG(TRACE) << "Step 0:"
				<< aview(
						_m, m_view.a0, m_view.a1, m_view.b0, m_view.b1
						);
		
			for(std::size_t i = m_view.a0; i < m_view.a1; i++) {
				CHEM_LOGV(9) << "Step 0." << i << "/" << m_view.b0 - m_view.a0 << ", current echelon:"
					<< aview(_m, m_view.a0, m_view.a1, m_view.b0, m_view.b1);
				for(std::size_t j = i+1; j < m_view.a1; j++) {
					auto m_j_i = _m[j][i];
					auto m_i_i = _m[i][i];
					CHEM_LOGV(9) << "m[" << j << "] = m[" << j << "]-(" << _m[j][i] << "/" << _m[i][i] << ")*m[" << i << "]";
					// By definition, _m[j][i]-_m[i][i]/_m[i][i]*_m[j][i]=1.0
					_m[j][i] = 1.0;
					// Do not need to change coefficients from 0 to i+1.
					// They should all be 0 to build the echelon, but not used
					// in solving.
					for(std::size_t k = i+1; k < _m[i].size(); k++) {
						_m[j][k] = _m[j][k] - _m[i][k]/m_i_i*m_j_i;
					}
					CHEM_LOGV(9) << j << ": " << _m[j];
				}
			}

			CHEM_LOG(TRACE) << "Step 1:"
				<< aview(_m, m_view.a0, m_view.a1, m_view.b0, m_view.b1);
			
			X x = y_view.x;
			x[y_view.a1-1] = _m[m_view.a1-1][m_view.b1]/_m[m_view.a1-1][m_view.b1-1];

			for(int i = y_view.a1-2; i>=0; i--) {
				x[i] = _m[i][m_view.b1];
				for(std::size_t j = i+1; j < m_view.b1; j++) {
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
			return Gauss<M, X>::solve(m, y);
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
		struct Gauss<double, double> {
			static double solve(const double& f, const double& y);
		};

	template<>
		struct Gauss<std::complex<double>, std::complex<double>> {
			static std::complex<double> solve(
					const std::complex<double>& f, const std::complex<double>& y);
		};
}}
#endif

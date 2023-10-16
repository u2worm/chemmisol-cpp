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

	template<typename A>
		struct AView {
			const A& a;
			std::size_t ma0;
			std::size_t ma1;
			std::size_t mb0;
			std::size_t mb1;
			std::size_t ya0;
			std::size_t ya1;

			AView(const A& a,
					std::size_t ma0, std::size_t ma1, std::size_t mb0, std::size_t mb1,
					std::size_t ya0, std::size_t ya1)
				: a(a), ma0(ma0), ma1(ma1), mb0(mb0), mb1(mb1), ya0(ya0), ya1(ya1) {
				}
		};

	template<typename A>
		AView<A> aview(const A& a,
					std::size_t ma0, std::size_t ma1, std::size_t mb0, std::size_t mb1,
					std::size_t ya0, std::size_t ya1) {
			return {a, ma0, ma1, mb0, mb1, ya0, ya1};
		}

	template<typename A>
		inline MAKE_LOGGABLE(AView<A>, a_view, os) {
			for(std::size_t i = a_view.ma0; i < a_view.ma1; i++) {
				os << std::endl << std::setw((int) std::log10(a_view.ma1)+1)
					<< i << ": " <<
					xview(a_view.a[i], a_view.mb0, a_view.mb1)
					<< xview(
							a_view.a[a_view.ya0-a_view.ma0+i],
							a_view.a[a_view.ya0-a_view.ma0+i].size()-1,
							a_view.a[a_view.ya0-a_view.ma0+i].size());
			}
			return os;
		}


	template<typename M, typename X>
	struct Gauss {
		typedef X result_type;

		static X solve(const M& m, const X& y);
	};

	template<typename M, typename X>
		struct Gauss<MView<M>, XView<X>> {
			typedef X result_type;

			static X solve(const MView<M>& m, const XView<X>& y);
		};

	
	template<typename M, typename X>
		X Gauss<M, X>::solve(const M& m, const X& y) {
			return Gauss<MView<M>, XView<X>>::solve(mview(m), xview(y));
		}

	template<typename M, typename X>
		X Gauss<MView<M>, XView<X>>::solve(const MView<M>& m_view, const XView<X>& y_view) {
			CHEM_LOG(TRACE) << "[GAUSS START]";
			auto _m = augment(m_view, y_view);
			CHEM_LOG(TRACE) << "Step 0:"
				<< aview(
						_m, m_view.a0, m_view.a1, m_view.b0, m_view.b1,
						y_view.a0, y_view.a1
						);
		
			for(std::size_t i = 0; i < m_view.a1-m_view.a0; i++) {
				CHEM_LOGV(9) << "Step 0." << i << "/" << m_view.a1 - m_view.a0 << ", current echelon:"
					<< aview(_m, m_view.a0, m_view.a1, m_view.b0, m_view.b1,
							y_view.a0, y_view.a1);
				for(std::size_t j = i+1; j < m_view.a1-m_view.a0; j++) {
					auto m_j_i = _m[m_view.a0+j][m_view.b0+i];
					auto m_i_i = _m[m_view.a0+i][m_view.b0+i];
					CHEM_LOGV(9) << "m[" << j << "] = m[" << j << "]-(" << m_j_i << "/" << m_i_i << ")*m[" << i << "]";
					// By definition, _m[j][i]-_m[i][i]/_m[i][i]*_m[j][i]=1.0
					_m[m_view.a0+j][m_view.b0+i] = 1.0;
					// Do not need to change coefficients from 0 to i+1.
					// They should all be 0 to build the echelon, but not used
					// in solving.
					for(std::size_t k = i+1; k < m_view.b1-m_view.b0; k++) {
						_m[m_view.a0+j][m_view.b0+k] =
							_m[m_view.a0+j][m_view.b0+k] -
							_m[m_view.a0+i][m_view.b0+k]/m_i_i*m_j_i;
					}
					{
						std::size_t _j = y_view.a0+j;
						std::size_t k = _m[_j].size()-1;
						_m[_j][k] = _m[_j][k] - _m[y_view.a0+i][k]/m_i_i*m_j_i;
					}
					CHEM_LOGV(9) << j << ": " << _m[j];
				}
			}

			CHEM_LOG(TRACE) << "Step 1:"
				<< aview(_m, m_view.a0, m_view.a1, m_view.b0, m_view.b1,
						y_view.a0, y_view.a1);
			
			X x = y_view.x;
			x[y_view.a1-1] = _m[y_view.a1-1][_m[y_view.a1-1].size()-1]
				/_m[m_view.a1-1][m_view.b1-1];

			CHEM_LOGV(9) << "x[" << y_view.a1-1 << "] = " << x[y_view.a1-1];
			for(std::size_t i = 1; i < y_view.a1-y_view.a0; i++) {
				x[y_view.a1-1-i] = _m[y_view.a1-1-i][_m[y_view.a1-1-i].size()-1];
				CHEM_LOGV(9) << "x[" << y_view.a1-1-i << "] = " << x[y_view.a1-1-i];
				for(std::size_t j = 1; j <= i; j++) {
					x[y_view.a1-1-i] -= _m[m_view.a1-1-i][m_view.b1-1-i+j] *
						x[y_view.a1-1-i+j];
				}
				x[y_view.a1-1-i] /= _m[m_view.a1-1-i][m_view.b1-1-i];
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
		typename Gauss<M, X>::result_type solve(const M& m, const X& y) {
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

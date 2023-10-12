#ifndef CHEMMISOL_LINEAR_H
#define CHEMMISOL_LINEAR_H

#include <array>
#include <vector>
#include <initializer_list>
#include <iostream>
#include <cmath>
#include <complex>
#include <iomanip>
#include "../logging.h"

/**
 * @file chemmisol/math/linear.h
 * File containing linear algebra features.
 */

namespace chemmisol {
	/*
	 * Fixed size array linear algebra
	 */

	/**
	 * Vector type hosted by a fixed size std::array.
	 */
	template<typename T, std::size_t N>
		using X = std::array<T, N>;

	/**
	 * Matrix type hosted by a fixed size std::array.
	 *
	 * The matrix contains N rows and P columns.
	 */
	template<typename T, std::size_t N, std::size_t P = N>
		using M = std::array<std::array<T, P>, N>;

	/**
	 * Vector type hosted by a dynamic size std::vector.
	 */
	template<typename T>
		using VecX = std::vector<T>;
	/**
	 * Matrix type hosted by a dynamic size std::vector.
	 */
	template<typename T>
		using VecM = std::vector<std::vector<T>>;

	template<typename M>
		struct MView {
			const M& m;
			const std::size_t a0;
			const std::size_t a1;
			const std::size_t b0;
			const std::size_t b1;

			MView(const M& m,
					std::size_t a0,
					std::size_t a1,
					std::size_t b0,
					std::size_t b1)
				: m(m), a0(a0), a1(a1), b0(b0), b1(b1) {
				}
			MView(const M& m) :
				MView(m, 0, 0, m.size(), m[m.size()-1].size()) {
				}
		};

	template<typename M>
		MView<M> mview(const M& m, 
					std::size_t a0,
					std::size_t a1,
					std::size_t b0,
					std::size_t b1) {
			return {m, a0, a1, b0, b1};
		}

	template<typename M>
		MView<M> mview(const M& m) {
			return {m};
		}

	template<typename X>
		struct XView {
			const X& x;
			const std::size_t a0;
			const std::size_t a1;

			XView(const X& x,
					std::size_t a0,
					std::size_t a1)
				: x(x), a0(a0), a1(a1) {
				}
			XView(const X& x) :
				XView(x, 0, x.size()) {
				}
		};

	template<typename X>
		XView<X> xview(const X& x,
				std::size_t a0,
				std::size_t a1) {
			return {x, a0, a1};
		}
	template<typename X>
		XView<X> xview(const X& x) {
			return {x};
		}

	template<typename X>
		inline MAKE_LOGGABLE(XView<X>, x_view, os) {
			os << "[";
			if(x_view.a1-x_view.a0 > 0)
				os << x_view.x[x_view.a0];
			for(std::size_t i = x_view.a0+1; i < x_view.a1; i++) {
				os << ", " << x_view.x[i];
			}
			os << "]";
			return os;
		}

	template<typename M>
		inline MAKE_LOGGABLE(MView<M>, m_view, os) {
			for(std::size_t i = m_view.a0; i < m_view.b0; i++) {
				os << std::endl << std::setw((int) std::log10(m_view.b0)+1) << i << ": " <<
					xview(m_view.m[i], m_view.a1, m_view.b1);
			}
			return os;
		}

	/**
	 * Computes `-x = [-x_0, ..., -x_N-1]`.
	 *
	 * @param x Vector of size N.
	 * @return Vector of size N.
	 */
	template<typename T, std::size_t N>
		X<T, N> operator-(const X<T, N>& x) {
			return -xview(x);
		}

	template<typename T, std::size_t N>
		X<T, N> operator-(const XView<X<T, N>>& x_view) {
			X<T, N> x1;
			for(std::size_t i = x_view.a0; i < x_view.a1; i++)
				x1[i] = -x_view.x[i];
			return x1;
		}

	/**
	 * Computes `x1+x2 = [x1_0+x2_0, ..., x1_N-1 + x2_N-1]`.
	 *
	 * @param x1 Vector of size N.
	 * @param x2 Vector of size N.
	 * @return Vector of size N.
	 */
	template<typename T, std::size_t N>
		X<T, N> operator+(const X<T, N>& x1, const X<T, N>& x2) {
			return xview(x1) + xview(x2);
		}

	template<typename T, std::size_t N>
		X<T, N> operator+(
				const XView<X<T, N>>& x_view1, const XView<X<T, N>>& x_view2) {
			X<T, N> x3;
			for(std::size_t i = x_view1.a0; i < x_view2.a1; i++)
				x3[i] = x_view1.x[i]+x_view2.x[i];
			return x3;
		}

	/**
	 * Computes m*x as a matrix product.
	 *
	 * @param m Matrix of size N*P.
	 * @param x Vector of size P.
	 * @return Vector of size N.
	 */
	template<typename T, std::size_t N, std::size_t P>
		X<T, N> operator*(const M<T, N, P>& m, const X<T, P>& x) {
			return mview(m) * xview(x);
		}

	template<typename T, std::size_t N, std::size_t P>
		X<T, N> operator*(
				const MView<M<T, N, P>>& m_view, const XView<X<T, P>>& x_view) {
			X<T, N> x1;
			for(std::size_t i = x_view.a0; i < x_view.a1; i++) {
				x1[i] = 0;
				for(std::size_t j = m_view.a1; j < m_view.b1; j++) {
					x1[i] += m_view.m[i][j] * x_view.x[j];
				}
			}
			return x1;
		}

	template<typename T, std::size_t N>
		X<T, N> operator*(const T& a, const X<T, N>& x) {
			X<T, N> result(x);
			for(auto& item : result)
				item *= a;
			return result;
		}

	template<typename T, std::size_t N>
		X<T, N> operator*(const T& a, const XView<X<T, N>>& x_view) {
			X<T, N> result;
			for(std::size_t i = x_view.a0; i < x_view.a1; i++) {
				result[i] = x_view.x[i]*a;
			}
			return result;
		}

	template<typename T, std::size_t N, std::size_t P>
		M<T, N, P> operator*(const T& a, const M<T, N, P>& m) {
			M<T, N, P> result(m);
			for(auto& row : result)
				for(auto& item : row)
					item *= a;
			return result;
		}

	template<typename T, std::size_t N, std::size_t P>
		M<T, N, P> operator*(const T& a, const MView<M<T, N, P>>& m_view) {
			M<T, N, P> result;
			for(std::size_t i = m_view.a0; i < m_view.b0; i++) {
				for(std::size_t j = m_view.a1; j < m_view.b1; j++) {
					result[i][j] = a * m_view.m[i][j];
				}
			}
			return result;
		}

	/**
	 * Computes the euclidean norm of the vector.
	 */
	template<typename T, std::size_t N>
		double norm(const X<T, N>& x) {
			T a = 0;
			for(auto& v : x) {
				a+=std::pow(v, 2);
			}
			return std::sqrt(a);
		}

	template<typename T, std::size_t N>
		double norm(const X<std::complex<T>, N>& x) {
			T a = 0;
			for(auto& v : x) {
				a+=std::pow(std::norm(v), 2);
			}
			return std::sqrt(a);
		}

	/**
	 * Computes the absolute value the vector.
	 */
	template<typename T, std::size_t N>
		X<T, N> abs(const X<T, N>& x) {
			X<T, N> abs_x;
			for(std::size_t i = 0; i < N; i++)
				abs_x[i] = std::abs(x[i]);
			return abs_x;
		}

	/**
	 * Implements the augment operation used by the [Gaussian elimination
	 * algortihm](https://en.wikipedia.org/wiki/Gaussian_elimination) for fixed
	 * size matrices and vectors.
	 *
	 * Returns a new matrix corresponding to the concatenation of m and x, where
	 * x represents the last column of the new matrix.
	 *
	 * @param m Matrix of size N*P.
	 * @param x Vector of size N.
	 * @return Matrix of size N*(P+1).
	 */
	template<typename T, std::size_t N, std::size_t P>
		M<T, N, P+1> augment(const M<T, N, P>& m, const X<T, N>& x) {
			return augment(mview(m), xview(x));
		};

	template<typename T, std::size_t N, std::size_t P>
		M<T, N, P+1> augment(
				const MView<M<T, N, P>>& m_view, const XView<X<T, N>>& x_view) {
			M<T, N, P+1> a;
			for(std::size_t i = m_view.a0; i < m_view.b0; i++) {
				for(std::size_t j = m_view.a1; j < m_view.b1; j++) {
					a[i][j] = m_view.m[i][j];
				}
			}
			for(std::size_t i = x_view.a0; i < x_view.a1; i++) {
				a[i][P] = x_view.x[i];
			}
			return a;
		};

	/**
	 * Fixed size vector stream output operator.
	 *
	 * The vector is serialized as "[x0, ..., xn]".
	 */
	template<typename T, std::size_t N>
		std::ostream& operator<<(std::ostream& o, const X<T, N>& x) {
			o << "[";
			for(std::size_t i = 0; i < N-1; i++)
				o << x[i] << ", ";
			if(N > 0)
				o << x[N-1];
			o << "]";
			return o;
		}

	/*
	 * Vector linear algebra
	 */

	/**
	 * Computes `-x = [-x_0, ..., -x_N-1]`.
	 *
	 * @param x Vector of size N.
	 * @return Vector of size N.
	 */
	template<typename T>
		VecX<T> operator-(const VecX<T>& x) {
			return -xview(x);
		}

	template<typename T>
		VecX<T> operator-(const XView<VecX<T>>& x_view) {
			VecX<T> x1(x_view.x.size());
			for(std::size_t i = x_view.a0; i < x_view.a1; i++)
				x1[i] = -x_view.x[i];
			return x1;
		}

	/**
	 * Computes `x1+x2 = [x1_0+x2_0, ..., x1_N-1 + x2_N-1]`.
	 *
	 * The behavior is unspecified if x1 and x2 are not the same size.
	 *
	 * @param x1 Vector of size N.
	 * @param x2 Vector of size N.
	 * @return Vector of size N.
	 */
	template<typename T>
		VecX<T> operator+(const VecX<T>& x1, const VecX<T>& x2) {
			return xview(x1)+xview(x2);
		}

	template<typename T>
		VecX<T> operator+(
				const XView<VecX<T>>& x_view1, const XView<VecX<T>>& x_view2) {
			VecX<T> x3(x_view1.x.size());
			for(std::size_t i = x_view1.a0; i < x_view1.a1; i++)
				x3[i] = x_view1.x[i]+x_view2.x[i];
			return x3;
		}

	/**
	 * Computes m*x as a matrix product.
	 *
	 * @param m Matrix of size N*P.
	 * @param x Vector of size P.
	 * @return Vector of size N.
	 */
	template<typename T>
		VecX<T> operator*(const VecM<T>& m, const VecX<T>& x) {
			return mview(m) * xview(x);
		}

	template<typename T>
		VecX<T> operator*(const MView<VecM<T>>& m_view, const XView<VecX<T>>& x_view) {
			VecX<T> x1(x_view.x.size());
			for(std::size_t i = x_view.a0; i < x_view.a1; i++) {
				x1[i] = 0;
				for(std::size_t j = m_view.a1; j < m_view.b1; j++) {
					x1[i] += m_view.m[i][j] * x_view.x[j];
				}
			}
			return x1;
		}

	template<typename T>
		VecX<T> operator*(const T& a, const VecX<T>& x) {
			VecX<T> result(x);
			for(auto& item : result)
				item *= a;
			return result;
		}

	template<typename T>
		VecX<T> operator*(const T& a, const XView<VecX<T>>& x_view) {
			VecX<T> result(x_view.x.size());
			for(std::size_t i = x_view.a0; i < x_view.a1; i++)
				result[i] = a * x_view.x[i];
			return result;
		}

	template<typename T>
		VecM<T> operator*(const T& a, const VecM<T>& m) {
			VecM<T> result(m);
			for(auto& row : result) {
				for(auto& item : row) {
					item *= a;
				}
			}
			return result;
		}

	template<typename T>
		VecM<T> operator*(const T& a, const MView<VecM<T>>& m_view) {
			VecM<T> result(m_view.m.size());
			for(std::size_t i = 0; i < result.size(); i++)
				result[i].resize(m_view.m[i].size());
			for(std::size_t i = m_view.a0; i < m_view.b0; i++) {
				for(std::size_t j = m_view.a1; j < m_view.b1; j++) {
					result[i][j] = a * m_view.m[i][j];
				}
			}
			return result;
		}

	/**
	 * Computes the euclidean norm of the vector.
	 */
	template<typename T>
		double norm(const VecX<T>& x) {
			double a = 0;
			for(auto& v : x) {
				a+=std::pow(v, 2);
			}
			return std::sqrt(a);
		}

	template<typename T>
		double norm(const VecX<std::complex<T>>& x) {
			double a = 0;
			for(auto& v : x) {
				a+=std::pow(std::norm(v), 2);
			}
			return std::sqrt(a);
		}

	/**
	 * Computes the absolute value the vector.
	 */
	template<typename T>
		VecX<T> abs(const VecX<T>& x) {
			VecX<T> abs_x(x);
			for(std::size_t i = 0; i < x.size(); i++)
				abs_x[i] = std::abs(x[i]);
			return abs_x;
		}

	/**
	 * Implements the augment operation used by the [Gaussian elimination
	 * algortihm](https://en.wikipedia.org/wiki/Gaussian_elimination) for
	 * dynamic size matrices and vectors.
	 *
	 * Returns a new matrix corresponding to the concatenation of m and x, where
	 * x represents the last column of the new matrix.
	 *
	 * The behavior is unspecified if the size of x is not equal to the count of
	 * rows in m.
	 *
	 * @param m Matrix of size N*P.
	 * @param x Vector of size N.
	 * @return Matrix of size N*(P+1).
	 */
	template<typename T>
		VecM<T>
		augment(const VecM<T>& m, const VecX<T>& x) {
			return augment(mview(m), xview(x));
		};

	template<typename T> VecM<T>
		augment(const MView<VecM<T>>& m_view, const XView<VecX<T>>& x_view) {
			VecM<T> a(m_view.m.size());
			for(std::size_t i = 0; i < a.size(); i++) {
				a[i].resize(m_view.m[i].size()+1);
			}

			for(std::size_t i = m_view.a0; i < m_view.b0; i++) {
				for(std::size_t j = m_view.a1; j < m_view.b1; j++) {
					a[i][j] = m_view.m[i][j];
				}
			}
			for(std::size_t i = x_view.a0; i < x_view.a1; i++) {
				a[i][m_view.m[i].size()] = x_view.x[i];
			}
			return a;
		};

	/**
	 * Dynamic size vector stream output operator.
	 *
	 * The vector is serialized as "[x0, ..., xn]".
	 */
	template<typename T>
		std::ostream& operator<<(std::ostream& o, const VecX<T>& x) {
			o << "[";
			for(std::size_t i = 0; i < x.size()-1; i++)
				o << x[i] << ", ";
			if(x.size() > 0)
				o << x[x.size()-1];
			o << "]";
			return o;
		}

	/*
	 * Scalar values algebra.
	 */

	/**
	 * Returns the norm of the scalar value x, i.e. its absolute value.
	 */
	template<typename T>
		double norm(const T& x) {
			return std::abs(x);
		}

	/**
	 * Returns the absolute value of the scalar value x.
	 */
	template<typename T>
		T abs(const T& x) {
			return std::abs(x);
		}
}
#endif

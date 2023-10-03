#ifndef CHEMMISOL_LINEAR_H
#define CHEMMISOL_LINEAR_H

#include <array>
#include <vector>
#include <initializer_list>
#include <iostream>
#include <cmath>
#include <complex>
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

	/**
	 * Computes m*x as a matrix multiplication.
	 *
	 * @param m Matrix of size N*P.
	 * @param x Vector of size P.
	 * @return Vector of size N.
	 */
	template<typename T, std::size_t N, std::size_t P>
		X<T, N> operator*(const M<T, N, P>& m, const X<T, P>& x) {
			X<T, N> x1;
			for(std::size_t i = 0; i < N; i++) {
				x1[i] = 0;
				for(std::size_t j = 0; j < P; j++) {
					x1[i] += m[i][j] * x[j];
				}
			}
			return x1;
		}

	/**
	 * Computes `-x = [-x_0, ..., -x_N-1]`.
	 *
	 * @param x Vector of size N.
	 * @return Vector of size N.
	 */
	template<typename T, std::size_t N>
		X<T, N> operator-(const X<T, N>& x) {
			X<T, N> x1;
			for(std::size_t i = 0; i < N; i++)
				x1[i] = -x[i];
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
			X<T, N> x3;
			for(std::size_t i = 0; i < N; i++)
				x3[i] = x1[i]+x2[i];
			return x3;
		}

	template<typename T, std::size_t N>
		X<T, N> operator*(const T& a, const X<T, N>& x) {
			X<T, N> result(x);
			for(auto& item : result)
				item *= a;
			return result;
		}

	template<typename T, std::size_t N, std::size_t P>
		M<T, N, P> operator*(const T& a, const M<T, N, P>& m) {
			M<T, N, P> result(m);
			for(auto& row : result) {
				for(auto& item : row) {
					item *= a;
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
	 *
	 * @tparam _M Fixed size matrix type.
	 */
	template<typename _M>
		M<
			typename _M::value_type::value_type,
			std::tuple_size<_M>::value,
			std::tuple_size<typename _M::value_type>::value+1>
		augment(
				const _M& m,
				const X<typename _M::value_type::value_type, std::tuple_size<_M>::value>& x
				) {
			M<
				typename _M::value_type::value_type,
				std::tuple_size<_M>::value,
				std::tuple_size<typename _M::value_type>::value+1> a;
				for(std::size_t i = 0; i < std::tuple_size<_M>::value; i++) {
					for(std::size_t j = 0; j < std::tuple_size<typename _M::value_type>::value; j++) {
						a[i][j] = m[i][j];
					}
					a[i][std::tuple_size<_M>::value] = x[i];
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
	 * Computes m*x as a matrix multiplication.
	 *
	 * @param m Matrix of size N*P.
	 * @param x Vector of size P.
	 * @return Vector of size N.
	 */
	template<typename T>
		VecX<T> operator*(const VecM<T>& m, const VecX<T>& x) {
			VecX<T> x1(x.size());
			for(std::size_t i = 0; i < x.size(); i++) {
				x1[i] = 0;
				for(std::size_t j = 0; j < m[i].size(); j++) {
					x1[i] += m[i][j] * x[j];
				}
			}
			return x1;
		}

	/**
	 * Computes `-x = [-x_0, ..., -x_N-1]`.
	 *
	 * @param x Vector of size N.
	 * @return Vector of size N.
	 */
	template<typename T>
		VecX<T> operator-(const VecX<T>& x) {
			VecX<T> x1(x.size());
			for(std::size_t i = 0; i < x.size(); i++)
				x1[i] = -x[i];
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
			VecX<T> x3(x1.size());
			for(std::size_t i = 0; i < x1.size(); i++)
				x3[i] = x1[i]+x2[i];
			return x3;
		}

	template<typename T>
		VecX<T> operator*(const T& a, const VecX<T>& x) {
			VecX<T> result(x);
			for(auto& item : result)
				item *= a;
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
			VecM<T> a(m.size());
			for(std::size_t i = 0; i < m.size(); i++) {
				a[i].resize(m[i].size()+1);
				for(std::size_t j = 0; j < m[i].size(); j++) {
					a[i][j] = m[i][j];
				}
				a[i][m[i].size()] = x[i];
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

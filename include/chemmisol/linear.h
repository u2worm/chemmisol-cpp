#include <array>
#include <vector>
#include <initializer_list>
#include <iostream>
#include <cmath>

namespace chemmisol {
	template<typename T>
		T log(const T& v) {
			return std::log10(v);
		}

	template<typename T>
		T ln(const T& v) {
			return std::log(v);
		}
	const double ln10 = ln(10);

	/*
	 * Fixed size array linear algebra
	 */

	template<typename T, int N>
		struct X : std::array<T, N> {
			typedef T value_type;
			static constexpr int n = N;

			X() = default;
			X(const std::initializer_list<T>& a) {
				std::size_t i = 0;
				for(auto& item : a)
					(*this)[i++] = item;
			}
		};

	template<typename T, int N, int P = N>
		struct M : std::array<std::array<T, P>, N> {
			typedef T coef_type;
			static constexpr int n = N;
			static constexpr int p = P;

			M() = default;
			M(const std::initializer_list<std::initializer_list<T>>& a) {
				std::size_t i = 0;
				for(auto& line : a) {
					std::size_t j = 0;
					for(auto& item : line) {
					(*this)[i][j++] = item;
					}
					++i;
				}
			}
		};

	template<typename T>
		using VecX = std::vector<T>;
	template<typename T>
		using VecM = std::vector<std::vector<T>>;

	template<typename T, int N, int P>
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

	template<typename T, int N>
		X<T, N> operator-(const X<T, N>& x) {
			X<T, N> x1;
			for(std::size_t i = 0; i < N; i++)
				x1[i] = -x[i];
			return x1;
		}

	template<typename T, int N>
		X<T, N> operator+(const X<T, N>& x1, const X<T, N>& x2) {
			X<T, N> x3;
			for(std::size_t i = 0; i < N; i++)
				x3[i] = x1[i]+x2[i];
			return x3;
		}

	template<typename T, int N>
		std::ostream& operator<<(std::ostream& o, const X<T, N>& x) {
			for(std::size_t i = 0; i < N-1; i++)
				o << x[i] << ", ";
			o << x[N-1];
			return o;
		}

	template<typename T, int N>
		double norm(const X<T, N>& x) {
			double a = 0;
			for(auto& v : x) {
				a+=std::pow(v, 2);
			}
			return std::sqrt(a);
		}

	template<typename T, int N>
		X<T, N> abs(const X<T, N>& x) {
			X<T, N> abs_x;
			for(std::size_t i = 0; i < N; i++)
				abs_x[i] = std::abs(x[i]);
			return abs_x;
		}


	template<typename _M>
		M<typename _M::coef_type, _M::n, _M::p+1>
		augment(const _M& m, const X<typename _M::coef_type, _M::n>& x) {
			M<typename _M::coef_type, _M::n, _M::p+1> a;
				for(std::size_t i = 0; i < _M::n; i++) {
					for(std::size_t j = 0; j < _M::n; j++) {
						a[i][j] = m[i][j];
					}
					a[i][_M::n] = x[i];
				}
				return a;
		};

	/*
	 * Vector linear algebra
	 */


	template<typename T>
		std::ostream& operator<<(std::ostream& o, const VecX<T>& x) {
			for(std::size_t i = 0; i < x.size()-1; i++)
				o << x[i] << ", ";
			o << x[x.size()-1];
			return o;
		}

	template<typename T>
		VecX<T> operator-(const VecX<T>& x) {
			VecX<T> x1(x.size());
			for(std::size_t i = 0; i < x.size(); i++)
				x1[i] = -x[i];
			return x1;
		}

	template<typename T>
		VecX<T> operator+(const VecX<T>& x1, const VecX<T>& x2) {
			VecX<T> x3(x1.size());
			for(std::size_t i = 0; i < x1.size(); i++)
				x3[i] = x1[i]+x2[i];
			return x3;
		}

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

	template<typename T>
		double norm(const VecX<T>& x) {
			double a = 0;
			for(auto& v : x) {
				a+=std::pow(v, 2);
			}
			return std::sqrt(a);
		}

	template<typename T>
		VecX<T> abs(const VecX<T>& x) {
			VecX<T> abs_x(x);
			for(std::size_t i = 0; i < x.size(); i++)
				abs_x[i] = std::abs(x[i]);
			return abs_x;
		}

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

}

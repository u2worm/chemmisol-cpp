#ifndef CHEMMISOL_ARITHMETICS_H
#define CHEMMISOL_ARITHMETICS_H

#include <complex>
#include <vector>
#include <iostream>
#include <numbers>

namespace chemmisol {
	constexpr double pi = 3.141592653589793238462643383279502884;

	template<typename T>
		std::vector<std::complex<T>> unit_roots(int n) {
			std::vector<std::complex<T>> r;
			std::complex<T> single_root = {
				std::cos(2*pi/n),
				std::sin(2*pi/n),
			};

			for(int i = 0; i < n; i++) {
				r.push_back(std::pow(single_root, i));
			}
			return r;
		}

	template<typename T>
		std::vector<std::complex<T>> roots(const std::complex<T>& c, int n) {
			std::vector<std::complex<T>> r;
			T theta = std::arg(c);
			std::complex<T> single_root = 
				std::polar(std::pow(std::abs(c), T(1.0)/n), theta/n);
			auto units = unit_roots<T>(n);
			for(auto unit : units) {
				r.push_back(single_root * unit);
			}
			return r;
		}
}
#endif

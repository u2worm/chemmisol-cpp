#include "chemmisol/math/gauss.h"

namespace chemmisol { namespace gauss {
	double Gauss<double,double>::solve(const double& f, const double& y) {
		return y/f;
	}

	std::complex<double> Gauss<std::complex<double>, std::complex<double>>::solve(
			const std::complex<double>& f, const std::complex<double>& y) {
		return y/f;
	}
}}

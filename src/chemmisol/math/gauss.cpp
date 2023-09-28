#include "chemmisol/math/gauss.h"

namespace chemmisol { namespace gauss {
	template<>
		double solve<double,double>(const double& f, const double& y) {
			return y/f;
		}

	template<>
		std::complex<double> solve<std::complex<double>, std::complex<double>>(
				const std::complex<double>& f, const std::complex<double>& y) {
			return y/f;
		}
}}

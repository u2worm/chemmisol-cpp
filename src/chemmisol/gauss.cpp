#include "chemmisol/gauss.h"

namespace chemmisol { namespace gauss {
	template<>
		double solve<double,double>(const double& f, const double& y) {
			return y/f;
		}
}}

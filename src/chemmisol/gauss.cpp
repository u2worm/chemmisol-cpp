#include "chemmisol/gauss.h"

namespace mineral { namespace gauss {
	template<>
		double solve<double,double>(const double& f, const double& y) {
			return y/f;
		}
}}

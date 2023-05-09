#include "gauss.h"

namespace gauss {
	template<>
		double solve<double,double>(const double& f, const double& y) {
			return y/f;
		}
}

#include "chemmisol/chemical/reaction.h"

namespace chemmisol {
	bool operator==(const Reagent& c1, const Reagent& c2) {
		return c1.name == c2.name
			&& c1.phase == c2.phase
			&& c1.coefficient == c2.coefficient;
	}
}

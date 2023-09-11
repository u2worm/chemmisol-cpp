#include "chemmisol/chemical/system.h"
#include <iomanip>

namespace chemmisol {
	std::ostream& operator<<(std::ostream& o, const ChemicalSpecies& s);

	bool operator==(const Reagent& c1, const Reagent& c2);

	void print_df(
			const ChemicalSystem& system,
			const solver::F& f,
			std::vector<std::vector<double>>& df,
			std::size_t unfixed_components);
}


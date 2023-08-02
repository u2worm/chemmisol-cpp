#include "chemmisol/chemical/species.h"
#include <iostream>

namespace chemmisol {
	std::ostream& operator<<(std::ostream& o, const Phase& phase) {
		switch(phase) {
			case SOLVENT:
				o << "SOLVENT";
				break;
			case AQUEOUS:
				o << "AQUEOUS";
				break;
			case MINERAL:
				o << "MINERAL";
		}
		return o;
	}

	const double AqueousSpecies::V = 1*l;

	double MineralSpecies::sites_count(
			double solid_concentration,
			double specific_surface_area,
			double site_concentration) {
		return AqueousSpecies::V*solid_concentration*specific_surface_area*site_concentration;
	}

}

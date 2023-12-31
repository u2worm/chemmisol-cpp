#include "chemmisol/tests/chemical/system/utils.h"

namespace chemmisol {
	std::ostream& operator<<(std::ostream& o, const ChemicalSpecies& s) {
		o << "{i:" << s.getIndex() << ", name:" << s.getName() << ", phase:" << s.getPhase() << "}";
		return o;
	}

	bool operator==(const Reagent& c1, const Reagent& c2) {
		return c1.name == c2.name
			&& c1.phase == c2.phase
			&& c1.coefficient == c2.coefficient;
	}

	void print_df(
			const ChemicalSystem& system,
			const solver::ReducedChemicalSystem<solver::X>& reduced_system,
			std::vector<std::vector<double>>& df,
			std::size_t unfixed_components) {

		std::ostringstream s;
		s << "  " << std::setw(6) << "";
		for(const auto& species : system.getSpecies()) {
			if(species->getName() != "H2O"
					&& species->getName() != "H+") {
				s << std::setw(16) << species->getName();
			}
		}
		CHEM_LOGV(5) << s.str();
		for(const auto& component : system.getComponents()) {
			s.str("");
			if(component->getSpecies()->getName() != "H2O"
					&& component->getSpecies()->getName() != "H+") {
				s << "  " << std::setw(6) << component->getSpecies()->getName();
				for(const auto& item : df[reduced_system.componentsIndexes()[component->getIndex()]])
					s << std::setw(16) << std::scientific << item;
				CHEM_LOGV(5) << s.str();
			}
		}
		for(const auto& reaction : system.getReactions()) {
			s.str("");
			s << "  " << std::setw(6) << reaction->getName();
			for(const auto& item : df[unfixed_components + reaction->getIndex()])
				s << std::setw(16) << std::scientific << item;
			CHEM_LOGV(5) << s.str();
		}
	}

	std::vector<std::minstd_rand::result_type> gen_seeds() {
		std::vector<int> gen;
		for(int i = 0; i < 10; i++)
			gen.push_back(i);
		std::seed_seq seq(gen.begin(), gen.end());

		std::vector<std::minstd_rand::result_type> seeds(10);
		seq.generate(seeds.begin(), seeds.end());

		return seeds;
	}
}


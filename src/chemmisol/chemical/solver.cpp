#include "chemmisol/chemical/system.h"
#include "chemmisol/math/newton.h"

namespace chemmisol {
	namespace solver {
		F::F(const ChemicalSystem& system)
			: system(system), init_concentrations(system.getSpecies().size()) {
				for(auto& species : system.getSpecies()) {
					init_concentrations[species->getIndex()] = species->concentration();
				}
			}

		X F::concentrations(const X& extent) const {
			// Init concentrations of all components
			std::vector<double> concentrations = init_concentrations;
			for(auto& reaction : system.getReactions()) {
				const std::vector<double>& coefficients
					= system.getReactionMatrix()[reaction->getIndex()];
				CHEM_LOG(TRACE) << "Calc concentrations from reaction " << reaction->getName();
				for(
						std::size_t species_index = 0;
						species_index < coefficients.size();
						species_index++) {

					CHEM_LOG(TRACE) << "  " <<
						system.getSpecies(species_index).getName() << ": " <<
						concentrations[species_index] << " + " <<
						coefficients[species_index]*extent[reaction->getIndex()] << " = " <<
						system.getSpecies(species_index)
								.concentration(
									concentrations[species_index],
									coefficients[species_index]*extent[reaction->getIndex()]
									);
					concentrations[species_index] =
						system.getSpecies(species_index)
								.concentration(
									concentrations[species_index],
									coefficients[species_index]*extent[reaction->getIndex()]
									);
				}
			}
			// Activities of all components resulting from the extent of all
			// reactions
			return concentrations;
		}
		X F::f(const X& extent) const {
			X concentrations = this->concentrations(extent);
			CHEM_LOG(TRACE) << "Init F concentrations: ";
			for(auto& species : system.getSpecies())
				CHEM_LOG(TRACE) << "  " << species->getName() << " " << concentrations[species->getIndex()];

			X f(extent.size());
			for(auto& reaction : system.getReactions()) {
				const std::vector<double>& coefficients
					= system.getReactionMatrix()[reaction->getIndex()];
				for(
						std::size_t species_index = 0;
						species_index < coefficients.size();
						species_index++) {
					auto& species = system.getSpecies(species_index);
					f[reaction->getIndex()] +=
						// Stoichiometric coefficient
						coefficients[species_index] * log(
								// Current activity of the species according
								// to the provided reaction extent
								species.activity(concentrations[species_index])
								);
				}
				// The previous sum should be equal to log(K)
				f[reaction->getIndex()] += reaction->getLogK();
			}
			return f;
		}
		M F::df(const X& extent) const {
			X current_concentrations = this->concentrations(extent);
			M jacobian(extent.size()); // reactions count

			for(auto& reaction : system.getReactions()) {
				jacobian[reaction->getIndex()].resize(extent.size());
				const std::vector<double>& coefficients_i
					= system.getReactionMatrix()[reaction->getIndex()];
				for(auto& dreaction : system.getReactions()) {
					const std::vector<double>& coefficients_j
						= system.getReactionMatrix()[dreaction->getIndex()];
					double dfi_dxj = 0;
					for(const auto& species : system.getSpecies()) {
						dfi_dxj += coefficients_i[species->getIndex()] *
							species->Dactivity(
									current_concentrations,
									coefficients_j[species->getIndex()]) /
							species->activity(
									current_concentrations[species->getIndex()]
									) / ln10;
					}
					jacobian[reaction->getIndex()][dreaction->getIndex()] = dfi_dxj;
				}
			}

			return jacobian;
		}

		X solve(const ChemicalSystem& system) {
			F f(system);
			// Initially all to 0
			X init_extent(system.getReactions().size());
			for(const auto& reaction : system.getReactions()) {
				init_extent[reaction->getIndex()]
					= system.getInitialGuessExtent(reaction->getName());
			}
			return Newton<X, M>(
					init_extent,
					[&f] (const X& x) {return f.f(x);},
					[&f] (const X& x) {return f.df(x);}
					).solve_iter(system.getMaxIteration());
		}
	}

	double GuessF::operator()(const double &extent) const {
		double f = 0;
		for(auto& reagent : reaction.getReagents()) {
			f += reagent.coefficient * log(system.getSpecies(reagent.name).activity(
					current_concentrations,
					reagent.coefficient * extent
					));
		}
		f += reaction.getLogK();
		return f;
	}
}

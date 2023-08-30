#include "chemmisol/chemical/system.h"
#include "chemmisol/math/newton.h"

namespace chemmisol {
	namespace solver {
		const std::size_t F::INVALID_INDEX = -1;

		F::F(const ChemicalSystem& system) :
			system(system),
			x_size(system.getSpecies().size()),
			f_x_size(system.getComponents().size()
					+ system.getReactions().size()),
			reaction_offset(system.getComponents().size()),
			components_indexes(system.getComponents().size()),
			species_indexes(system.getSpecies().size()),
			fixed_activities(system.getComponents().size()) {
				for(const auto& component : system.getComponents())
					components_indexes[component->getIndex()] = component->getIndex();
				for(const auto& species : system.getSpecies())
					species_indexes[species->getIndex()]
						= species->getIndex();

				std::vector<std::size_t> species_offsets(system.getSpecies().size());
				for(const auto& component : system.getComponents()) {
					if(component->isFixed()) {
						fixed_activities[component->getIndex()]
							= component->getSpecies()->activity();
						--x_size;
						--f_x_size;
						--reaction_offset;
						// This component index cannot be used
						components_indexes[component->getIndex()] = INVALID_INDEX;
						for(
								std::size_t i = component->getIndex()+1;
								i < components_indexes.size(); i++) {
							--components_indexes[i];
						}
						// This species index cannot be used
						species_indexes[component->getSpecies()->getIndex()] = INVALID_INDEX;
						for(
								std::size_t i = component->getSpecies()->getIndex()+1;
								i < species_indexes.size(); i++) {
							--species_indexes[i];
							++species_offsets[i];
						}
					}
				}

				revert_species_indexes.resize(x_size);
				for(const auto& species : system.getSpecies()) {
					if(species_indexes[species->getIndex()] != INVALID_INDEX)
						revert_species_indexes[species_indexes[species->getIndex()]]
							= species_indexes[species->getIndex()]+species_offsets[species->getIndex()];
				}
			}

		X F::reducedActivities() const {
			X x(x_size);
			for(const auto& species : system.getSpecies())
				if(species_indexes[species->getIndex()] != INVALID_INDEX)
					x[species_indexes[species->getIndex()]] = species->activity();
			return x;
		}

		X F::completeActivities(const X& reduced_activities) const {
			X complete_activities(system.getSpecies().size());
			for(const auto& component : system.getComponents())
				if(component->isFixed())
					complete_activities[component->getSpecies()->getIndex()]
						= component->getSpecies()->activity();
			for(std::size_t i = 0; i < reduced_activities.size(); i++) {
				complete_activities[revert_species_indexes[i]]
					= reduced_activities[i];
			}
			return complete_activities;
		}

		X F::f(const X& activities) const {
			X f_x(f_x_size);
			// Complete activities vector with an entry for each species, so
			// that it is compatible with ChemicalSystem::massConservationLaw()
			// and ChemicalSystem::distanceToEquilibrium().
			// TODO: improve this, too many copies for nothing.
			X complete_activities = completeActivities(activities);
			{
				X mass_conservation_results(system.getComponents().size());
				// Actual total quantity of each component
				system.massConservationLaw(
						complete_activities, mass_conservation_results);
				for(std::size_t i = 0; i < mass_conservation_results.size(); i++) {
					if(components_indexes[i] != INVALID_INDEX)
						f_x[components_indexes[i]] = mass_conservation_results[i];
				}
			}
			
			for(const auto& reaction : system.getReactions()) {
				f_x[reaction_offset+reaction->getIndex()]
					= system.distanceToEquilibrium(complete_activities, *reaction);
			}
			return f_x;
		}

		M F::df(const X& activities) const {
			M jacobian(f_x_size);
			std::vector<std::size_t> component_offsets;

			// Begin mass conservation law equations
			for(auto& component : system.getComponents()) {
				if(components_indexes[component->getIndex()] != INVALID_INDEX) {
					auto& d_f = jacobian[components_indexes[component->getIndex()]];
					d_f.resize(x_size);
					d_f[species_indexes[component->getSpecies()->getIndex()]] = 1.0;
				}
			}
			for(const auto& reaction : system.getReactions()) {
				// Species produced by this reaction
				const ChemicalSpeciesReagent& species
					= system.getSpeciesReagent(*reaction);
				// The index of the species produced by the reaction is
				// necessarily valid.
				for(const auto& reagent : system.getComponentReagents(*reaction)) {
					if(components_indexes[reagent.component->getIndex()] != INVALID_INDEX)
						jacobian
							[components_indexes[reagent.component->getIndex()]]
							[species_indexes[species.species->getIndex()]]
								= reagent.coefficient/(-species.coefficient);
				}
			}
			// End mass conservation law equations

			// Begin equilibrium equations
			for(const auto& reaction : system.getReactions()) {
				auto& d_f = jacobian[reaction_offset + reaction->getIndex()];
				d_f.resize(x_size);
				// dx_species: activity variable from which the current reaction
				// d_f is derived
				for(const auto& dx_species : system.getSpecies()) {
					// No derivative to compute for fixed species (that
					// correspond to invalid indexes)
					if(species_indexes[dx_species->getIndex()] != INVALID_INDEX) {
						d_f[species_indexes[dx_species->getIndex()]] = 0.0;
						double d_reactives = -reaction->getK();
						double d_products = 1.0;
						double species_coefficient_in_reactions;
						bool species_in_reaction = false;

						// Process the produced species
						const auto& reagent = system.getSpeciesReagent(*reaction);
						if(reagent.species == dx_species.get()) {
							// The current reagent is the variable from which
							// d_f is derived
							species_coefficient_in_reactions = reagent.coefficient;
							species_in_reaction = true;
						} else {
							// The produced reagent cannot be fixed (currently)
							// so it necessarily corresponds to a valid index
							if(activities[species_indexes[reagent.species->getIndex()]] > 0.0) {
								if(reagent.coefficient < 0.0) {
									d_products *= std::pow(
											activities[species_indexes[reagent.species->getIndex()]],
											-reagent.coefficient
											);
								} else {
									d_reactives *= std::pow(
											activities[species_indexes[reagent.species->getIndex()]],
											reagent.coefficient
											);
								}
							} else {
								if(reagent.coefficient < 0.0) {
									d_products = 0.0;
								} else {
									d_reactives = 0.0;
								}
							}
						}

						// Process reaction components
						for(const auto& reagent : system.getComponentReagents(*reaction)) {
							if(reagent.component->getSpecies() == dx_species.get()) {
								species_coefficient_in_reactions = reagent.coefficient;
								species_in_reaction = true;
							} else {
								double species_activity = 0.0;
								if(species_indexes[reagent.component->getSpecies()->getIndex()] != INVALID_INDEX) {
									species_activity = activities[species_indexes[reagent.component->getSpecies()->getIndex()]];
								} else {
									species_activity = fixed_activities[reagent.component->getIndex()];
								}
								if(species_activity > 0.0) {
									if(reagent.coefficient < 0.0) {
										d_products *= std::pow(
												species_activity,
												-reagent.coefficient
												);
									} else {
										d_reactives *= std::pow(
												species_activity,
												reagent.coefficient
												);
									}
								} else {
									if(reagent.coefficient < 0.0) {
										d_products = 0.0;
									} else {
										d_reactives = 0.0;
									}
								}
							}
						}
						// If the species is not in the reaction,
						// d_f/dx_species = 0.0
						if(species_in_reaction) {
							std::size_t index = species_indexes[dx_species->getIndex()];
							double activity = activities[index];
							// d x^a/dx = a * x^(a-1)
							if(species_coefficient_in_reactions < 0)
								d_f[index] = -species_coefficient_in_reactions
									* std::pow(
											activity,
											-species_coefficient_in_reactions-1
											) * d_products;
							else
								d_f[index] = species_coefficient_in_reactions
									* std::pow(
											activity,
											species_coefficient_in_reactions-1
											) * d_reactives;
						}
					}
				}
			}
			// End equilibrium equations
			return jacobian;
		}

		X solve(const ChemicalSystem& system) {
			F f(system);
			// Initial activities in the system
			X reduced_activities = f.reducedActivities();
			reduced_activities = AbsoluteNewton<X, M>(
					reduced_activities,
					[&f] (const X& x) {return f.f(x);},
					[&f] (const X& x) {return f.df(x);}
					).solve_iter(system.getMaxIteration());
			return f.completeActivities(reduced_activities);
		}
	}
}

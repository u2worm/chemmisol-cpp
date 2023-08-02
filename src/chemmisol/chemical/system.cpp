#include "chemmisol/chemical/system.h"
#include "chemmisol/math/regula_falsi.h"

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cassert>

namespace chemmisol {
	/*
	 *const double ChemicalSystem::bulk_density = 1.17*gram/cm3;
	 *const double ChemicalSystem::V = std::pow(1*cm, 3);
	 *const double ChemicalSystem::mineral_weight  = V*bulk_density;
	 */

	ChemicalSystem::ChemicalSystem(
					double solid_concentration,
					double specific_surface_area,
					double site_concentration,
					const std::string& surface_component
					) :
		solid_concentration(solid_concentration),
		specific_surface_area(specific_surface_area),
		site_concentration(site_concentration) {
			addComponent(surface_component, MINERAL, 1);
		}

	ChemicalSystem::ChemicalSystem(const ChemicalSystem& other) :
				component_index(other.component_index),
				species_index(other.species_index),
				reaction_index(other.reaction_index),
				species(other.species.size()),
				components(other.components.size()),
				reactions(other.reactions.size()),
				compiled_reactions(other.compiled_reactions.size()),
				max_iteration(other.max_iteration),
				solid_concentration(other.solid_concentration),
				specific_surface_area(other.specific_surface_area),
				site_concentration(other.site_concentration) {
					for(const auto& component : other.getComponents()) {
						if(component->isFixed()) {
							this->fixComponent(
									component->getSpecies()->getName(),
									component->getSpecies()->getPhase(),
									component->getSpecies()->concentration(),
									component->getSpecies()->getIndex(),
									component->getIndex()
									);
						} else {
							this->addComponent(
									component->getSpecies()->getName(),
									component->getSpecies()->getPhase(),
									component->getSpecies()->concentration(),
									component->getSpecies()->getIndex(),
									component->getIndex()
									);
						}
					}
					for(const auto& species : other.getSpecies()) {
						if(!this->species[species->getIndex()]) {
							// If false, the unique_ptr is empty so the species
							// has not been added as a Component
							this->addSpecies(
									species->getName(),
									species->getPhase(),
									species->concentration(),
									species->getIndex()
									);
						}
					}
					for(const auto& reaction : other.getReactions()) {
						this->addReaction(
								reaction->getName(),
								reaction->getLogK(),
								reaction->getReagents(),
								reaction->getIndex()
								);
					}
					for(const auto& reaction : other.compiled_reactions) {
						this->compiled_reactions[reaction.reaction->getIndex()]
							= {this->reactions[reaction.reaction->getIndex()].get()};
						CompiledReaction& compiled_reaction
							= this->compiled_reactions[reaction.reaction->getIndex()];
						compiled_reaction.produced_species = {
							reaction.produced_species.coefficient,
							this->species[reaction.produced_species.species->getIndex()].get()
						};
						for(auto& component : reaction.components) {
							compiled_reaction.components.emplace_back(
									component.coefficient,
									this->components[component.component->getIndex()].get()
									);
						}
					}
					this->reaction_matrix = other.reaction_matrix;
			}

	void ChemicalSystem::setUp() {
		initReactionMatrix();
		compiled_reactions.clear();
		compiled_reactions.resize(reactions.size());
		for(const auto& reaction : reactions)
			compile(reaction.get());
	}

	void ChemicalSystem::addSpecies(ChemicalSpecies* species, std::size_t index) {
		this->species.resize(std::max(this->species.size(), index+1));
		this->species[index].reset(species);
		this->species_by_name[species->getName()] = species;
	}

	void ChemicalSystem::addComponent(Component* component,
			std::size_t species_index, std::size_t component_index) {
		addSpecies(component->getSpecies(), species_index);
		components.resize(std::max(components.size(), component_index+1));
		components[component_index].reset(component);
		components_by_name[component->getSpecies()->getName()] = component;
	}

	void ChemicalSystem::addReaction(Reaction* reaction, std::size_t index) {
		reactions.resize(std::max(reactions.size(), index+1));
		reactions[index].reset(reaction);
		reactions_by_name[reaction->getName()] = reaction;
	}

	void ChemicalSystem::addReaction(
			std::string name,
			double K,
			std::vector<Reagent> reactives,
			std::size_t index) {
		addReaction(new Reaction(
				name, index, K, reactives
				), index);
	}
	void ChemicalSystem::addReaction(
			std::string name,
			double K,
			std::vector<Reagent> reactives
			) {
		addReaction(name, K, reactives, reaction_index++);
	}

	void ChemicalSystem::fixComponent(
			const std::string& name,
			double concentration
			) {
		fixComponent(name, AQUEOUS, concentration);
	}

	void ChemicalSystem::fixComponent(
			const std::string& name,
			Phase phase,
			double concentration,
			std::size_t species_index,
			std::size_t component_index
			) {
		FixedChemicalSpecies* fixed_species;

		switch(phase) {
			case AQUEOUS:
				fixed_species = new FixedAqueousSpecies(
						name, species_index,
						concentration
						);
				break;
			case MINERAL:
				fixed_species = new FixedMineralSpecies(
						name, species_index,
						solid_concentration, specific_surface_area,
						site_concentration,
						concentration
						);
				break;
			case SOLVENT:
				fixed_species = new Solvent(name, species_index);
				break;
			default:
				// Should not append
				fixed_species = nullptr;
		}
		FixedComponent* fixed_component = new FixedComponent(
				fixed_species, component_index,
				fixed_species->ChemicalSpecies::quantity());

		auto existing_component = components_by_name.find(name);
		if (existing_component != components_by_name.end()) {
			components[component_index].reset(fixed_component);
			species[species_index].reset(fixed_species);
			existing_component->second = fixed_component;
			species_by_name.find(fixed_species->getName())->second = fixed_species;
		} else {
			addComponent(fixed_component, species_index, component_index);
		}
	}
	void ChemicalSystem::fixComponent(
			const std::string& name,
			Phase phase,
			double concentration
			) {

		auto existing_component = components_by_name.find(name);
		if (existing_component != components_by_name.end()) {
			fixComponent(
					name, phase, concentration, 
					existing_component->second->getSpecies()->getIndex(),
					existing_component->second->getIndex()
					);
		} else {
			fixComponent(
					name, phase, concentration, 
					species_index++, component_index++
					);
		}
	}

	void ChemicalSystem::addComponent(const std::string& name, double concentration) {
		addComponent(name, AQUEOUS, concentration);
	}

	void ChemicalSystem::addComponent(
			const std::string& name, Phase phase, double concentration) {
		addComponent(name, phase, concentration, species_index++, component_index++);
		}

	void ChemicalSystem::addComponent(
			const std::string& name, Phase phase, double concentration,
			std::size_t species_index, std::size_t component_index) {
		ChemicalSpecies* species;
		switch(phase) {
			case AQUEOUS:
				species = new AqueousSpecies(name, species_index, concentration);
				break;
			case MINERAL:
				species = new MineralSpecies(
						name, species_index,
						solid_concentration, specific_surface_area, site_concentration,
						concentration);
				break;
			case SOLVENT:
				species = new Solvent(name, species_index);
		}
		addComponent(
				new Component(species, component_index, species->quantity()),
				species_index, component_index
				);
	}

	void ChemicalSystem::addSpecies(
			const std::string& name, Phase phase) {
		addSpecies(name, phase, 0, species_index++);
	}

	void ChemicalSystem::addSpecies(
			const std::string& name, Phase phase, double concentration,
			std::size_t species_index) {
		ChemicalSpecies* species;
		switch(phase) {
			case AQUEOUS:
				species = new AqueousSpecies(name, species_index, concentration);
				break;
			case MINERAL:
				species = new MineralSpecies(
						name, species_index,
						solid_concentration, specific_surface_area, site_concentration,
						concentration);
				break;
			case SOLVENT:
				species = new Solvent(name, species_index);
		}
		addSpecies(species, species_index);
	}

	void ChemicalSystem::initPH(double pH) {
		addComponent("H+", std::pow(10, -pH)*mol/l);
		//addComponent("OH-", std::pow(10, -pH)*mol/l);
	}

	void ChemicalSystem::fixPH(double pH) {
		fixComponent(
					"H+", AQUEOUS,
					std::pow(10, -pH)*mol/l
				);
		//addComponent(
					//"OH-", AQUEOUS,
					//std::pow(10, getReaction("OH-").getLogK())
						/// std::pow(10, -pH)*mol/l
				//);
		//initPH(pH);
		// A fake reaction that produces (or consumes) as much H+ as needed to
		// fix the pH to the provided value
		//addReaction("pH", -pH, {{"H+", -1}});
	}

	double ChemicalSystem::getPH() const {
		return -log(getComponent("H+").getSpecies()->concentration());
	}

	const Reaction& ChemicalSystem::getReaction(const std::string& name) const {
		return *reactions_by_name.find(name)->second;
	}

	const Component& ChemicalSystem::getComponent(const std::string& name) const {
		return *components_by_name.find(name)->second;
	}

	const Component& ChemicalSystem::getComponent(const std::size_t& id) const {
		return *components[id];
	}

	const ChemicalSpecies& ChemicalSystem::getSpecies(const std::string& name) const {
		return *species_by_name.find(name)->second;
	}

	const ChemicalSpecies& ChemicalSystem::getSpecies(const std::size_t& id) const {
		return *species[id];
	}


	void ChemicalSystem::compile(const Reaction* reaction) {
		compiled_reactions[reaction->getIndex()] = {reaction};
		CompiledReaction& compiled_reaction = compiled_reactions[reaction->getIndex()];
		for(const auto& reagent : reaction->getReagents()) {
			if(reagent.phase != SOLVENT) {
				const auto& component = components_by_name.find(reagent.name);
				if(component != components_by_name.end()) {
					compiled_reaction.components.emplace_back(
							reagent.coefficient,
							// Access in components since components_by_name
							// gives only a const access to components
							components[component->second->getIndex()].get()
							);
				} else {
					// Access in species since species_by_name gives only a
					// const access to species
					auto& species = this->species[
						species_by_name.find(reagent.name)->second->getIndex()
					];
					compiled_reaction.produced_species = {reagent.coefficient, species.get()};
				}
			}
		}
	}

	void ChemicalSystem::initReactionMatrix() {
		for(auto& reaction : reactions) {
			for(auto& species : reaction->getReagents()) {
				if(components_by_name.find(species.name) == components_by_name.end()) {
					// All components are already added to the system using
					// addCompoent(). All other reagents are considered as
					// compound species.
					if(species_by_name.find(species.name) == species_by_name.end()) {
						// Adds the species only if it was not added from an
						// other reaction
						addSpecies(species.name, species.phase);
					}
				}
			}
		}

		reaction_matrix.resize(reactions.size());
		for(const auto& reaction : reactions) {
			reaction_matrix[reaction->getIndex()].resize(species.size());
			for(const auto& reactive : reaction->getReagents()) {
				reaction_matrix[reaction->getIndex()][getSpecies(reactive.name).getIndex()]
					= reactive.coefficient;
			}
		}
	}

	const std::unordered_map<std::string, double>& ChemicalSystem::guessInitialExtents() {
		for(const auto& reaction : getReactions()) {
			if(initial_guess_extents.find(reaction->getName())
					== initial_guess_extents.end()) {
				// No manually specified guess extent: automatically compute one

				double reactants_product = 1.0;
				for(const auto& reagent : reaction->getReagents()) {
					if(reagent.coefficient > 0.0) {
						reactants_product *=
							std::pow(getSpecies(reagent.name).activity(), reagent.coefficient);
					}
				}
				if(reactants_product == 0.0) {
					// TODO: support for this situation?
					throw std::logic_error(
							"Some reactants are missing: the problem is ill formed."
							);
				}
			}
		}
		// Concentrations used in the incremental guess process
		std::vector<double> guessed_concentrations(getSpecies().size());
		for(auto& component : getSpecies())
			guessed_concentrations[component->getIndex()]
				= component->concentration();

		// The overall idea of the guess process is that we have an initial
		// quantity of reactives, specified by the user, and some null products.
		// We need to:
		// 1. initialize the products quantity to a non null value.
		// 2. initialize it to a value as close as possible to the final
		// equilibrium so the Newton method can converge.
		//
		// The logK() values are the first hint to find activities. Considering
		// that all free products concentrations are null, reactions
		// with a big positive logK() value will try to consume as much
		// reactive as possible to maximise the products quantity and establish
		// the equilibrium, while reactions with a big negative logK() value
		// will only need to consumme a negligible quantity of reactives to
		// establish the equilibrium.
		//
		// In order to estimate the final equilibrium, we start by ordering
		// reactions by their logK() value, so that reactions that need a lot of
		// reactives will be handled first. Then the algorithm consists in
		// solving each reaction as if they where alone in the system. At each
		// reaction, the estimated quantity of reactives consumed is updated in
		// guessed_concentrations, so that the first reaction with big logK()
		// value have a lot of reactives at the beginning, and the last one can
		// still get close to the equilibrium from what's left.

		std::vector<const Reaction*> sorted_reactions;
		{

			// Sort reactions by logK() values
			for(auto& reaction : getReactions()) {
				double products_product = 1.0;
				for(const auto& reagent : reaction->getReagents()) {
					if(reagent.coefficient < 0.0) {
						products_product *= getSpecies(reagent.name).activity();
					}
				}
				if(products_product == 0.0) {
					// The guess algorithm can only be used when at least one
					// product is missing
					sorted_reactions.push_back(reaction.get());
				}
			}

			std::sort(sorted_reactions.begin(), sorted_reactions.end(),
					[] (const Reaction*& r1, const Reaction*& r2) {
					// > instead of < so that the list is sorted in descending
					// order
					return r1->getLogK() > r2->getLogK();
					});
		}

		for(auto& reaction : sorted_reactions) {
			CHEM_LOG(INFO) <<
				"Try to guess extent for reaction " << reaction->getName() <<
				" (log K = " << reaction->getLogK() << ")";

			// Limiting reactive finding algorithm
			//
			// The limiting reactive is used to define the maximum possible
			// extent of the reaction so that all reagent concentrations stay
			// positive.
			const ChemicalSpecies* limiting_reactive = nullptr;
			double limiting_reactive_coefficient = 0;
			double smallest_limiting_factor = std::numeric_limits<double>::infinity();

			for(auto& reagent : reaction->getReagents()) {
				if(reagent.coefficient > 0.0) {
					// Consider only reactives
					const auto& component = components_by_name.find(reagent.name);
					if(
							// The reagent is not a component
							component == components_by_name.end()
							// Or if it is a component, it is not fixed
							|| !(component->second->isFixed())) {
						const ChemicalSpecies& species = getSpecies(reagent.name);
						// A fixed reactive cannot be limiting
						double limiting_factor = species.quantity(
								guessed_concentrations[species.getIndex()]
								) / reagent.coefficient;
						if(limiting_factor < smallest_limiting_factor) {
							limiting_reactive = &species;
							limiting_reactive_coefficient = reagent.coefficient;
							smallest_limiting_factor = limiting_factor;
						}
					}
				} 
			}

			// See the GuessF definition for more details.
			// Return log(Q(x)) - logK where x is the extent of the reaction and
			// Q is the reaction quotient with products at the numerator
			GuessF f = GuessF(*this, *reaction, guessed_concentrations);

			double max_N;
			if(limiting_reactive != nullptr) {
				CHEM_LOG(DEBUG) << "  Limiting reactive: " << limiting_reactive_coefficient
					<< " " << limiting_reactive->getName() << " " << 
					limiting_reactive_coefficient * smallest_limiting_factor << " mol/l";

				max_N = smallest_limiting_factor;
			} else {
				// In this case, there is no limiting factor. This for example
				// the case for the reaction OH-: H2O <-> H+ + OH- where H2O is
				// a solvent and H+ is fixed with the pH. Using a value of max_N
				// = 1mol/l still allows to find the corresponding OH-
				// concentration using the RegulaFalsi solver.
				max_N = 1*mol/l;
			}

			// Uses the RegulaFalsi method to find the root of f(x)=log(Q(x)) - logK
			// Solutions are searched in ]-N, 0[ where N is the limiting factor.
			// Due to the convention that products are specified with negative
			// coefficients and that we find a solution that necessarily form
			// products, extents considered are negative.
			double guess_extent = RegulaFalsi<double>(
					// Value just below 0
					- std::numeric_limits<double>::min(),
					// Value just above -N
					-(max_N) * (1 - std::numeric_limits<double>::epsilon()),
					// Function to optimize
					f
					).solve_iter(10000);
			// About std::numeric_limits<double>::epsilon():
			// epsilon is the next float that can be represented after 1.0, so
			// that it is the smallest quantity such that:
			// 1.0 - epsilon < 1.0 < 1.0 + epsilon
			// Multiplying this inequality by N gives the "just above -N" value.

			CHEM_LOG(INFO) << "  Guessed extent: " << guess_extent;

			// Set up the guessed value
			initial_guess_extents[reaction->getName()] = guess_extent;
			// Updates concentrations from the guessed extent, so that the
			// system remains consistent for next reactions.
			for(auto& reagent : reaction->getReagents()) {
				auto& species = getSpecies(reagent.name);
				guessed_concentrations[species.getIndex()]
					= species.concentration(
							guessed_concentrations[species.getIndex()],
							reagent.coefficient * guess_extent
							);
			}
		}
		return initial_guess_extents;
	}

	void ChemicalSystem::proceed(const Reaction& reaction, double extent) {
		CompiledReaction& compiled_reaction
			= compiled_reactions[reaction.getIndex()];
		for(auto& component : compiled_reaction.components) {
			component.component->getSpecies()->incrementConcentration(
					component.coefficient * extent 
					);
		}
		compiled_reaction.produced_species.species->incrementConcentration(
				compiled_reaction.produced_species.coefficient * extent
				);
	}

	void ChemicalSystem::solveEquilibrium() {
		setUp();
		guessInitialExtents();

		using namespace solver;
		solver::X dx = solve(*this);
		for(auto& reaction : reactions) {
			proceed(*reaction.get(), dx[reaction->getIndex()]);
		}
	}

	double ChemicalSystem::reactionQuotient(const std::string& name) const {
		double quotient = 1;
		auto& reaction = getReaction(name);
		for(auto& reactive : reaction.getReagents()) {
			quotient *= std::pow(getSpecies(reactive.name).activity(), reactive.coefficient);
		}
		// 1/quotient due to the following conventions:
		// - products are given with **negative** stoichiometric coefficients
		// - products from the numerous of the reaction quotient
		return 1/quotient;
	}

	void ChemicalSystem::massConservationLaw(std::vector<double>& x) const {
		for(auto& component : getComponents()) {
			x[component->getIndex()]
				= component->getSpecies()->quantity();
		}
		for(auto& reaction : compiled_reactions) {
			for(const ComponentReagent& reagent
					: reaction.components) {
				x[reagent.component->getIndex()] +=
					reagent.coefficient / (-reaction.produced_species.coefficient)
					* reaction.produced_species.species->quantity();
			}
		}
	}
}

#include "chemmisol/chemical/system.h"
#include "chemmisol/math/regula_falsi.h"

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cassert>

namespace chemmisol {

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
		switch(phase) {
			case AQUEOUS:
			case MINERAL:
				{
					ChemicalSpecies* species = nullptr;
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
						default:
							// Cannot occur
							break;
					}
					addComponent(
							new Component(species, component_index, species->quantity()),
							species_index, component_index
							);
				}
				break;
			case SOLVENT:
				Solvent* solvent = new Solvent(name, species_index);
				addComponent(
						new FixedComponent(
							solvent, component_index, solvent->ChemicalSpecies::quantity()
							),
						species_index, component_index
						);
		}
	}

	void ChemicalSystem::addSolvent(const std::string& name) {
		addComponent(name, SOLVENT, 0.0); // Concentration is ignored when
										  // phase=SOLVENT
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

	const std::vector<ComponentReagent>& ChemicalSystem::getComponentReagents(
			const Reaction& reaction) const {
		return compiled_reactions[reaction.getIndex()].components;
	}
	const ChemicalSpeciesReagent& ChemicalSystem::getSpeciesReagent(
			const Reaction& reaction) const {
		return compiled_reactions[reaction.getIndex()].produced_species;
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

	double ChemicalSystem::distanceToEquilibrium(
			const std::vector<double>& activities,
			const Reaction& reaction) const {
		double reactives = reaction.getK();
		double products = 1.0;
		auto& compiled_reaction = compiled_reactions[reaction.getIndex()];
		for(const auto& reagent : compiled_reaction.components) {
			double activity = activities[reagent.component->getSpecies()->getIndex()];
			if(activity > 0.0) {
				if(reagent.coefficient > 0) {
					reactives *= std::pow(activity, reagent.coefficient);
				} else {
					products *= std::pow(activity, -reagent.coefficient);
				}
			}
		}
		double activity = activities[compiled_reaction.produced_species.species->getIndex()];
		if(activity > 0.0)
			// This coefficient is necessarily negative
			products *= std::pow(
					activity,
					-compiled_reaction.produced_species.coefficient);
		return products - reactives;
	}

	double ChemicalSystem::distanceToEquilibrium(const Reaction& reaction) const {
		std::vector<double> activities(species.size());
		for(const auto& species : this->species)
			activities[species->getIndex()] = species->activity();
		return distanceToEquilibrium(activities, reaction);
	}

	void ChemicalSystem::solveEquilibrium() {
		setUp();

		using namespace solver;
		solver::X activities = solve(*this);
		CHEM_LOG(TRACE) << "Solved activities: " << activities;
		for(std::size_t index = 0; index < activities.size(); index++) {
			species[index]->setActivity(activities[index]);
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

	void ChemicalSystem::massConservationLaw(
			const std::vector<double>& activities,
			std::vector<double>& result) const {
		for(auto& component : getComponents()) {
			result[component->getIndex()]
				= activities[component->getIndex()];
		}
		for(auto& reaction : compiled_reactions) {
			for(const ComponentReagent& reagent
					: reaction.components) {
				if(!reagent.component->isFixed())
					result[reagent.component->getIndex()] +=
						reagent.coefficient / (-reaction.produced_species.coefficient)
						* activities[reaction.produced_species.species->getIndex()];
			}
		}

	}

	void ChemicalSystem::massConservationLaw(std::vector<double>& result) const {
		std::vector<double> activities(species.size());
		for(const auto& species : this->species)
			activities[species->getIndex()] = species->activity();
		massConservationLaw(activities, result);
	}
}

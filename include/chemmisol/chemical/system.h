#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>
#include <memory>
#include "chemmisol/math/linear.h"
#include "solver.h"
#include "../logging.h"

namespace chemmisol {
	struct ComponentReagent {
		double coefficient;
		Component* component;

		ComponentReagent() = default;
		ComponentReagent(
				double coefficient, Component* component)
			: coefficient(coefficient), component(component) {
			}
	};

	struct ChemicalSpeciesReagent {
		double coefficient;
		ChemicalSpecies* species;

		ChemicalSpeciesReagent() = default;
		ChemicalSpeciesReagent(
				double coefficient, ChemicalSpecies* species)
			: coefficient(coefficient), species(species) {
			}
	};

	/**
	 * A ChemicalSystem is defined by a set of Components that interact
	 * according to defined Reactions.
	 *
	 * A ChemicalSystem can be used to define both a pure solution system (where
	 * all species are aqueous) as well as a mineral adsorption model. See
	 * constructors for more information.
	 */
	class ChemicalSystem {
		private:
			struct CompiledReaction {
				const Reaction* reaction;
				ChemicalSpeciesReagent produced_species;
				std::vector<ComponentReagent> components;

				CompiledReaction() = default;
				CompiledReaction(const Reaction* reaction)
					: reaction(reaction) {
					}
			};

			std::size_t component_index = 0;
			std::size_t species_index = 0;
			std::size_t reaction_index = 0;

			std::vector<std::unique_ptr<ChemicalSpecies>> species;
			std::unordered_map<std::string, const ChemicalSpecies*> species_by_name;
			std::vector<std::unique_ptr<Component>> components;
			std::unordered_map<std::string, const Component*> components_by_name;
			std::vector<std::unique_ptr<Reaction>> reactions;
			std::unordered_map<std::string, const Reaction*> reactions_by_name;
			std::vector<std::vector<double>> reaction_matrix;

			std::vector<CompiledReaction> compiled_reactions;

			void addSpecies(ChemicalSpecies* component, std::size_t index);
			void addComponent(
					Component* component,
					std::size_t species_index,
					std::size_t component_index);
			void addReaction(Reaction* reaction, std::size_t index);

			std::size_t max_iteration = 200;

			// Adsorption model parameters
			double solid_concentration;
			double specific_surface_area;
			double site_concentration;

			void addComponent(
					const std::string& name,
					Phase phase,
					double concentration,
					std::size_t species_index,
					std::size_t component_index
					);
			void fixComponent(
					const std::string& name,
					Phase phase,
					double concentration,
					std::size_t species_index,
					std::size_t component_index
					);

			void addSpecies(
					const std::string& name,
					Phase phase,
					double concentration,
					std::size_t index
					);

			void addReaction(
					std::string name,
					double K,
					std::vector<Reagent> reactives,
					std::size_t index);

			void compile(const Reaction* reaction);

			/**
			 * Initializes the reaction matrix from all components and reactions
			 * added to the system until now.
			 *
			 * Automatically called by solveEquilibrium(), but might be called
			 * by the user.
			 */
			void initReactionMatrix();

		protected:
			void addSpecies(
					const std::string& name
					);

			void addSpecies(
					const std::string& name,
					Phase phase
					);

		public:
			ChemicalSystem(const ChemicalSystem& other);

			/**
			 * Defines a pure solution model, with only AQUEOUS components.
			 *
			 * Trying to force an addComponent() call with MINERAL phase will
			 * yield undefined behaviors.
			 */
			ChemicalSystem() = default;

			/**
			 * Defines an adsorption model that can handle both AQUEOUS and
			 * MINERAL species. See MineralComponent for more detailed
			 * information about how mineral parameters are defined and used.
			 *
			 * The surface component corresponds to the "free sites" species. It
			 * is automatically added as MineralComponent with a molar fraction
			 * of 1, so that all sites are free at the initial state.
			 *
			 * Chemmisol currently does not support the initialization of
			 * surface components with a non null initial molar fraction.
			 *
			 * @param solid_concentration Quantity of mineral in suspension in
			 * the solution, usually expressed in g/l.
			 * @param specific_surface_area Surface of the solid in contact with
			 * the solution per unit of mass, usually expressed in m2/g.
			 * @param site_concentration Quantity of sites per unit of surface
			 * in contact with the solution, usually expressed as sites/nm2.
			 * @param surface_component Name of the free site surface component
			 * (usually =SOH).
			 */
			// TODO: possibility to initialize surface components with a non
			// null initial molar fraction.
			ChemicalSystem(
					double solid_concentration,
					double specific_surface_area,
					double site_concentration,
					const std::string& surface_component
					);

			void setUp();

			/**
			 * Adds a new Reaction to the system.
			 *
			 * See Reaction::Reaction() for the precise meaning of each
			 * parameter.
			 *
			 * If reagents have not been added with addComponent() before
			 * solveEquilibrium() or initReactionMatrix() is called, missing
			 * components are added to the system with a null initial
			 * concentration.
			 *
			 * @param name Reaction name
			 * @param log_K Equilibrium constant
			 * @param reagents Reagents of the reaction
			 */
			void addReaction(
					std::string name,
					double log_K,
					std::vector<Reagent> reagents
					);

			/**
			 * Adds a new AqueousComponent to the chemical system, assuming that
			 * the Component is AQUEOUS by default.
			 *
			 * @param name Name of the component
			 * @param concentration Initial concentration
			 */
			void addComponent(
					const std::string& name,
					double concentration
					);

			/**
			 * Adds a new Component to the chemical system, depending on the
			 * provided phase, the following Component is instantiated:
			 * - AQUEOUS: AqueousComponent
			 * - MINERAL: MineralComponent
			 *
			 * @param name Name of the component
			 * @param concentration Initial concentration for AQUEOUS species,
			 * initial molar fraction for MINERAL species
			 */
			void addComponent(
					const std::string& name,
					Phase phase,
					double concentration
					);

			void fixComponent(
					const std::string& name,
					double concentration
					);

			void fixComponent(
					const std::string& name,
					Phase phase,
					double concentration
					);

			void addSolvent(const std::string& name);

			/**
			 * Initializes the pH of the chemical system, setting a
			 * concentration of 10^-pH for HO- and H+ species.
			 *
			 * Notice that this only consists in a start value, that is very
			 * likely to vary when the equilibrium is solved, except if fixPH()
			 * was called.
			 *
			 * @param pH Initial pH
			 */
			void initPH(double pH);

			/**
			 * Adds a constraint to the system so that at equilibrium the pH
			 * should be strictly equal to the specified value.
			 *
			 * @param pH Fixed pH
			 */
			void fixPH(double pH);

			double getPH() const;

			/**
			 * Gets the reaction named `name`.
			 *
			 * The behavior of the method is unspecified if the name does not
			 * correspond to any reaction added with addReaction().
			 *
			 * @param name Name of the reaction
			 */
			const Reaction& getReaction(const std::string& name) const;

			/**
			 * Gets the chemical species named `name`.
			 *
			 * The behavior of the method is unspecified if the name does not
			 * correspond to any chemical species added explicitly with
			 * addComponent() (called by the user or initialized by default by
			 * initReactionMatrix()).
			 *
			 * @param name Name of the component
			 */
			const ChemicalSpecies& getSpecies(const std::string& name) const;

			/**
			 * Gets the component with the specified index, that can be
			 * retrieved from an existing component with Component::getIndex().
			 *
			 * The behavior is unspecified if the index does not correspond to
			 * any component added with addComponent() (called by the user or
			 * initialized by default by initReactionMatrix()).
			 *
			 * @param index Index of the component
			 */
			const ChemicalSpecies& getSpecies(const std::size_t& id) const;

			/**
			 * Gets the component named `name`.
			 *
			 * The behavior of the method is unspecified if the name does not
			 * correspond to any component added explicitly with addComponent()
			 * (called by the user or initialized by default by
			 * initReactionMatrix()).
			 *
			 * @param name Name of the component
			 */
			const Component& getComponent(const std::string& name) const;

			/**
			 * Gets the component with the specified index, that can be
			 * retrieved from an existing component with Component::getIndex().
			 *
			 * The behavior is unspecified if the index does not correspond to
			 * any component added with addComponent() (called by the user or
			 * initialized by default by initReactionMatrix()).
			 *
			 * @param index Index of the component
			 */
			const Component& getComponent(const std::size_t& id) const;

			/**
			 * Returns references to all the chemical_species in the system.
			 */
			const std::vector<std::unique_ptr<Component>>& getComponents() const {
				return components;
			}

			/**
			 * Returns references to all the chemical_species in the system.
			 */
			const std::vector<std::unique_ptr<ChemicalSpecies>>& getSpecies() const {
				return species;
			}

			/**
			 * Returns references to all reactions in the system.
			 */
			const std::vector<std::unique_ptr<Reaction>>& getReactions() const {
				return reactions;
			}

			const std::vector<ComponentReagent>& getComponentReagents(
					const Reaction& reaction) const;
			const ChemicalSpeciesReagent& getSpeciesReagent(
					const Reaction& reaction) const;

			/**
			 * Returns the reaction matrix describing the system.
			 *
			 * Each row corresponds to a Reaction of the system, each column
			 * corresponds to a Component of the system, and the coefficient at
			 * (i, j) correspond to the stoichiometric coefficient of the
			 * Component j in the reaction i. Notice that the coefficient might
			 * be null.
			 *
			 * @par Example
			 *
			 * Considering the following chemical system:
			 *
			 *     HO- : H2O       <-> H+ + HO-
			 *     NaCl: Na+ + Cl- <-> NaCl
			 *     NaOH: Na+ + H2O <-> NaOH + H+
			 *
			 * The reaction matrix can be defined as follows:
			 *
			 * |      | H2O | HO- | H+ | NaCl | Na+ | Cl- | NaOH |
			 * |------|-----|-----|----|------|-----|-----|------|
			 * | HO-  | 1   | -1  | -1 | 0    | 0   | 0   | 0    |
			 * | NaCl | 0   | 0   | 0  | -1   | 1   | 1   | 0    |
			 * | NaOH | 1   | 0   | -1 | 0    | 1   | 0   | -1   |
			 *
			 * The actual indexes of each reaction and component correspond to
			 * Reaction::getIndex() and Component::getIndex().
			 */
			const std::vector<std::vector<double>>& getReactionMatrix() const {
				return reaction_matrix;
			}

			void proceed(const Reaction& reaction, double extent);

			double distanceToEquilibrium(
					const std::vector<double>& activities,
					const Reaction& reaction) const;
			double distanceToEquilibrium(const Reaction& reaction) const;

			/**
			 * Solves the equilibrium of the system using the Newton method (see
			 * solver::solve() and solver::F()).
			 *
			 * Concentrations of all components are updated accordingly upon
			 * return.
			 */
			void solveEquilibrium();

			/**
			 * Computes the reaction quotient of the reaction named name.
			 *
			 * By convention, the products of the reaction (i.e. reagents with a
			 * negative coefficient) form the numerous of the quotient.
			 *
			 *
			 * The behavior of the method is unspecified if the name does not
			 * correspond to any component added with addComponent() (called by
			 * the user or initialized by default by initReactionMatrix()).
			 *
			 * @param name Reaction name
			 */
			double reactionQuotient(const std::string& name) const;

			void massConservationLaw(
					const std::vector<double>& activities,
					std::vector<double>& result) const;

			void massConservationLaw(std::vector<double>& result) const;

			/**
			 * Returns the maximum count of iterations allowed for the
			 * equilibrium solver (200 by default).
			 */
			std::size_t getMaxIteration() const {
				return max_iteration;
			}

			/**
			 * Sets the maximum count of iterations allowed for the
			 * equilibrium solver (200 by default).
			 */
			void setMaxIteration(const std::size_t& max_iteration) {
				this->max_iteration = max_iteration;
			}
	};

}

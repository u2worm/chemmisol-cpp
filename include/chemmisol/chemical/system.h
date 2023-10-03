#ifndef CHEMMISOL_SYSTEM_H
#define CHEMMISOL_SYSTEM_H

#include "chemmisol/math/linear.h"
#include "reaction.h"

/**
 * @file chemmisol/chemical/system.h
 *
 * Defines chemical system features.
 */

/**
 * @example basic_chemical_system/main.cpp
 *
 * Basic chemical system definition and equilibrium solver usage example.
 */

namespace chemmisol {
	namespace solver {
		class Solver;
	}

	/**
	 * Compiled reagent corresponding to a ChemicalComponent.
	 *
	 * The ComponentReagent is bound to a concrete ChemicalComponent instance
	 * within a ChemicalSystem.
	 */
	struct ComponentReagent {
		/**
		 * Stoichiometric coefficient associated to the reagent in the reaction.
		 */
		double coefficient;
		/**
		 * Component associated to this reagent within a ChemicalSystem.
		 */
		ChemicalComponent* component;

		ComponentReagent() = default;
		/**
		 * Defines a ComponentReagent.
		 *
		 * @param coefficient Stoichiometric coefficient
		 * @param component Associated component
		 */
		ComponentReagent(
				double coefficient, ChemicalComponent* component)
			: coefficient(coefficient), component(component) {
			}
	};

	/**
	 * Compiled reagent corresponding to a ChemicalSpecies that is not a
	 * component, i.e. the produced species of a reaction.
	 *
	 * The ChemicalSpeciesReagent is bound to a concrete ChemicalSpecies
	 * instance within a ChemicalSystem.
	 */
	struct ChemicalSpeciesReagent {
		/**
		 * Stoichiometric coefficient associated to the reagent in the reaction.
		 */
		double coefficient;
		/**
		 * Component associated to this reagent within a ChemicalSystem.
		 */
		ChemicalSpecies* species;

		ChemicalSpeciesReagent() = default;

		/**
		 * Defines a ChemicalSpeciesReagent.
		 *
		 * @param coefficient Stoichiometric coefficient
		 * @param species Associated chemical species
		 */
		ChemicalSpeciesReagent(
				double coefficient, ChemicalSpecies* species)
			: coefficient(coefficient), species(species) {
			}
	};

	/**
	 * Exception thrown when an invalid species is added to a chemical system.
	 */
	class InvalidSpecies : public std::exception {
		private:
			/**
			 * Reference to system to which an invalid species was added.
			 */
			const ChemicalSystem* chemical_system;
			/**
			 * Name of the invalid species.
			 */
			std::string name;
			/**
			 * Phase of the invalid species.
			 */
			Phase phase;
			
		protected:
			/**
			 * Defines an InvalidSpecies exception.
			 *
			 * @param chemical_system Chemical system to which an invalid
			 * species was added.
			 * @param name Name of the invalid species.
			 * @param phase Phase of the invalid species.
			 */
			InvalidSpecies(
					const ChemicalSystem* chemical_system,
					const std::string& name,
					Phase phase) :
				chemical_system(chemical_system),
				name(name), phase(phase) {
				}

		public:
			/**
			 * Returns a reference to the system to which an invalid component
			 * was added.
			 */
			const ChemicalSystem& getChemicalSystem() const {
				return *chemical_system;
			}

			/**
			 * Name of the invalid species.
			 */
			const std::string& getName() const {
				return name;
			}

			/**
			 * Phase of the invalid species.
			 */
			Phase getPhase() const {
				return phase;
			}
	};

	/**
	 * Exception thrown when trying to add a mineral species to a chemical
	 * system where the mineral sites count is not properly defined.
	 */
	class InvalidMineralSpeciesWithUndefinedSitesCount : public InvalidSpecies {
		private:
			std::string message;

		public:
			/**
			 * Defines an InvalidMineralSpeciesWithUndefinedSitesCount
			 * exception.
			 *
			 * @param chemical_system chemical system with an invalid sites
			 * count.
			 * @param name Name of the mineral species that cannot be added to
			 * the system.
			 */
			InvalidMineralSpeciesWithUndefinedSitesCount(
					const ChemicalSystem* chemical_system,
					const std::string& name
					);

			/**
			 * Returns a message that contains suggestions about how to solve
			 * the issue.
			 */
			const char* what() const noexcept override {
				return message.c_str();
			}
	};

	/**
	 * A ChemicalSystem is defined as a set of Components that interact
	 * according to Reactions in order to produce ChemicalSpecies.
	 *
	 * A ChemicalSystem can be used to define both a pure solution system (where
	 * all species are aqueous) as well as a mineral adsorption model.
	 *
	 * Once components and reactions are specified, the solveEquilibrium()
	 * method can be used to run a solver and set the activities of all chemical
	 * species so that the [law of conservation of
	 * mass](https://en.wikipedia.org/wiki/Conservation_of_mass) and [chemical
	 * equilibriums](https://en.wikipedia.org/wiki/Chemical_equilibrium) are
	 * satisfied.
	 *
	 * All quantities must be specified using the chemmisol unit system, or
	 * directly in core units. See #UNITS().
	 *
	 * @see ChemicalComponent
	 * @see ChemicalSpecies
	 * @see Reaction
	 */
	class ChemicalSystem {
		private:
			/**
			 * Reactions are compiled when setUp() is called. Contrary to raw
			 * Reactions, CompiledReactions are associated to concrete
			 * ChemicalComponent instances and the produced species is formally
			 * identified.
			 */
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

			/* species indexes */
			std::vector<std::unique_ptr<ChemicalSpecies>> species;
			std::unordered_map<std::string, const ChemicalSpecies*> species_by_name;

			/* components indexes */
			std::vector<std::unique_ptr<ChemicalComponent>> components;
			std::unordered_map<std::string, const ChemicalComponent*> components_by_name;

			/* reaction indexes */
			std::vector<std::unique_ptr<Reaction>> reactions;
			std::unordered_map<std::string, const Reaction*> reactions_by_name;

			/* Internal structures initialized by the setUp() method. */
			std::vector<std::vector<double>> reaction_matrix;
			std::vector<CompiledReaction> compiled_reactions;

			void addSpecies(ChemicalSpecies* component, std::size_t index);

			void addComponent(
					ChemicalComponent* component,
					std::size_t species_index,
					std::size_t component_index);

			const Reaction& addReaction(Reaction* reaction, std::size_t index);

			std::size_t max_iteration = 200;

			// Adsorption model parameters
			double solid_concentration = 0;
			double specific_surface_area = 0;
			double site_concentration = 0;

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

			const Reaction& addReaction(
					const std::string& name,
					double K,
					const std::vector<Reagent>& reagents,
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

			template<typename T>
				void _massConservationLaw(
						const std::vector<T>& activities,
						std::vector<T>& result) const;

			template<typename T>
				T _distanceToEquilibrium(
						const std::vector<T>& activities,
						const Reaction& reaction) const;

		public:
			/**
			 * ChemicalSystem copy constructor.
			 */
			ChemicalSystem(const ChemicalSystem& other);

			/*
			 * TODO: ChemicalSystem copy assignment operator.
			 */
			const ChemicalSystem& operator=(const ChemicalSystem& other) = delete;

			/*
			 * TODO: ChemicalSystem move constructor.
			 */
			ChemicalSystem(ChemicalSystem&& other) = delete;

			/*
			 * TODO: ChemicalSystem move assignment operator.
			 */
			const ChemicalSystem& operator=(ChemicalSystem&& other) = delete;

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
			 * The surface complex corresponds to the "free sites" species. It
			 * is automatically added as MineralComponent with a total molar
			 * fraction of 1.
			 *
			 * @param solid_concentration Quantity of mineral in suspension in
			 * the solution, usually expressed in g/l.
			 * @param specific_surface_area Surface of the solid in contact with
			 * the solution per unit of mass, usually expressed in m2/g.
			 * @param site_concentration Quantity of sites per unit of surface
			 * in contact with the solution, usually expressed as entities/nm2.
			 */
			// TODO: possibility to initialize surface components with a non
			// null initial molar fraction.
			ChemicalSystem(
					double solid_concentration,
					double specific_surface_area,
					double site_concentration
					);

			/**
			 * Sets up the system so that it is ready to solve.
			 * 
			 * The method is automatically called by solveEquilibrium(), but can
			 * be called by the user to inspect the state of the system before
			 * solving.
			 *
			 * Reagents of reactions that have not been added with addSpecies(),
			 * addComponent() or fixComponent() are automatically added as
			 * chemical species with an initially null concentration.
			 *
			 * @throws EmptyReagents if a reaction with an empty reagents list
			 * is processed.
			 * @throws MissingProducedSpeciesInReaction if a reaction is missing
			 * a produced species.
			 * @throws TooManyProducedSpeciesInReaction if several produced
			 * species seems to be specified and no unique produced species can
			 * be identified.
			 */
			void setUp();

			/**
			 * Adds a new Reaction to the system.
			 *
			 * See Reaction for the precise meaning of each parameter.
			 *
			 * @param name Reaction name.
			 * @param log_K Equilibrium constant.
			 * @param reagents Reagents of the reaction.
			 * @return Reference to the new Reaction.
			 *
			 * @par Examples
			 * \ref basic_chemical_system/main.cpp
			 */
			const Reaction& addReaction(
					const std::string& name,
					double log_K,
					const std::vector<Reagent>& reagents
					);

			/**
			 * Adds a new AqueousComponent to the chemical system, assuming that
			 * the Component is AQUEOUS by default.
			 *
			 * @param name Name of the component.
			 * @param total_concentration Initial total concentration.
			 *
			 * @par Examples
			 * \ref basic_chemical_system/main.cpp
			 */
			void addComponent(
					const std::string& name,
					double total_concentration
					);

			/**
			 * Adds a new Component to the chemical system, depending on the
			 * provided phase, the following Component is instantiated:
			 * - AQUEOUS: AqueousComponent
			 * - MINERAL: MineralComponent
			 *
			 * @param name Name of the component.
			 * @param phase Chemical phase of the component.
			 * @param total_concentration Initial total concentration for
			 * AQUEOUS species, initial total molar fraction for MINERAL
			 * species.
			 *
			 * @par Examples
			 * \ref basic_chemical_system/main.cpp
			 *
			 * @throws InvalidMineralSpeciesWithUndefinedSitesCount if the phase
			 * is MINERAL but the sites quantity of the system is null. See the
			 * ChemicalSystem(double, double, double, const std::string&) to
			 * ensure the sites quantity is properly defined.
			 */
			void addComponent(
					const std::string& name,
					Phase phase,
					double total_concentration
					);
			/**
			 * Adds a new AqueousComponent to the chemical system, assuming that
			 * the Component is AQUEOUS by default, and fixes it to the
			 * specified concentration.
			 *
			 * @param name Name of the component.
			 * @param total_concentration Fixed total concentration.
			 */
			void fixComponent(
					const std::string& name,
					double total_concentration
					);

			/**
			 * Adds a new Component to the chemical system, depending on the
			 * provided phase, the following Component is instantiated:
			 * - AQUEOUS: AqueousComponent
			 * - MINERAL: MineralComponent
			 * The concentration of the Component is then fixed to the
			 * specified concentration.
			 *
			 * @param name Name of the component.
			 * @param phase Chemical phase of the component.
			 * @param total_concentration Fixed concentration for AQUEOUS species,
			 * fixed molar fraction for MINERAL species.
			 *
			 * @throws InvalidMineralSpeciesWithUndefinedSitesCount if the phase
			 * is MINERAL but the sites quantity of the system is null. See the
			 * ChemicalSystem(double, double, double, const std::string&) to
			 * ensure the sites quantity is properly defined.
			 */
			void fixComponent(
					const std::string& name,
					Phase phase,
					double total_concentration
					);

			/**
			 * Adds a new Solvent to the system.
			 *
			 * @param name Name of the solvent.
			 */
			void addSolvent(const std::string& name);

			/**
			 * Initializes the pH of the chemical system, adding a component
			 * with the specified name and  a concentration of 10^-pH.
			 *
			 * Notice that this only consists in a start value, that is very
			 * likely to vary when the equilibrium is solved, except if fixPH()
			 * was called.
			 *
			 * @param pH Initial pH.
			 * @param h_component_name Name of the component used to compute the
			 * pH.
			 */
			void initPH(double pH, const std::string& h_component_name);

			/**
			 * Initializes the pH in the default H+ component.
			 *
			 * @see initPH(double pH, const std::string& h_component_name)
			 */
			void initPH(double pH) {
				initPH(pH, "H+");
			}

			/**
			 * Fixed the pH of the system to the specified value in the
			 * component with the specified name, that is created if it does not
			 * exist yet.
			 *
			 * @param pH Fixed pH value.
			 * @param h_component_name Name of the component used to compute the
			 * pH.
			 */
			void fixPH(double pH, const std::string& h_component_name);

			/**
			 * Fixes the pH in the default H+ component.
			 *
			 * @see fixPH(double pH, const std::string& h_component_name)
			 */
			void fixPH(double pH) {
				fixPH(pH, "H+");
			};

			/**
			 * Returns the pH of the chemical system from the component with the
			 * specified name.
			 */
			double getPH(const std::string& h_component_name) const;

			/**
			 * Returns the pH of the chemical system from the default H+
			 * component.
			 */
			double getPH() const {
				return getPH("H+");
			}

			/**
			 * Returns the quantity of mineral sites currently defined in the
			 * system. Might be null if the system is aqueous.
			 *
			 * @returns quantity of sites in mol
			 */
			double sitesQuantity() const {
				return MineralSpecies::sites_quantity(
						solid_concentration, specific_surface_area, site_concentration
						);
			}

			/**
			 * Sets the total quantity of the specified component.
			 *
			 * @warning
			 * The quantity of species produced from the specified component,
			 * including the species associated to the component, are not
			 * updated by this method. The solveEquilibrium() method must be run
			 * to update species quantities according the new equilibrium state
			 * and the mass conservation law.
			 */
			void setTotalQuantity(const ChemicalComponent& component, double quantity);

			/**
			 * Sets the total quantity of the specified component to a quantity
			 * that corresponds to the quantity of its associated species at the
			 * provided concentration.
			 *
			 * @warning
			 * The quantity of species produced from the specified component,
			 * including the species associated to the component, are not
			 * updated by this method. The solveEquilibrium() method must be run
			 * to update species quantities according the new equilibrium state
			 * and the mass conservation law.
			 */
			void setTotalConcentration(const ChemicalComponent& component, double concentration);

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
			 * correspond to any chemical species in the system. Valid species
			 * include species added explicitly by \addComponent, \fixComponent
			 * or \addSolvent, and species added implicitly as reagents of
			 * reactions by the setUp() method.
			 *
			 * @param name Name of the component
			 */
			const ChemicalSpecies& getSpecies(const std::string& name) const;

			/**
			 * Gets the component with the specified index, that can be
			 * retrieved from an existing component with Component::getIndex().
			 *
			 * The behavior of the method is unspecified if the index does not
			 * correspond to any chemical species in the system. Valid species
			 * include species added explicitly by \addComponent, \fixComponent
			 * or \addSolvent, and species added implicitly as reagents of
			 * reactions by the setUp() method.
			 *
			 * @param index Index of the component
			 */
			const ChemicalSpecies& getSpecies(const std::size_t& index) const;

			/**
			 * Gets the component named `name`.
			 *
			 * The behavior of the method is unspecified if the name does not
			 * correspond to any component in the system. Valid components
			 * include components added explicitly with \addComponent,
			 * \fixComponent or \addSolvent.
			 *
			 * @param name Name of the component.
			 */
			const ChemicalComponent& getComponent(const std::string& name) const;

			/**
			 * Gets the component with the specified index, that can be
			 * retrieved from an existing component with Component::getIndex().
			 *
			 * The behavior of the method is unspecified if the index does not
			 * correspond to any component in the system. Valid components
			 * include components added explicitly with \addComponent,
			 * \fixComponent or \addSolvent.
			 *
			 * @param index Index of the component.
			 */
			const ChemicalComponent& getComponent(const std::size_t& index) const;

			/**
			 * Checks if the species that corresponds to the specified name is
			 * associated to a component.
			 *
			 * The behavior of the method is unspecified if the name does not
			 * correspond to any chemical species in the system. Valid species
			 * include species added explicitly by \addComponent, \fixComponent
			 * or \addSolvent, and species added implicitly as reagents of
			 * reactions by the setUp() method.
			 *
			 * @param species_name Name of the species to check.
			 * @return True if the species is associated to a component.
			 */
			bool isComponent(const std::string& species_name) const;

			/**
			 * Checks if the species that corresponds to the specified name is
			 * associated to a component.
			 *
			 * @param species Species to check.
			 * @return True if the species is associated to a component.
			 */
			bool isComponent(const ChemicalSpecies& species) const;

			/**
			 * Returns references to all the ChemicalComponents available in the system.
			 */
			const std::vector<std::unique_ptr<ChemicalComponent>>& getComponents() const {
				return components;
			}

			/**
			 * Returns references to all the ChemicalSpecies available in the system.
			 */
			const std::vector<std::unique_ptr<ChemicalSpecies>>& getSpecies() const {
				return species;
			}

			/**
			 * Returns references to all the Reactions added to the system.
			 */
			const std::vector<std::unique_ptr<Reaction>>& getReactions() const {
				return reactions;
			}

			/**
			 * Gets the list of compiled component reagents of the specified
			 * reaction.
			 *
			 * The behavior of the method is unspecified if setUp() has not
			 * been called since the specified reaction was added to the system.
			 *
			 * @param reaction Reference to a reaction added to the system.
			 * @return List of reagents bound to ChemicalComponents of this
			 * chemical system.
			 */
			const std::vector<ComponentReagent>& getComponentReagents(
					const Reaction& reaction) const;

			/**
			 * Gets the compiled produced species of the specified reaction.
			 *
			 * The behavior of the method is unspecified if setUp() has not
			 * been called since the specified reaction was added to the system.
			 *
			 * @param reaction Reference to a reaction added to the system.
			 * @return Produced species of the reaction, bound to a
			 * ChemicalSpecies of this chemical system.
			 */
			const ChemicalSpeciesReagent& getSpeciesReagent(
					const Reaction& reaction) const;

			/**
			 * Returns the reaction matrix describing the system, that should
			 * have been previously computed by the setUp() method.
			 *
			 * Each row corresponds to a Reaction of the system, each column
			 * corresponds to a ChemicalSpecies of the system, and the
			 * coefficient at (i, j) correspond to the stoichiometric
			 * coefficient of the ChemicalSpecies j in the Reaction i. Notice
			 * that the coefficients might be 0.
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
			 * Reaction::getIndex() and ChemicalSpecies::getIndex().
			 *
			 * @return Matrix representing the reactions of the chemical system.
			 */
			const std::vector<std::vector<double>>& getReactionMatrix() const {
				return reaction_matrix;
			}

			/**
			 * Updates the quantities of species in the system resulting from
			 * the specified [extent of the
			 * reaction](https://en.wikipedia.org/wiki/Extent_of_reaction).
			 *
			 * It is the responsibility of the user to ensure that extent is in
			 * a valid range, according to limiting factors of the reaction.
			 *
			 * According to the Chemmisol convention about stoichiometric
			 * coefficients (see Reaction), a positive extent increases the
			 * quantity of **reactants**, and a negative extent increases the
			 * quantity of **products** (since products are associated to
			 * negative coefficients, by convention).
			 *
			 * @param reaction Reaction to proceed.
			 * @param extent Required extent of the the reaction (positive or
			 * negative).
			 */
			void proceed(const Reaction& reaction, double extent);

			/**
			 * Computes the distance between the system state represented by the
			 * specified activities and the equilibrium state of the specified
			 * reaction.
			 *
			 * Considering the following reaction equation (see Reaction for
			 * more details) with an equilibirum constant K:
			 *
			 *     n A + m B - k C - l D <-> 0
			 *
			 * the distance to the equilibrium is computed as:
			 * 
			 *     d = K * a(C)^k * a(D)^l - (a(A)^n + a(B)^m)
			 *
			 * where a(S) denotes the activity of the species S, with a(S)^n set
			 * to 0 if a(S)=0.
			 *
			 * The distance to the equilibrium thus constitutes a polynomial
			 * equation where variables represent the activity of each species,
			 * and the equilibrium of the system constitutes roots of the
			 * polynom. The total degree of this polynom can be queried with the
			 * [degree()](@ref degree(const Reaction&) const) method. In this
			 * case, the total degree would be `max(k+l, n+m)`.
			 *
			 * Such distance might be positive or negative, depending on which
			 * side of the equilibrium the current state is.
			 *
			 * @param activities Activities considered to compute the distance
			 * to equilibrium. Activities must be specified so that
			 * activities[`species`->getIndex()] is equal to the activity of
			 * `species`, for any `species` corresponding to a reagent of the
			 * reaction.
			 * @param reaction Reaction in this chemical system.
			 * @return Distance to the equilibrium state of the specified reaction.
			 */
			double distanceToEquilibrium(
					const std::vector<double>& activities,
					const Reaction& reaction) const;

			std::complex<double> distanceToEquilibrium(
					const std::vector<std::complex<double>>& activities,
					const Reaction& reaction) const;

			/**
			 * Returns the total degree of the polynomial equation representing
			 * the "distance to equilibrium equation".
			 *
			 * See [distanceToEquilibrium()](@ref distanceToEquilibrium(const
			 * std::vector<double>&, const Reaction&) const) for an example.
			 */
			double degree(const Reaction& reaction) const;

			/**
			 * Computes the distance between the current state of the system and
			 * the equilibrium state of the specified reaction.
			 *
			 * See distanceToEquilibrium(const std::vector<double>&, const Reaction&)
			 *
			 * @param reaction Reaction in this chemical system.
			 * @return Distance to the equilibrium state of the specified reaction.
			 */
			double distanceToEquilibrium(const Reaction& reaction) const;

			/**
			 * Solves the equilibrium of the system using the
			 * [default_solver](@ref solver::default_solver).
			 *
			 * Concentrations of all components are updated accordingly upon
			 * return.
			 *
			 * @throws EmptyReagents if a reaction with an empty reagents list
			 * is processed.
			 * @throws MissingProducedSpeciesInReaction if a reaction is missing
			 * a produced species.
			 * @throws TooManyProducedSpeciesInReaction if several produced
			 * species seems to be specified and no unique produced species can
			 * be identified.
			 */
			void solveEquilibrium();

			/**
			 * Solves the equilibrium of the system using the provided solver.
			 *
			 * Concentrations of all components are updated accordingly upon
			 * return.
			 *
			 * @throws EmptyReagents if a reaction with an empty reagents list
			 * is processed.
			 * @throws MissingProducedSpeciesInReaction if a reaction is missing
			 * a produced species.
			 * @throws TooManyProducedSpeciesInReaction if several produced
			 * species seems to be specified and no unique produced species can
			 * be identified.
			 */
			void solveEquilibrium(const solver::Solver& solver);

			/**
			 * Computes the reaction quotient of the specified reaction.
			 *
			 * By convention, the products of the reaction (i.e. reagents with a
			 * **negative** coefficient) form the numerous of the quotient.
			 *
			 * Returns a NaN value if the activity of one of the reactants is
			 * null.
			 *
			 * @param reaction Reaction in this chemical system.
			 * @return Reaction quotient of the reaction.
			 */
			double reactionQuotient(const Reaction& reaction) const;

			/**
			 * Computes the reaction quotient of the reaction corresponding to
			 * the specified name.
			 *
			 * The behavior of the method is unspecified if the name does not
			 * correspond to any reaction in the system. Valid reactions include
			 * reactions added explicitly with \addReaction.
			 *
			 * @see reactionQuotient(const Reaction&)
			 */
			double reactionQuotient(const std::string& name) const;

			/**
			 * Computes the value of the [law of conservation of
			 * mass](https://en.wikipedia.org/wiki/Conservation_of_mass) for a
			 * system state corresponding to the provided species activity.
			 *
			 * The actual total quantity of a component is defined as  the sum of the
			 * quantity of species associated to the component itself and the quantities
			 * of all species compound from this component, weighted by the
			 * stoichiometric coefficient associated to the component in the
			 * reaction producing exactly one unit of the compound. See the
			 * ChemicalComponent documentation for a detailed example.
			 *
			 * @par Example
			 * \parblock
			 *
			 * Let's define the same system as the example of the ChemicalSystem
			 * documentation, with the following components: H+ and PO4.
			 * 
			 * The following reactions can be defined to introduce the OH- and H3PO4
			 * produced species as follows:
			 *
			 *     H2O       <-> OH- + H+
			 *     PO4 + 3H+ <-> H3PO4
			 *
			 * The **actual** (or "output") total quantity N of the PO4 and H+
			 * components are then defined as follows:
			 *
			 *     N(PO4) = n(PO4) + n(H3PO4)
			 *     N(H+) = n(H+) - n(OH-) + 3 * n(H3PO4)
			 *
			 * Then, let's define Q(C) as the **input** total quantity of the
			 * component C, that can be set with
			 * ChemicalComponent::setTotalQuantity(). The law of conservation of
			 * mass states that at any time in the system, the following
			 * relations should hold:
			 *
			 *     Q(PO4) = N(PO4)
			 *     Q(H+)  = N(H+)
			 *
			 * and so:
			 *
			 *     n(PO4) + n(H3PO4) - Q(PO4) = 0
			 *     n(H+) - n(OH-) + 3 * n(H3PO4) - Q(H+) = 0
			 *
			 * Those constraints are thus added to the equation system solved by
			 * the solveEquilibrium() method, next to the chemical equilibrium
			 * constraints for each reaction.
			 *
			 * The massConservationLaw() method then returns a vector such that
			 * `result[i] = N(C[i])-Q(C[i])` where `C[i]` represents the
			 * component for index i, except if the component `C[i]` is
			 * **fixed**. In this case, result[i] = 0 (see the note below).
			 *
			 * \endparblock
			 *
			 * @note
			 * \parblock
			 *
			 * The usage of the mass conservation law is only relevant for
			 * components that are not fixed. Indeed, the notion of "fixed
			 * component" cannot exist in a closed real system, since the extent
			 * of any reaction will modify the quantity of fixed species,
			 * including solvents. In order to apply the law of conservation of
			 * mass to a fixed component, it is necessary to implicitly define a
			 * "buffer" (for example, a [pH buffer
			 * solution](https://en.wikipedia.org/wiki/Buffer_solution)) that
			 * can consume and produce any desired quantity of the species
			 * associated to the fixed component to keep the quantity of the
			 * fixed species to the desired value . The law of mass conservation
			 * is then rather expressed as a [mass
			 * balance](https://en.wikipedia.org/wiki/Mass_balance) equation of
			 * the form `Input + Generation = Output + Accumulation +
			 * Consumption`. For example, the mass conservation of H2O in the
			 * above system should then be expressed as:
			 *
			 *     N(H2O) + H2O production = n(H2O) + n(OH-) + H2O consumption
			 *
			 * so that _H2O production_ and _H2O consumption_ are set to any value that
			 * guarantees that n(H2O) is fixed to the desired value, whatever the
			 * quantity of OH- (and so the extent of the H20 <-> H+ + OH- reaction) is.
			 *
			 * Rather than trying to simulate the action of a buffer, using a
			 * fixed component implies the existence of a perfect buffer, what is
			 * mathematically represented by the fact that in any case
			 * `Q(C)-N(C)=0` for any fixed species. In consequence, no
			 * constraint is added to the solved system for fixed species and
			 * `result[index]=0` for any index corresponding to a fixed
			 * component.
			 *
			 * \endparblock
			 *
			 * @param activities Activities of all species in the system, so
			 * that `activities[i]` corresponds to the activity of the species
			 * with index i.
			 * @param result Vector in which massConservationLaw results for all
			 * components in the system are written, so that `result[i]`
			 * corresponds to the result of the mass conservation law for the
			 * component with index i.
			 */
			void massConservationLaw(
					const std::vector<double>& activities,
					std::vector<double>& result) const;

			void massConservationLaw(
					const std::vector<std::complex<double>>& activities,
					std::vector<std::complex<double>>& result) const;

			/**
			 * Computes the value of the [law of conservation of
			 * mass](https://en.wikipedia.org/wiki/Conservation_of_mass) for the
			 * current system state.
			 *
			 * See [massConservationLaw(activities, result)](@ref massConservationLaw(const std::vector<double>&, std::vector<double>&) const)
			 * for detailed information.
			 *
			 * @param result Vector in which massConservationLaw results for all
			 * components in the system are written, so that `result[i]`
			 * corresponds to the result of the mass conservation law for the
			 * component with index i.
			 */
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

	template<typename T>
		void ChemicalSystem::_massConservationLaw(
				const std::vector<T>& activities,
				std::vector<T>& result) const {
			for(auto& component : getComponents()) {
				if(component->isFixed()) {
					result[component->getIndex()] = 0.0;
				} else {
					result[component->getIndex()]
						= component->getSpecies()->quantity(
								activities[component->getSpecies()->getIndex()]
								);
				}
			}
			for(auto& reaction : compiled_reactions) {
				for(const ComponentReagent& reagent
						: reaction.components) {
					if(!reagent.component->isFixed()) {
						result[reagent.component->getIndex()] +=
							reagent.coefficient / (-reaction.produced_species.coefficient)
							* reaction.produced_species.species->quantity(
									activities[reaction.produced_species.species->getIndex()]
									);
					}
				}
			}
			for(auto& component : getComponents()) {
				if(!component->isFixed()) {
					result[component->getIndex()]
						-= component->getTotalQuantity();
				}
			}
		}

	template<typename T>
		T ChemicalSystem::_distanceToEquilibrium(
				const std::vector<T>& activities,
				const Reaction& reaction) const {
			T reactives = {reaction.getK()};
			T products {1.0};
			auto& compiled_reaction = compiled_reactions[reaction.getIndex()];
			for(const auto& reagent : compiled_reaction.components) {
				T activity = activities[reagent.component->getSpecies()->getIndex()];
				if(activity != 0.0) {
					if(reagent.coefficient > 0) {
						reactives *= std::pow(activity, reagent.coefficient);
					} else {
						products *= std::pow(activity, -reagent.coefficient);
					}
				} else {
					if(reagent.coefficient > 0) {
						reactives = 0.0;
					} else {
						products = 0.0;
					}
				}
			}
			T activity = activities[compiled_reaction.produced_species.species->getIndex()];
			if(activity != 0.0)
				// This coefficient is necessarily negative
				products *= std::pow(
						activity,
						-compiled_reaction.produced_species.coefficient);
			else
				products = 0.0;
			return products - reactives;
		}
}
#endif /*CHEMMISOL_SYSTEM_H*/

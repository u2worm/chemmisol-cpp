#ifndef CHEMMISOL_REACTION_H
#define CHEMMISOL_REACTION_H

#include "species.h"
#include <cmath>

/**
 * @file chemmisol/chemical/reaction.h
 *
 * Chemical reactions related features.
 */

namespace chemmisol {

	/**
	 * Describes a reagent of a reaction.
	 */
	struct Reagent {
		/**
		 * Name of the reagent.
		 */
		std::string name;
		/**
		 * Phase of the reagent.
		 */
		Phase phase;
		/**
		 * Stoichiometric coefficient associated to the reagent in the described
		 * reaction.
		 */
		double coefficient;

		Reagent() = default;
		/**
		 * Defines a Reagent.
		 *
		 * @param name Name of the reagent.
		 * @param phase Phase of the reagent.
		 * @param coefficient Stoichiometric coefficient associated to the
		 * reagent in the described reaction.
		 */
		Reagent(
				const std::string& name,
				Phase phase,
				double coefficient)
			: name(name), phase(phase), coefficient(coefficient) {
			}

		/**
		 * Defines a Reagent.
		 *
		 * The phase is set to AQUEOUS by default, or to SOLVENT if the
		 * reagent name is H2O.
		 *
		 * @param name Name of the reagent.
		 * @param coefficient Stoichiometric coefficient associated to the
		 * reagent in the described reaction.
		 */
		Reagent(
				const std::string& name,
				double coefficient)
			: name(name), phase(name == "H2O" ? SOLVENT : AQUEOUS), coefficient(coefficient) {
			}
	};

	class ChemicalSystem;
	class Reaction;

	/**
	 * Exception thrown when a ChemicalSystem tries to process an ill formed
	 * reaction.
	 *
	 * A reaction is considered valid when its list of reagents is not empty and
	 * contains only components except one species (the "produced species") that
	 * cannot be a component.
	 */
	class InvalidReaction : public std::exception {
		protected:
			/**
			 * Reference to the chemical system that tried to handle the ill
			 * formed reaction.
			 */
			const ChemicalSystem* chemical_system;
			/**
			 * Reference to the ill formed reaction.
			 */
			const Reaction* invalid_reaction;

		public:
			/**
			 * Defines an InvalidReaction exception.
			 *
			 * @param chemical_system chemical system that tried to handle the
			 * ill formed reaction.
			 * @param invalid_reaction ill formed reaction.
			 */
			InvalidReaction(
					const ChemicalSystem* chemical_system,
					const Reaction* invalid_reaction)
				: chemical_system(chemical_system), invalid_reaction(invalid_reaction) {
				}

			/**
			 * Returns a reference to the chemical system that tried to handle
			 * the ill formed reaction.
			 */
			const ChemicalSystem& getChemicalSystem() const {
				return *chemical_system;
			}

			/**
			 * Returns reference to the ill formed reaction.
			 */
			const Reaction& getInvalidReaction() const {
				return *invalid_reaction;
			}
	};

	/**
	 * Exception thrown when the reagents list of a reaction is empty.
	 */
	class EmptyReagents : public InvalidReaction {
		private:
			std::string message;
		public:
			/**
			 * Defines an EmptyReagents exception.
			 *
			 * @param chemical_system chemical system that tried to handle the
			 * ill formed reaction.
			 * @param invalid_reaction ill formed reaction.
			 */
			EmptyReagents(
					const ChemicalSystem* chemical_system,
					const Reaction* invalid_reaction);

			/**
			 * Returns a message that contains suggestions about how to solve
			 * the issue, such as components defined in the chemical system.
			 */
			const char* what() const noexcept override {
				return message.c_str();
			}
	};

	/**
	 * Exception thrown when the produced species seems to be missing in a
	 * reaction definition.
	 */
	class MissingProducedSpeciesInReaction : public InvalidReaction {
		private:
			std::string message;
		public:
			/**
			 * Defines an MissingProducedSpeciesInReaction exception.
			 *
			 * @param chemical_system chemical system that tried to handle the
			 * ill formed reaction.
			 * @param invalid_reaction ill formed reaction.
			 */
			MissingProducedSpeciesInReaction(
					const ChemicalSystem* chemical_system,
					const Reaction* invalid_reaction);

			/**
			 * Returns a message that contains suggestions about how to solve
			 * the issue, such as current reaction reagents and components
			 * defined in the chemical system.
			 */
			const char* what() const noexcept override {
				return message.c_str();
			}
	};

	/**
	 * Exception thrown when more than one produced species (i.e. a species that
	 * does not correspond to a component) is found in the definition of the
	 * reaction.
	 *
	 * The reaction system should either be rewritten, or components should be
	 * defined in the ChemicalSystem.
	 */
	class TooManyProducedSpeciesInReaction : public InvalidReaction {
		private:
			std::string message;

		public:
			/**
			 * Defines an TooManyProducedSpeciesInReaction exception.
			 *
			 * @param chemical_system chemical system that tried to handle the
			 * ill formed reaction.
			 * @param invalid_reaction ill formed reaction.
			 */
			TooManyProducedSpeciesInReaction(
					const ChemicalSystem* chemical_system,
					const Reaction* invalid_reaction);
			/**
			 * Returns a message that contains suggestions about how to solve
			 * the issue, such as current reaction reagents and components
			 * defined in the chemical system.
			 */
			const char* what() const noexcept override {
				return message.c_str();
			}
	};

	/**
	 * A chemical reaction is a process that transform _reactants_ into
	 * _products_.
	 *
	 * A reaction can be fully described by:
	 * - a list of reactives and products (reagents)
	 * - stoichiometric coefficients
	 * - a log K value
	 *
	 * Let's consider the following example reaction:
	 *
	 *     n A + m B <-> k C + l D
	 *
	 * A and B are reactants, C and D are products and n, m, k and l are
	 * corresponding stoichiometric coefficients.
	 *
	 * By convention, the reaction is rewritten as
	 *
	 *     n A + m B - k C - l D <-> 0
	 *
	 * so that stoichiometric coefficients of reactants are **positive** and
	 * stoichiometric coefficients of products are **negative**.
	 *
	 * The [equilibrium
	 * constant](https://en.wikipedia.org/wiki/Chemical_equilibrium) K is then
	 * defined according to the [law of mass
	 * action](https://en.wikipedia.org/wiki/Law_of_mass_action) as:
	 *
	 *    ([C]^k * [D]^l) / ([A]^n * [B]^m) = K
	 *
	 * where the bracket notation denotes the activity of each component **at
	 * equilibrium**. By convention, the products form the numerator.
	 *
	 * This should follow the conventions used by the reference
	 * [VMinteq](https://vminteq.com/) software, so that values used in the
	 * VMinteq database can be reused as is in Chemmisol.
	 *
	 * Within this same software, the name of each reaction also defines its
	 * produced species, to which a coefficient of -1 is associated. This
	 * convention comes from the idea that a species is actually a compound of
	 * base species ("components"), that corresponds to the reactants of the
	 * reaction that form the produced species. For example, in VMinteq, to
	 * define the Na+ + Cl- <-> NaCl reaction, you only need to create a
	 * reaction called "NaCl" with two components, Na+ and Cl-, both with a
	 * coefficient of +1, in the idea that NaCl is compound by the basic Na+ and
	 * Cl- ions.
	 *
	 * In Chemmisol the produced species must be specified explicitly as a
	 * reagent, and the name of each reaction is arbitraty, even if we recommend
	 * to stick with the VMinteq convention to name reactions. The name of each
	 * reaction must also be unique within a ChemicalSystem.
	 *
	 * @par Examples
	 *
	 *     // H2O <-> OH- + H+
	 *     Reaction oh("OH-",
	 *         -13.997 // log K value, with H+ and OH- as products
	 *         { // Reactives
	 *             {"OH-", -1}, // Product
	 *             {"H+", -1},  // Product
	 *             {"H2O", 1}   // Reactant
	 *         });
	 *
	 *     // Na+ + Cl- <-> NaCl
	 *     Reaction nacl("NaCl",
	 *         -0.3 // log K value, with NaCl as product
	 *         { // Reactives
	 *             {"NaCl", -1}, // Product
	 *             {"Na+", -1},  // Reactant
	 *             {"Cl-", 1}    // Reactant
	 *         });
	 *
	 * @see ChemicalSystem::addComponent(const std::string &, Phase, double)
	 * @see ChemicalSystem::addReaction(std::string, double, std::vector< Reagent >)
	 */
	class Reaction {
		private:
			std::string name;
			std::size_t index;
			double K;
			double log_K;
			std::vector<Reagent> reagents;

		public:
			/**
			 * Defines a new reaction.
			 *
			 * @param name Name of the reaction
			 * @param index Index used to retrieve the reaction in data
			 * structures used internally by the ChemicalSystem. The index can
			 * also be used to uniquely identify each reaction of a system.
			 * @param log_K Base 10 logarithm of the equilibrium constant. By
			 * convention, K should correspond to the reaction quotient at
			 * equilibrium with products (i.e. species with **negative**
			 * coefficients) at the numerous. For example, using this
			 * convention, the [self-ionization of water
			 * reaction](https://en.wikipedia.org/wiki/Self-ionization_of_water)
			 * is described with coefficients of -1 for HO- and H+ species, a
			 * coefficient of +1 for H2O, and a log_K value of -14.
			 * @param reagents Defines reagents of the reaction and associates
			 * a stoichiometric coefficient to each, where coefficients of
			 * products are negative and coefficient of reactants are positive
			 * by convention. Reagents are identified by a component name, that
			 * should be a valid argument for the ChemicalSystem::getComponent()
			 * method.
			 */
			Reaction(
					const std::string& name, std::size_t index, double log_K,
					const std::vector<Reagent>& reagents
					)
				: name(name), index(index), K(std::pow(10, log_K)), log_K(log_K), reagents(reagents) {
				}

			/**
			 * Name of the reaction, usually the name of its produced species by
			 * convention.
			 */
			const std::string& getName() const {
				return name;
			}

			/**
			 * Index used to retrieve the reaction in data structures used
			 * internally by the ChemicalSystem. The index can also be used to
			 * uniquely identify each reaction of a system.
			 */
			std::size_t getIndex() const {
				return index;
			}

			/**
			 * Equilibrium constant of this reaction.
			 */
			double getK() const {
				return K;
			}

			/**
			 * Base 10 logarithm of the equilibrium constant of this reaction.
			 */
			double getLogK() const {
				return log_K;
			}

			/**
			 * Returns the reagents of the reaction.
			 */
			const std::vector<Reagent>& getReagents() const {
				return reagents;
			}
	};

}
#endif /*CHEMMISOL_REACTION_H*/

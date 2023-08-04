#include <string>
#include "species.h"
#include <cmath>

namespace chemmisol {

	struct Reagent {
		std::string name;
		Phase phase;
		double coefficient;

		Reagent() = default;
		Reagent(
				const std::string& name,
				Phase phase,
				double coefficient)
			: name(name), phase(phase), coefficient(coefficient) {
			}
		Reagent(
				const std::string& name,
				double coefficient)
			: name(name), phase(name == "H2O" ? SOLVENT : AQUEOUS), coefficient(coefficient) {
			}
	};

	bool operator==(const Reagent& c1, const Reagent& c2);

	/**
	 * A chemical reaction is a process that transform _reactants_ into
	 * _products_.
	 *
	 * A reaction can be fully described by:
	 * - a list of reactives and products
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
	 * [VMinteq](https://vminteq.com/) software, so that values observed in the
	 * VMinteq database can be reused as is in Chemmisol.
	 *
	 * Within this same software, the name of each reaction also defines its
	 * main product, to which a coefficient of -1 is associated. This convention
	 * comes from the idea that a species is actually a compound of base
	 * species, that corresponds to the reactants of the reaction that form the
	 * main product. For example, in VMinteq, to define the Na+ + Cl- <-> NaCl
	 * reaction, you only need to create a reaction called "NaCl" with two
	 * components, Na+ and Cl-, both with a coefficient of +1, in the idea that
	 * NaCl is compound by the basic Na+ and Cl- ions.
	 *
	 * Such conventions might however seem quite awkward for non experts, and
	 * relies on implicit things. That's why in Chemmisol we chose to specify
	 * even the product of the reaction as an explicit product, even if we
	 * recommend to stick with the VMinteq convention to name reactions.
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
			 * @param Index Index used to retrieve the reaction in data
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
			 * @param reactives Defines reagents of the reaction and associates
			 * a stoichiometric coefficient to each, where coefficients of
			 * products are negative and coefficient of reactants are positive
			 * by convention. Reagents are identified by a component name, that
			 * should be a valid argument for the ChemicalSystem::getComponent()
			 * method. Notice that coefficients for unspecified components are
			 * assumed to be null.
			 */
			Reaction(
					const std::string& name, std::size_t index, double log_K,
					const std::vector<Reagent>& reagents
					)
				: name(name), index(index), K(std::pow(10, log_K)), log_K(log_K), reagents(reagents) {
				}

			/**
			 * Name of the reaction, usually the name of its main product by
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

			double getK() const {
				return K;
			}

			/**
			 * Base 10 logarithm of the equilibrium constant.
			 */
			double getLogK() const {
				return log_K;
			}

			/**
			 * Reagents of the reaction.
			 */
			const std::vector<Reagent>& getReagents() const {
				return reagents;
			}
	};

}

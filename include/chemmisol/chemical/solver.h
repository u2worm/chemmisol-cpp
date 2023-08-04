#include "reaction.h"

namespace chemmisol {
	class ChemicalSystem;

	/**
	 * Namespace containing solver features.
	 */
	namespace solver {
		/**
		 * Vector type used by the Newton method.
		 *
		 * Each value corresponds to an **extent** of each reaction of the
		 * system, so that `X[reaction.getIndex()]` corresponds to the extent of
		 * a given `reaction`.
		 */
		typedef std::vector<double> X;
		/**
		 * Matrix type used by the Newton method.
		 *
		 * This is notably used to represent the Jacobian matrix of equilibrium
		 * equations representing the system.
		 */
		typedef std::vector<std::vector<double>> M;

		/**
		 * Function considered by the Newton method to find the
		 * equilibrium state of the ChemicalSystem. See the definition of f()
		 * for detailed explanation on the solved equation system.
		 */
		class F {
			private:
				static const std::size_t INVALID_INDEX;
				const ChemicalSystem& system;
				std::size_t x_size;
				std::size_t f_x_size;
				std::size_t reaction_offset;
				std::vector<std::size_t> components_indexes;
				std::vector<std::size_t> species_offsets;
				std::vector<std::size_t> species_indexes;
				std::vector<std::size_t> revert_species_indexes;
				std::vector<double> fixed_activities;

			public:

				/**
				 * Returns the resulting concentrations in the current chemical
				 * system that result from the provided extents.
				 */
				//X concentrations(const X& extents) const;

				/**
				 * Initializes a function F to find the equilibrium of the
				 * provided chemical system.
				 */
				F(const ChemicalSystem& system);

				X reducedActivities() const;
				X completeActivities(const X& reduced_activities) const;

				const std::vector<std::size_t>& componentsIndexes() const {
					return components_indexes;
				}

				const std::vector<std::size_t>& speciesIndexes() const {
					return species_indexes;
				}

				/**
				 * Returns the value of f for the provided extents.
				 *
				 * @par Example
				 *
				 * Let's consider an example chemical system where two reactions
				 * occur:
				 *
				 *     H2O       <-> H+ + HO-
				 *     Na+ + Cl- <-> NaCl
				 *     Na+ + H2O <-> NaOH + H+
				 *
				 * The [law of mass
				 * action](https://en.wikipedia.org/wiki/Law_of_mass_action)
				 * states that at equilibrium the following relation should
				 * hold:
				 *
				 *     [H+]*[HO-] / [H2O]          = K1
				 *     [Na+]*[Cl-] / [NaCl]        = K2
				 *     [Na+]*[H2O] / ([NaOH]*[H+]) = K3
				 *
				 * where `[H2O]=1` (activity of the solvent) and the brackets
				 * notation denotes the current activity of each component.
				 *
				 * However, when the system is **not** at equilibrium, the
				 * previous relation does **not** hold. In consequence, the
				 * purpose of the problem is to find an _extent_ of the 3
				 * reactions considered (`X=[x1, x2, x3]`) that brings back the
				 * system to the equilibrium. The extent x2 can for example be
				 * interpreted as "when the reaction 1 advances of x2 mole, 1
				 * mole of H2O is consumed to produce 1 mole of H+ and 1 mole of
				 * HO-". The count of each products and components produced for
				 * a unitary extent are called "stoichiometric coefficients".
				 * Notice that extents might be positive or negative.
				 *
				 * The problem then consists in finding x1, x2 and x3 such that:
				 * 
				 *    [H+ + x1 + x2]*[HO- + x1]                     = K1
				 *    [Na+ - x2 - x3]*[Cl- - x2]/[NaCl + x2]        = K2
				 *    [Na+ - x2 - x3]/([NaOH + x3]*[H+ + x1 + x2])  = K3
				 *
				 * Considering the fact that activity are usually exponential
				 * factors with a high amplitude (usually from 10-16 to 10-6),
				 * and to **highly** simplify the derivative of f() (see df()),
				 * the final equation actually solved is as follows:
				 *
				 * Finding `X=[x1, x2, x3]` such that `f(X) = [0, 0, 0]` where
				 *
				 *            log([H+ + x1 + x2]) + log([HO- + x1]) - log(K1)
				 *     f(X) = log([Na+ - x2 - x3]) + log([Cl- - x2]) - log([NaCl + x2]) - log(K2)
				 *            log([Na+ - x2 - x3]) - log([NaOH + x3]) - log([H+ + x1 + x2]) - log(K3)
				 * A few extra information:
				 * - the new concentration of each component considering each
				 *   extent (e.g. H+ + x1 + x2) is computed with the
				 *   concentrations() method, that is itself based on the
				 *   generic Component::concentration() method.
				 * - the activity of each component might be computed very
				 *   differently depending on the type of the component
				 *   (Solvent, AqueousComponent, MineralComponent...). This is
				 *   handled by the corresponding Component::activity()
				 *   implementations.
				 */
				X f(const X& extents) const;

				/**
				 * Computes the [Jacobian
				 * matrix](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant)
				 * of f() at the given extents.
				 *
				 * The computation of each derivative is handled by the
				 * Component::Dactivity() method. Expect for
				 * ElectrostaticComponent that is very particular, the
				 * computation of the matrix is trivial when the system is
				 * expressed in the log form (see f() documentation) since
				 * `dlog(f(x))/dx = 1/x * df(x)/dx * 1/ln(10)` where `f(x)` is
				 * the activity of a given component, that generally corresponds
				 * to a simple sum.
				 *
				 * For example:
				 *
				 *     dlog([H+ + x1 + x2])/dx1 = 1/[H+ + x1 + x2] * 1/V * 1/ln(10)
				 *
				 * since
				 *
				 *     d[H+ + x1 + x2]/dx1 = d ((H+ + x1 + x2)/V)/dx1 = 1/V
				 *
				 * where V is the volume of the solution and the activity of an
				 * aqueous species is defined as its concentration C such that
				 * `C=n/V` where n is the quantity of the species in mole.
				 */
				M df(const X& extents) const;
		};

		/**
		 * Finds and returns reactions extents that will bring back the system
		 * to equilibrium, finding the root of F with the Newton method.
		 */
		X solve(const ChemicalSystem& system);
	}
}

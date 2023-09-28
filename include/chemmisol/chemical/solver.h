#ifndef CHEMMISOL_SOLVER_H
#define CHEMMISOL_SOLVER_H

#include "reaction.h"

/**
 * @file chemmisol/chemical/solver.h
 *
 * Implements chemical equilibrium solver features.
 */

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
		 * Function used in the Newton method to find the equilibrium state of
		 * the ChemicalSystem. See the definition of f() for detailed
		 * explanation on the solved equation system.
		 */
		class F {
			public:
				/**
				 * Special index value used to specify that no entry exist for
				 * some species.
				 */
				static const std::size_t INVALID_INDEX;
			private:
				const ChemicalSystem& system;
				std::size_t x_size;
				std::size_t f_x_size;
				std::size_t reaction_offset;
				std::vector<std::size_t> components_indexes;
				std::vector<std::size_t> species_indexes;
				std::vector<std::size_t> revert_species_indexes;
				std::vector<double> fixed_activities;

			public:

				/**
				 * Initializes a function F to find the equilibrium of the
				 * provided chemical system.
				 *
				 * @param system Chemical system to solve.
				 */
				F(const ChemicalSystem& system);

				/**
				 * Returns the reduced vector of activities.
				 *
				 * The reduced vector corresponds to an initial vector A such
				 * that A[i] equals the activities of species with index i, from
				 * which all entries that corresponds to species associated to
				 * fixed components have been removed.
				 *
				 * In consequence, it is **not** guaranteed that A[i]
				 * corresponds to the activity of the species with index i, but
				 * A[speciesIndexes()[i]] can be used to retrieve the activity
				 * of species i if it's not associated to a fixed component.
				 *
				 * The complete vector of activities can be rebuilt from the
				 * reduced activities using the completeActivities() method.
				 *
				 * @return Vector containing activities of not fixed species in
				 * the solved system.
				 */
				X reducedActivities() const;

				/**
				 * Rebuilds a complete vector activity from a
				 * reducedActivities() vector, so that the built vector A is
				 * such that for any species A[i] equals the activity of species
				 * with index i.
				 *
				 * Activities for species that are not fixed are taken from the
				 * specified reduced_activities vector, and activities for
				 * species associated to fixed components are taken from the
				 * solved chemical system.
				 *
				 * @param reduced_activities Reduced vector of activities,
				 * without entries for species associated to fixed components.
				 * @return Vector containing activities of all species.
				 */
				X completeActivities(const X& reduced_activities) const;

				/**
				 * Maps the indexes of components to indexes in the f() vector.
				 *
				 * If componentsIndexes()[i] is equal to #INVALID_INDEX, no entry
				 * is available for the component with index i in the f() vector
				 * (i.e. the component i is fixed).
				 *
				 * Else, f(a)[speciesIndexes()[i]] returns the result of the
				 * \massConservationLaw for the component with index i
				 * considering the current activities a. 
				 */
				const std::vector<std::size_t>& componentsIndexes() const {
					return components_indexes;
				}

				/**
				 * Offset used to retrieve distance to equilibrium of reactions
				 * in the f() vector, so that f(a)[reactionOffset()+i]
				 * corresponds to the distance to equilibrium for the reaction
				 * with index i.
				 *
				 * @note
				 * By construction, the reaction offset is equal to the count of
				 * law of conservation of mass equations, i.e. the count of not
				 * fixed components.
				 */
				std::size_t reactionOffset() const {
					return reaction_offset;
				}

				/**
				 * Maps the indexes of species to indexes in the
				 * reducedActivities() vector.
				 *
				 * If speciesIndex()[i] is equal to #INVALID_INDEX, no entry is
				 * available for the species with index i in the reduced vector
				 * (i.e. the species i is associated to a fixed component).
				 *
				 * Else, A[speciesIndexes()[i]] returns the activity of species
				 * with index i in the reducedActivity() vector A. 
				 */
				const std::vector<std::size_t>& speciesIndexes() const {
					return species_indexes;
				}

				/**
				 * Returns the value of f for the provided reduced activities.
				 *
				 * The values of f are of the same size as the vector of reduced
				 * activities, and as the following form:
				 *
				 *     m_0
				 *     ...
				 *     m_n
				 *     e_0
				 *     ...
				 *     e_p
				 *
				 * where:
				 * - \parblock
				 *   m_i entries correspond to the law of conservation of mass
				 *   applied to **not fixed** components. This represents n
				 *   equations. The result for the component with index i can be
				 *   retrieved with f(a)[componentsIndexes()[i]].
				 *   @see \massConservationLaw
				 *   \endparblock
				 * - \parblock
				 *   e_i entries correspond to the distance to equilibrium for
				 *   each reaction of the system to solve. This represents m
				 *   equations. The result for the reaction with index i can be
				 *   retrieved with f(a)[reactionOffset()+i].
				 *   @see \distanceToEquilibrium
				 *   \endparblock
				 *
				 *
				 * Considering the definitions of the [law of conservation of
				 * mass](https://en.wikipedia.org/wiki/Conservation_of_mass) and
				 * of the [law of mass
				 * action](https://en.wikipedia.org/wiki/Law_of_mass_action),
				 * finding the equilibrium state of the chemical system
				 * corresponds to find activities such that `f(activities)=0`.
				 *
				 * @param reduced_activities Current activities for species not
				 * associated to a fixed component.
				 * @return Vector value of f, result of the mass conservation
				 * law applied to not fixed components and distances to
				 * equilibrium of each reaction.
				 */
				X f(const X& reduced_activities) const;

				/**
				 * Computes the [Jacobian
				 * matrix](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant)
				 * of f() for the specified reduced_activities.
				 *
				 * Since the expressions of both the \massConservationLaw and
				 * \distanceToEquilibrium are polynomial functions of the
				 * reduced activities, it is straightforward to compute the
				 * partial derivatives of f().
				 */
				M df(const X& extents) const;
		};

		class Solver {
			public:
				virtual X solve(const ChemicalSystem& system) const = 0;

				virtual ~Solver() {
				}
		};

		class AbsoluteNewton : public Solver {
			public:
				/**
				 * Finds and returns activities that correspond to the system to
				 * equilibrium, finding the root of F with the AbsoluteNewton
				 * method.
				 *
				 * @param system Chemical system to solve.
				 * @return A complete vector of activities A, such that A[i]
				 * corresponds to the found activity of species with index i.
				 */
				X solve(const ChemicalSystem& system) const override;
		};

		/**
		 * Default equilibrium solver.
		 */
		extern AbsoluteNewton default_solver;
	}
}
#endif /*CHEMMISOL_SOLVER_H*/

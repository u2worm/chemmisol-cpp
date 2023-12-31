#ifndef CHEMMISOL_EQUILIBRIUM_SOLVER_H
#define CHEMMISOL_EQUILIBRIUM_SOLVER_H

#include "system.h"
#include "chemmisol/math/homotopy.h"
#include <list>
#include <random>
#include <complex>

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
		typedef std::vector<std::complex<double>> CX;
		/**
		 * Matrix type used by the Newton method.
		 *
		 * This is notably used to represent the Jacobian matrix of equilibrium
		 * equations representing the system.
		 */
		typedef std::vector<std::vector<double>> M;
		typedef std::vector<std::vector<std::complex<double>>> CM;

		template<typename M, typename X>
		struct ChemicalLinearSolver {
			std::size_t n_c;
			std::size_t n_s;

			ChemicalLinearSolver(std::size_t n_c, std::size_t n_s)
				: n_c(n_c), n_s(n_s) {
				}

			X solve(const M& df, const X& f) const {
				CHEM_LOGV(9) << "Prepare";
				CHEM_LOGV(9) << "df=" << mview(df);
				CHEM_LOGV(9) << "f=" << xview(f);

				const MView<M> D2 = mview(df, n_c, n_s, n_c, n_s);
				const M inv_D2 = inv_diag(D2);
				const MView<M> inv_D2_view = mview(inv_D2, n_c, n_s, n_c, n_s);
				const MView<M> C = mview(df, 0, n_c, n_c, n_s);
				const MView<M> B = mview(df, n_c, n_s, 0, n_c);
				const MView<M> D1 = mview(df, 0, n_c, 0, n_c);

				const M C_inv_D2 = C * mview(inv_D2, n_c, n_s, n_c, n_s);
				const MView<M> C_inv_D2_view = mview(C_inv_D2, 0, n_c, n_c, n_s);
				const M A = D1 - mview(C_inv_D2_view*B, 0, n_c, 0, n_c);
				const X Y = xview(f, 0, n_c) - xview(
						C_inv_D2_view * xview(f, n_c, n_s), 0, n_c);

				const X x1 = gauss::Gauss<MView<M>, XView<X>>::solve(
						mview(A, 0, n_c, 0, n_c), xview(Y, 0, n_c));

				const X x2 = inv_D2_view*xview(
						xview(f, n_c, n_s)-
						xview(B*xview(x1, 0, n_c), n_c, n_s),
						n_c, n_s
						);

				X x(n_s);
				for(std::size_t i = 0; i < n_c; i++)
					x[i] = x1[i];
				for(std::size_t i = n_c; i < n_s; i++)
					x[i] = x2[i];
				CHEM_LOGV(9) << "x=" << x;
/*
 *                {
 *                    M _df = df;
 *                    X _f = f;
 *
 *                    X _x = gauss::Gauss<M, X>::solve(_df, _f);
 *                    CHEM_LOGV(9) << "Method 3 x=" << _x;
 *                    return _x;
 *                }
 */
				return x;
			}
		};

		template<typename X>
		class ReducedChemicalSystem {
			public:
				/**
				 * Special index value used to specify that no entry exist for
				 * some species.
				 */
				static const std::size_t INVALID_INDEX;
			private:
				const ChemicalSystem& _system;
				std::size_t n_components;
				std::size_t n_species;
				std::size_t x_size;
				std::size_t f_x_size;
				std::size_t reaction_offset;
				std::vector<std::size_t> components_indexes;
				std::vector<std::size_t> species_indexes;
				std::vector<std::size_t> revert_species_indexes;
				X fixed_activities;

			public:
				ReducedChemicalSystem(const ChemicalSystem& system);

				const ChemicalSystem& system() const {
					return _system;
				}

				std::size_t nComponents() const {
					return n_components;
				}

				std::size_t nSpecies() const {
					return n_species;
				}

				std::size_t xSize() const {
					return x_size;
				}

				std::size_t fxSize() const {
					return f_x_size;
				}

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

				const X& fixedActivities() const {
					return fixed_activities;
				}

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


		};

		template<typename X>
		const std::size_t ReducedChemicalSystem<X>::INVALID_INDEX = -1;

		template<typename X>
		ReducedChemicalSystem<X>::ReducedChemicalSystem(const ChemicalSystem& system) :
			_system(system),
			n_components(system.getComponents().size()),
			n_species(system.getSpecies().size()),
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
						--n_components;
						--n_species;
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


		template<typename X>
		X ReducedChemicalSystem<X>::reducedActivities() const {
			X x(x_size);
			for(const auto& species : _system.getSpecies())
				if(species_indexes[species->getIndex()] != INVALID_INDEX)
					x[species_indexes[species->getIndex()]] = species->activity();
			return x;
		}

		template<typename X>
		X ReducedChemicalSystem<X>::completeActivities(const X& reduced_activities) const {
			X complete_activities(_system.getSpecies().size());
			for(const auto& component : _system.getComponents())
				if(component->isFixed())
					complete_activities[component->getSpecies()->getIndex()]
						= component->getSpecies()->activity();
			for(std::size_t i = 0; i < reduced_activities.size(); i++) {
				complete_activities[revert_species_indexes[i]]
					= reduced_activities[i];
			}
			return complete_activities;
		}

		/**
		 * Function used in the Newton method to find the equilibrium state of
		 * the ChemicalSystem. See the definition of f() for detailed
		 * explanation on the solved equation system.
		 */
		template<typename X, typename M>
		class F {
			private:
				typedef typename X::value_type T;
				const ReducedChemicalSystem<X>& reduced_system;
				const ChemicalSystem& system;
			public:
				/**
				 * Initializes a function F to find the equilibrium of the
				 * provided chemical system.
				 *
				 * @param system Chemical system to solve.
				 */
				F(
						const ReducedChemicalSystem<X>& reduced_system,
						const ChemicalSystem& system) :
					reduced_system(reduced_system), system(system) {
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
				M df(const X& reduced_activities) const;
		};

		template<typename X, typename M>
		X F<X, M>::f(const X& activities) const {
			X f_x(reduced_system.fxSize());
			// Complete activities vector with an entry for each species, so
			// that it is compatible with ChemicalSystem::massConservationLaw()
			// and ChemicalSystem::distanceToEquilibrium().
			// TODO: improve this, too many copies for nothing.
			X complete_activities = reduced_system.completeActivities(activities);
			{
				X mass_conservation_results(system.getComponents().size());
				// Actual total quantity of each component
				system.massConservationLaw(
						complete_activities, mass_conservation_results);
				const auto& components_indexes = reduced_system.componentsIndexes();
				for(std::size_t i = 0; i < mass_conservation_results.size(); i++) {
					if(components_indexes[i] != ReducedChemicalSystem<X>::INVALID_INDEX)
						f_x[components_indexes[i]] = mass_conservation_results[i];
				}
			}
			
			for(const auto& reaction : system.getReactions()) {
				f_x[reduced_system.reactionOffset()+reaction->getIndex()]
					= system.distanceToEquilibrium(complete_activities, *reaction);
			}
			return f_x;
		}

		template<typename X, typename M>
		M F<X, M>::df(const X& activities) const {
			M jacobian(reduced_system.fxSize());
			std::vector<std::size_t> component_offsets;

			const auto& components_indexes = reduced_system.componentsIndexes();
			const auto& species_indexes = reduced_system.speciesIndexes();
			// Begin mass conservation law equations
			for(auto& component : system.getComponents()) {
				if(components_indexes[component->getIndex()]
						!= ReducedChemicalSystem<X>::INVALID_INDEX) {
					auto& d_f = jacobian[components_indexes[component->getIndex()]];
					d_f.resize(reduced_system.xSize());
					d_f[species_indexes[component->getSpecies()->getIndex()]]
						= component->getSpecies()->dQdA();
				}
			}
			for(const auto& reaction : system.getReactions()) {
				// Species produced by this reaction
				const ChemicalSpeciesReagent& species
					= system.getSpeciesReagent(*reaction);
				// The index of the species produced by the reaction is
				// necessarily valid.
				for(const auto& reagent : system.getComponentReagents(*reaction)) {
					if(components_indexes[reagent.component->getIndex()]
							!= ReducedChemicalSystem<X>::INVALID_INDEX) {
						jacobian
							[components_indexes[reagent.component->getIndex()]]
							[species_indexes[species.species->getIndex()]]
								= reagent.coefficient/(-species.coefficient)*species.species->dQdA();
					}
				}
			}
			// End mass conservation law equations

			// Begin equilibrium equations
			const auto& fixed_activities = reduced_system.fixedActivities();
			for(const auto& reaction : system.getReactions()) {
				auto& d_f = jacobian[reduced_system.reactionOffset() + reaction->getIndex()];
				d_f.resize(reduced_system.xSize());
				// dx_species: activity variable from which the current reaction
				// d_f is derived
				for(const auto& dx_species : system.getSpecies()) {
					// No derivative to compute for fixed species (that
					// correspond to invalid indexes)
					if(species_indexes[dx_species->getIndex()]
							!= ReducedChemicalSystem<X>::INVALID_INDEX) {
						d_f[species_indexes[dx_species->getIndex()]] = 0.0;
						T d_reactives = -reaction->getK();
						T d_products = 1.0;
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
							if(activities[species_indexes[reagent.species->getIndex()]] != 0.0) {
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
								T species_activity = 0.0;
								if(species_indexes[reagent.component->getSpecies()->getIndex()]
										!= ReducedChemicalSystem<X>::INVALID_INDEX) {
									species_activity = activities[species_indexes[reagent.component->getSpecies()->getIndex()]];
								} else {
									species_activity = fixed_activities[reagent.component->getIndex()];
								}
								if(species_activity != 0.0) {
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
							T activity = activities[index];
							// d x^a/dx = a * x^(a-1)
							if(species_coefficient_in_reactions < 0) {
								d_f[index] = -species_coefficient_in_reactions
									* d_products;
								if(species_coefficient_in_reactions < -1) {
									// Prevents std::pow(0, 0) issues
									d_f[index] *= std::pow(
											activity,
											-species_coefficient_in_reactions-1
											);
								}
							} else {
								d_f[index] = species_coefficient_in_reactions
									* d_reactives;
								if(species_coefficient_in_reactions > 1) {
									// Prevents std::pow(0, 0) issues
									d_f[index] *= std::pow(
											activity,
											species_coefficient_in_reactions-1
											);
								}
							}
						}
					}
				}
			}
			// End equilibrium equations
			return jacobian;
		}

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

		class G {
			private:
				const ReducedChemicalSystem<CX>& reduced_system;
				const ChemicalSystem& system;
				const X degrees;
				const CX a;
				const CX b;


				template<typename R>
					static CX gen_random(R& random_generator, std::size_t n);
				static X compute_degrees(
						const ReducedChemicalSystem<CX>& reduced_system,
						const ChemicalSystem& system
						);

			public:
				template<typename R>
				G(
						const ReducedChemicalSystem<CX>& reduced_system,
						R& random_generator);

				CX g(const CX& reduced_activities) const;
				CM dg(const CX& reduced_activities) const;
				std::vector<CX> initValues() const;
		};

		template<typename R>
			CX G::gen_random(R& random_generator, std::size_t n) {
				std::uniform_real_distribution<double> rd(-1, 1);
				CX x(n);
				for(auto& _x : x)
					_x = {rd(random_generator), rd(random_generator)};
				return x;
			}

		template<typename R>
			G::G(
					const ReducedChemicalSystem<CX>& reduced_system,
					R& random_generator) :
				reduced_system(reduced_system), system(reduced_system.system()),
				degrees(compute_degrees(reduced_system, system)),
				a(gen_random(random_generator, reduced_system.xSize())),
				b(gen_random(random_generator, reduced_system.xSize())) {
				}

		template<typename R>
		class HomotopyContinuation : public Solver {
			private:
				mutable R random_gen;
				std::size_t homotopy_n;
				std::size_t local_solver_n;
			public:
				HomotopyContinuation(
						const R& random_gen,
						std::size_t homotopy_n,
						std::size_t local_solver_n) :
					random_gen(random_gen),
					homotopy_n(homotopy_n), local_solver_n(local_solver_n) {
					}

				R& randomGenerator() {
					return random_gen;
				}

				X solve(const ChemicalSystem& system) const override;
		};

		template<typename R>
		X HomotopyContinuation<R>::solve(const ChemicalSystem& system) const {
			ReducedChemicalSystem<CX> reduced_system(system);
			F<CX, CM> f(reduced_system, system);
			G g(reduced_system, random_gen);

			ChemicalLinearSolver<CM, CX> linear_solver(
					reduced_system.nComponents(),
					reduced_system.nSpecies()
					);
			Newton<CX, CM, I<CX>, ChemicalLinearSolver<CM, CX>> local_solver(
					linear_solver
					);
			Homotopy<CX, CM, Newton<CX, CM, I<CX>, ChemicalLinearSolver<CM, CX>>>
				homotopy_solver(local_solver);

			CHEM_LOG(TRACE) << "[HOMOTOPY] Init values: " << g.initValues();
			std::vector<SolverResult<CX>> solutions = homotopy_solver.solve(
					g.initValues(),
					[&f] (const CX& x) {return f.f(x);},
					[&f] (const CX& x) {return f.df(x);},
					[&g] (const CX& x) {return g.g(x);},
					[&g] (const CX& x) {return g.dg(x);},
					homotopy_n, local_solver_n
					);
			CX closest_solution(reduced_system.xSize());
			double current_minimum = std::numeric_limits<double>::infinity();
			auto sol_it = solutions.begin();
			while(sol_it != solutions.end()) {
				CHEM_LOG(TRACE) << "[HOMOTOPY] Possible solution: " << *sol_it;
				if(sol_it->isFinite()) {
					CX quantity_to_minimize;;
					for(const auto& c : sol_it->x) {
						if(c.real() < 0)
							quantity_to_minimize.push_back(c);
						else
							quantity_to_minimize.push_back({0.0, c.imag()});
					}
					auto _norm = norm(quantity_to_minimize);
					if(_norm < current_minimum) {
						closest_solution = sol_it->x;
						current_minimum = _norm;
					}
				}
				++sol_it;
			}
			CX complete_activities
				= reduced_system.completeActivities(closest_solution);
			X solution(complete_activities.size());
			for(std::size_t i = 0; i < complete_activities.size(); i++) {
				solution[i] = abs(complete_activities[i].real());
			}
			return solution;
		}
		/**
		 * Default equilibrium solver.
		 */
		extern AbsoluteNewton default_solver;
	}
}
#endif /*CHEMMISOL_SOLVER_H*/

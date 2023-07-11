#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>
#include <memory>
#include "linear.h"
#include "units.h"
#include "logging.h"

namespace chemmisol {
	class ChemicalSystem;

	enum Phase {
		SOLVENT, AQUEOUS, MINERAL
	};

	std::ostream& operator<<(std::ostream& o, const Phase& phase);

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
				const ChemicalSystem& system;
				std::vector<double> init_concentrations;

			public:
				/**
				 * Returns the resulting concentrations in the current chemical
				 * system that result from the provided extents.
				 */
				X concentrations(const X& extents) const;

				/**
				 * Initializes a function F to find the equilibrium of the
				 * provided chemical system.
				 */
				F(const ChemicalSystem& system);

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

	/**
	 * An abstract chemical component that might be involved in a reaction.
	 */
	class Component {
		public:
			/**
			 * Solution volume. Currently fixed to 1 liter.
			 *
			 * The volume is only required for internal calculation, but can be
			 * chosen arbitrarily since the quantity is always represented as
			 * concentrations for the user.
			 */
			static const double V;
		private:
			std::string name;
			std::size_t index;

		protected:
			/**
			 * Defines a simple component.
			 *
			 * @param name Name of the component
			 * @param index Index used to retrieve the component in data
			 * structures used internally by the ChemicalSystem. The index can
			 * also be used to uniquely identify each component of a system.
			 */
			Component(const std::string& name, std::size_t index)
				: name(name), index(index) {
				}

		public:
			/**
			 * Name of the component.
			 */
			const std::string& getName() const {
				return name;
			}

			/**
			 * Index used to retrieve the component in data
			 * structures used internally by the ChemicalSystem. The index can
			 * also be used to uniquely identify each component of a system.
			 */
			std::size_t getIndex() const {
				return index;
			}

			virtual Phase getPhase() const = 0;

			virtual bool isFixed() const {
				return false;
			}

			/**
			 * Increments the concentration of the current component,
			 * considering that its absolute quantity is incremented by
			 * `extent`.
			 *
			 * For components in solution, this is typically `C = C + extent/V`.
			 *
			 * @param extent Quantity of the current component added to the
			 * system
			 */
			virtual void incrementConcentration(double extent) = 0;

			/**
			 * Computes the concentration that would result if a quantity of
			 * `extent` was added to the current component at
			 * `current_concentration`.
			 *
			 * Notice that `current_concentration` might be different from
			 * concentration(), and that this call does **not** modify the
			 * concentration of the current component.
			 *
			 * @param current_concentration Concentration of the component
			 * @param extent Quantity of component
			 * @return Computed concentration of the component
			 */
			virtual double concentration(double current_concentration, double extent) const = 0;

			/**
			 * Returns the current concentration of the component.
			 *
			 * The theoretical definition of the concentration might vary
			 * depending on the nature of the component (aqueous, gaz,
			 * mineral...).
			 *
			 * @return Current concentration of the component
			 */
			virtual double concentration() const = 0;

			/**
			 * Returns the current activity of the component.
			 *
			 * The theoretical definition of the activity might vary depending
			 * on the nature of the component (aqueous, gaz, mineral...).
			 *
			 * @return Current activity of the component
			 */
			double activity() const {
				return activity(concentration());
			}

			/**
			 * Computes the activity of the current component corresponding to
			 * the specified `concentration`.
			 *
			 * Notice that `concentration` might be different from
			 * concentration(), and that this call does **not** modify the
			 * activity of the current component.
			 *
			 * @param concentration Concentration of the component
			 * @return Computed activity of the component
			 */
			virtual double activity(double concentration) const = 0;

			double quantity() const {
				return quantity(concentration());
			}

			virtual double quantity(double concentration) const = 0;

			/**
			 * Computes the activity that would result if a quantity of
			 * `extent` was added to the current component, considering
			 * `current_concentrations` as the concentrations of all components
			 * in the system.
			 *
			 * Notice that the concentration of the current component specified
			 * in `current_concentrations` might be different from
			 * concentration(), and that this call does **not** modify the
			 * activity of the current component.
			 *
			 * @param current_concentrations Concentrations of all components.
			 * The concentration of the current component can be retrieved with
			 * current_concentrations[getIndex()].
			 * @param extent Quantity of component
			 * @return Computed activity of the component
			 */
			virtual double activity(
					const std::vector<double>& current_concentrations,
					double extent) const = 0;

			/**
			 * Computes dc_k/dx_j where c_k denotes the concentration of the
			 * current component, and x_j denotes the extent of reaction j.
			 *
			 * See solver::F::df() for more detailed explanation about how this
			 * value is used and computed.
			 *
			 * @param current_concentrations Current concentrations of all
			 * components. The concentration of the current component can be
			 * retrieved with current_concentrations[getIndex()].
			 * @param d_coef Stoichiometric coefficient of the current component
			 * in the reaction k. Might be positive, negative or 0.
			 */
			virtual double Dactivity(
					const std::vector<double>& current_concentrations,
					double d_coef
					) const = 0;
			
			virtual ~Component() {
			}
	};

	class FixedComponent : public Component {
		private:
			double C;
			double A;
			double Q;

		protected:
			FixedComponent(
					const std::string& name, std::size_t index,
					double C, double A, double Q
					)
				: Component(name, index), C(C), A(A), Q(Q) {
				}

		public:
			bool isFixed() const override {
				return true;
			}

			void incrementConcentration(double) override {
			}

			double concentration(double, double) const override{
				return C;
			}

			double concentration() const override {
				return C;
			}

			double quantity(double) const override {
				return Q;
			}

			double activity(double) const override {
				return A;
			}

			double activity(
					const std::vector<double>&,
					double) const override {
				return A;
			}

			double Dactivity(
					const std::vector<double>&,
					double
					) const override {
				return 0;
			}
	};

	/**
	 * Solvent component.
	 *
	 * A solvent is a component always in excess, so that its activity is always
	 * 1.0 and its concentration does not matter. This is typically the case for
	 * water (H2O component).
	 */
	class Solvent : public FixedComponent {
		public:
			/**
			 * Defines a solvent.
			 */
			Solvent(const std::string& name, std::size_t id)
				: FixedComponent(name, id, 1*mol/l, 1.0, 1*mol/l*V) {
				}

			Phase getPhase() const override {
				return SOLVENT;
			}
	};

	/**
	 * Aqueous Component implementation.
	 *
	 * The concentration of an aqueous component is defined as `C=n/V` where n
	 * is the component quantity and V is the volume of the solution. A
	 * concentration is usually expressed in mol/L (molar).
	 *
	 * The activity of an aqueous component is defined as `C/C0` where C0 is the
	 * standard state concentration defined as 1mol/L.
	 */
	class AqueousComponent : public Component {
		private:
			double C;
		
		public:
			/**
			 * Defines a new AqueousComponent and initializes its concentration
			 * to C. It is the responsibility of the user to ensure unit
			 * consistencies. Predefined units can be used for this purpose, for
			 * example:
			 *
			 *     using namespace mineral;
			 *     AqueousComponent component("Na+", i, 0.1*mol/l);
			 */
			AqueousComponent(
					const std::string& name, std::size_t id,
					double C)
				: Component(name, id), C(C) {
				}

			Phase getPhase() const override {
				return AQUEOUS;
			}

			void incrementConcentration(double extent) override {
				C+=extent/V;
			}

			double concentration(double current_concentration, double extent) const override{
				return current_concentration + extent/V;
			}

			double concentration() const override {
				return C;
			}

			double quantity(double concentration) const override {
				return concentration*V;
			}

			double activity(double current_concentration) const override {
				return current_concentration/(1*mol/l);
			}

			double activity(
					const std::vector<double>& current_concentrations,
					double extent) const override {
				return activity(current_concentrations[getIndex()] + (extent/V));
			}

			double Dactivity(
					const std::vector<double>&,
					double d_coef
					) const override {
				// Note: this is 0 if d_coef=0, i.e. if the reaction k does not
				// use this component.
				return (d_coef/V)/(1*mol/l);
			}
	};

	class FixedAqueousComponent : public FixedComponent {
		public:
			FixedAqueousComponent(
					const std::string& name, std::size_t id,
					double C)
				: FixedComponent(name, id, C, C/(1*mol/l), C*V) {
				}

			Phase getPhase() const override {
				return AQUEOUS;
			}
	};

	/**
	 * Mineral Component implementation.
	 *
	 * A mineral component is a component adsorbed on a mineral (i.e. solid)
	 * phase to form _surface complexes_.
	 *
	 * In order to define this concept, we need to introduce the notion of _free
	 * sites_. A free site is a metal hydroxide component integrated in a
	 * mineral phase, at the interface with the surrounding solution. Free sites
	 * are typically noted `=SOH`, where S is a metal ion such as Fe or Al.
	 *
	 * Mineral components can then be formed fixing aqueous components on free
	 * sites. Here is some example reactions (charges are not represented except
	 * for H+ for the sake of simplicity), that form the `=SOH2` and `=SOPO3`
	 * MineralComponents:
	 *
	 *     =SOH + H+ <-> =SOH2
	 *     =SOH + PO4 + H+ <-> =SOPO3 + H2O
	 *
	 * An AqueousComponent can also occupy several free sites to form bidentate
	 * surface complexes:
	 *
	 *     2=SOH + PO4 + 2H+ <-> (=SOH)2 PO3 + 2H2O
	 *
	 * Such feature is easily represented by stoichiometric coefficients, and
	 * does not require any particular consideration as long as mineral
	 * components activities are properly defined.
	 *
	 * Indeed, the concentration of MineralComponents is rather considered as
	 * molar fractions (without unit) rather than in mole concentrations
	 * (mol/l). To define this fraction, we need to define N, the total quantity
	 * of possible free sites. Notice that this is different from the currently
	 * free sites: it is a constant representing **all** the surface sites,
	 * occupied and unoccupied. The molar fraction of the MineralComponent used
	 * as concentration() can then be defined as `n/N` where n denotes the
	 * current quantity of adsorbed MineralComponent. The activity of the
	 * MineralComponent can then be defined as strictly equal to this
	 * concentration.
	 *
	 * See [Hiemstra 1996](https://doi.org/10.1006/jcis.1996.0242) for more
	 * detailed and theoretical considerations about surface complexes.
	 */
	class MineralComponent : public Component {
		private:
			double fraction; // S/N
			double N;

		public:
			static double sites_count(
					double solid_concentration,
					double specific_surface_area,
					double site_concentration);
			/**
			 * Defines a MineralComponent.
			 *
			 * N is calculated as `V * solid_concentration *
			 * specific_surface_area * site_concentration`. Where V is the
			 * arbitrary and constant volume of the solution defined
			 * AqueousSpecies::V.
			 *
			 * The molar fraction of the component is initialized from the
			 * specified `fraction` parameter.
			 *
			 * @param name Name of the component
			 * @param index Unique index of the component
			 * @param solid_concentration Quantity of mineral in suspension in
			 * the solution, usually expressed in g/l.
			 * @param specific_surface_area Surface of the solid in contact with
			 * the solution per unit of mass, usually expressed in m2/g.
			 * @param site_concentration Quantity of sites per unit of surface
			 * in contact with the solution, usually expressed as sites/nm2.
			 *
			 * It is the responsibility of the user to ensure unit
			 * consistencies. Predefined units can be used for this purpose, for
			 * example:
			 *
			 *     using namespace mineral;
			 *     MineralComponent component(
			 *         "=SOH2", i,
			 *         0.1, // molar fraction, without unit
			 *         2.5 * g/l, // solid concentration
			 *         24.2 * m2/g, // specific surface area
			 *         0.8 * entities/nm2 // site concentration
			 *     );
			 *
			 * Notice that the previous example is actually **consistent**, even
			 * if the specific surface area is specified as **m2** /g and the
			 * site concentration is specified in entities/ **nm2**, since the
			 * purpose of the `*unit` computation is to convert each numerical
			 * value to the standard and consistent unit system used internally.
			 */
			// TODO: it is not efficient to specify solid parameters for ALL
			// MineralComponents
			MineralComponent(
					const std::string& name, std::size_t index,
					double solid_concentration,
					double specific_surface_area,
					double site_concentration,
					double fraction)
				: Component(name, index), fraction(fraction),
				N(sites_count(solid_concentration, specific_surface_area, site_concentration)) {
				}
			
			Phase getPhase() const override {
				return MINERAL;
			}

			void incrementConcentration(double extent) override {
				fraction+=extent/N;
			}

			double concentration(double current_concentration, double extent) const override{
				return current_concentration + extent/N;
			}
			double concentration() const override {
				return fraction;
			}

			double quantity(double fraction) const override {
				return fraction * N;
			}

			double activity(double current_concentration) const override {
				return current_concentration;
			}

			double activity(
					const std::vector<double>& current_concentrations,
					double extent) const override {
				return activity(current_concentrations[getIndex()]+(extent/N));
			}

			double Dactivity(
					const std::vector<double>&,
					double d_coef
					) const override {
				// Note: this is 0 if d_coef=0, i.e. if the reaction k does not
				// use this component.
				return d_coef/N;
			}

	};

	class FixedMineralComponent : public FixedComponent {
		public:
			FixedMineralComponent(
					const std::string& name, std::size_t id,
					double solid_concentration,
					double specific_surface_area,
					double site_concentration,
					double F
					)
				: FixedComponent(name, id, F, F,
						MineralComponent::sites_count(
							solid_concentration, specific_surface_area,
							site_concentration
							)*F) {
				}

			Phase getPhase() const override {
				return MINERAL;
			}
	};

	struct ReactionComponent {
		std::string name;
		Phase phase;
		double coefficient;

		ReactionComponent() = default;
		ReactionComponent(
				const std::string& name,
				Phase phase,
				double coefficient)
			: name(name), phase(phase), coefficient(coefficient) {
			}
		ReactionComponent(
				const std::string& name,
				double coefficient)
			: name(name), phase(name == "H2O" ? SOLVENT : AQUEOUS), coefficient(coefficient) {
			}
	};

	bool operator==(const ReactionComponent& c1, const ReactionComponent& c2);

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
			double log_K;
			std::vector<ReactionComponent> reagents;

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
					const std::vector<ReactionComponent>& reagents
					)
				: name(name), index(index), log_K(log_K), reagents(reagents) {
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

			/**
			 * Base 10 logarithm of the equilibrium constant.
			 */
			double getLogK() const {
				return log_K;
			}

			/**
			 * Reagents of the reaction.
			 */
			const std::vector<ReactionComponent> getReagents() const {
				return reagents;
			}
	};

	class GuessF {
		private:
			const ChemicalSystem& system;
			const Reaction& reaction;
			const std::vector<double>& current_concentrations;

		public:
			GuessF(
					const ChemicalSystem& system,
					const Reaction& reaction,
					const std::vector<double>& current_concentrations
					) :
				system(system), reaction(reaction),
				current_concentrations(current_concentrations) {
				}
			double operator()(const double& extent) const;
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
			std::size_t component_index = 0;
			std::size_t reaction_index = 0;

			std::vector<std::unique_ptr<Component>> components;
			std::unordered_map<std::string, const Component*> components_by_name;
			std::vector<std::unique_ptr<Reaction>> reactions;
			std::unordered_map<std::string, const Reaction*> reactions_by_name;
			std::vector<std::vector<double>> reaction_matrix;

			void addComponent(Component* component);
			void fixComponent(FixedComponent* component);
			void addReaction(Reaction* reaction);

			std::size_t max_iteration = 200;

			// Adsorption model parameters
			double solid_concentration;
			double specific_surface_area;
			double site_concentration;

			std::unordered_map<std::string, double> initial_guess_extents;


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
					std::vector<ReactionComponent> reagents
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
			 * Gets the component named `name`.
			 *
			 * The behavior of the method is unspecified if the name does not
			 * correspond to any component added with addComponent() (called by
			 * the user or initialized by default by initReactionMatrix()).
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
			 * Returns references to all the components in the system.
			 */
			const std::vector<std::unique_ptr<Component>>& getComponents() const {
				return components;
			}

			/**
			 * Returns references to all reactions in the system.
			 */
			const std::vector<std::unique_ptr<Reaction>>& getReactions() const {
				return reactions;
			}

			void setInitialGuessExtent(const std::string& reaction, double guess) {
				initial_guess_extents[reaction] = guess;
				CHEM_LOG(INFO) << "User specified extent guess for " << reaction
					<< ": " << guess;
			}

			double getInitialGuessExtent(const std::string& reaction) const {
				auto guess = initial_guess_extents.find(reaction);
				if(guess != initial_guess_extents.end())
					return guess->second;
				return 0;
			}

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

			/**
			 * Initializes the reaction matrix from all components and reactions
			 * added to the system until now.
			 *
			 * Automatically called by solveEquilibrium(), but might be called
			 * by the user.
			 */
			void initReactionMatrix();

			const std::unordered_map<std::string, double>& guessInitialExtents();

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

/*
 *    class Soil {
 *        public:
 *        static const double bulk_density;
 *        static const double V;
 *        static const double mineral_weight;
 *
 *        private:
 *
 *        double _K1;
 *        double _K2;
 *        double _K3;
 *        double _initPH;
 *
 *        double H;
 *        double P;
 *        double C;
 *
 *        double N;
 *        double S;
 *        double SH;
 *        double SP;
 *        double SC;
 *
 *        double concentration(const double& n) const;
 *
 *        private:
 *        ChemicalSystem(
 *                double K1, double K2, double K3, double initPH,
 *                double H, double P, double C,
 *                double N,
 *                double S, double SH, double SP, double SC
 *                );
 *
 *        public:
 *        static ChemicalSystem equilibrium(
 *                double K1, double K2, double K3,
 *                double pH, double solution_P, double solution_C,
 *                double mineral_N
 *                );
 *        static ChemicalSystem defaultEquilibrium();
 *        static ChemicalSystem Devau2011Control();
 *
 *        //static ChemicalSystem soilParameters(
 *        //double K1, double K2, double K3,
 *        //double pH, double soil_P, double soil_C,
 *        //double mineral_N
 *        //);
 *        //static ChemicalSystem defaultSoil();
 *
 *
 *        void incrementP(double P);
 *        void setPH(double pH);
 *
 *        //X reactionQuotient() const;
 *        void distanceToEquilibrium() const;
 *
 *        double cH() const;
 *        double cP() const;
 *        double cC() const;
 *        double cS() const;
 *        double cSH() const;
 *        double cSP() const;
 *        double cSC() const;
 *
 *        double nH() const;
 *        double nP() const;
 *        double nC() const;
 *        double nS() const;
 *        double nSH() const;
 *        double nSP() const;
 *        double nSC() const;
 *
 *        double K1() const;
 *        double K2() const;
 *        double K3() const;
 *        double initPH() const;
 *
 *        double pH() const;
 *    };
 */
}

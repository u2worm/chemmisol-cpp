#ifndef CHEMMISOL_SPECIES_H
#define CHEMMISOL_SPECIES_H

#include <string>
#include <vector>
#include "../units.h"

namespace chemmisol {
	enum Phase {
		SOLVENT, AQUEOUS, MINERAL
	};

	std::ostream& operator<<(std::ostream& o, const Phase& phase);

	class ChemicalSpecies {
		private:
			std::string name;
			std::size_t index;

		protected:
			/**
			 * Defines a simple chemical species.
			 *
			 * @param name Name of the species
			 * @param index Index used to retrieve the species in data
			 * structures used internally by the ChemicalSystem. The index can
			 * also be used to uniquely identify each species of a system.
			 */
			ChemicalSpecies(const std::string& name, std::size_t index)
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

			virtual void setActivity(double activity) = 0;

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
			
			virtual ~ChemicalSpecies() {
			}
	};

	
	class FixedChemicalSpecies : public ChemicalSpecies {
		private:
			double C;
			double A;
			double Q;

		protected:
			FixedChemicalSpecies(
					const std::string& name, std::size_t index,
					double C, double A, double Q
					)
				: ChemicalSpecies(name, index), C(C), A(A), Q(Q) {
				}

		public:
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

			void setActivity(double) override {
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
	 * Aqueous Component implementation.
	 *
	 * The concentration of an aqueous component is defined as `C=n/V` where n
	 * is the component quantity and V is the volume of the solution. A
	 * concentration is usually expressed in mol/L (molar).
	 *
	 * The activity of an aqueous component is defined as `C/C0` where C0 is the
	 * standard state concentration defined as 1mol/L.
	 */
	class AqueousSpecies : public ChemicalSpecies {
		private:
			double C;
		
		public:
			/**
			 * Solution volume. Currently fixed to 1 liter.
			 *
			 * The volume is only required for internal calculation, but can be
			 * chosen arbitrarily since the quantity is always represented as
			 * concentrations for the user.
			 */
			static const double V;

			/**
			 * Defines a new AqueousComponent and initializes its concentration
			 * to C. It is the responsibility of the user to ensure unit
			 * consistencies. Predefined units can be used for this purpose, for
			 * example:
			 *
			 *     using namespace mineral;
			 *     AqueousComponent component("Na+", i, 0.1*mol/l);
			 */
			AqueousSpecies(
					const std::string& name, std::size_t id,
					double C)
				: ChemicalSpecies(name, id), C(C) {
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

			void setActivity(double activity) override {
				C = activity * 1*mol/l;
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

	class FixedAqueousSpecies : public FixedChemicalSpecies {
		public:
			FixedAqueousSpecies(
					const std::string& name, std::size_t id,
					double C)
				: FixedChemicalSpecies(name, id, C, C/(1*mol/l), C*AqueousSpecies::V) {
				}

			Phase getPhase() const override {
				return AQUEOUS;
			}
	};

	/**
	 * Solvent component.
	 *
	 * A solvent is a component always in excess, so that its activity is always
	 * 1.0 and its concentration does not matter. This is typically the case for
	 * water (H2O component).
	 */
	class Solvent : public FixedAqueousSpecies {
		public:
			/**
			 * Defines a solvent.
			 */
			Solvent(const std::string& name, std::size_t id)
				: FixedAqueousSpecies(name, id, 1*mol/l) {
				}

			Phase getPhase() const override {
				return SOLVENT;
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
	class MineralSpecies : public ChemicalSpecies {
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
			MineralSpecies(
					const std::string& name, std::size_t index,
					double solid_concentration,
					double specific_surface_area,
					double site_concentration,
					double fraction)
				: ChemicalSpecies(name, index), fraction(fraction),
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

			void setActivity(double activity) override {
				fraction = activity;
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

	class FixedMineralSpecies : public FixedChemicalSpecies {
		public:
			FixedMineralSpecies(
					const std::string& name, std::size_t id,
					double solid_concentration,
					double specific_surface_area,
					double site_concentration,
					double F
					)
				: FixedChemicalSpecies(name, id, F, F,
						MineralSpecies::sites_count(
							solid_concentration, specific_surface_area,
							site_concentration
							)*F) {
				}

			Phase getPhase() const override {
				return MINERAL;
			}
	};

	class ChemicalComponent {
		private:
			ChemicalSpecies* species;
			std::size_t index;
			double total_quantity;

		public:
			ChemicalComponent(
					ChemicalSpecies* species,
					std::size_t index, double total_quantity)
				: species(species), index(index), total_quantity(total_quantity) {
				}

			ChemicalSpecies* getSpecies() {
				return species;
			}

			const ChemicalSpecies* getSpecies() const {
				return species;
			}

			std::size_t getIndex() const {
				return index;
			}

			virtual bool isFixed() const {
				return false;
			}

			double getTotalQuantity() const {
				return total_quantity;
			}

			void setTotalQuantity(double quantity) {
				total_quantity = quantity;
			}

			virtual ~ChemicalComponent() {
			}
	};

	class FixedChemicalComponent : public ChemicalComponent {
		public:
			FixedChemicalComponent(
					FixedChemicalSpecies* species,
					std::size_t index, double total_quantity)
				: ChemicalComponent(species, index, total_quantity) {
				}

			virtual bool isFixed() const override {
				return true;
			}
	};

}
#endif /*CHEMMISOL_SPECIES_H*/

#ifndef CHEMMISOL_SPECIES_H
#define CHEMMISOL_SPECIES_H

#include <vector>
#include "../units.h"

namespace chemmisol {
	/**
	 * Possible phases of a chemical species.
	 */
	enum Phase {
		SOLVENT, AQUEOUS, MINERAL
	};

	/**
	 * Phases stream output operator.
	 */
	std::ostream& operator<<(std::ostream& o, const Phase& phase);

	/**
	 * Describes a chemical species contained in a chemical system. A quantity()
	 * is assigned to each chemical species in a system, typically in mole, and
	 * a concentration() and an activity() corresponding to this quantity is
	 * also defined, depending on the phase and type of the chemical species.
	 */
	class ChemicalSpecies {
		private:
			std::string name;
			std::size_t index;

		protected:
			/**
			 * Defines a simple chemical species.
			 *
			 * @param name Name of the species.
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
			 * Index used to retrieve the species in data structures used
			 * internally by the ChemicalSystem. The index can also be used to
			 * uniquely identify each chemical species of a system.
			 */
			std::size_t getIndex() const {
				return index;
			}

			/**
			 * Phase of the chemical species.
			 */
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
			 * The implementation is equivalent to a call to
			 * `activity(concentration())`, so that the implementation of the
			 * activity(double concentration) method defines the value returned
			 * by the activity() method.
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

			/**
			 * Returns the current quantity of the chemical species.
			 *
			 * The quantity must be understood as a count of items, for example
			 * the count of molecules corresponding to the chemical species in
			 * the system, typically expressed in mole.
			 *
			 * The implementation is equivalent to a call to
			 * `quantity(concentration())`, so that the implementation of the
			 * quantity(double concentration) method defines the value returned
			 * by the quantity() method.
			 *
			 * @return Quantity of the chemical species in the chemical system
			 */
			double quantity() const {
				return quantity(concentration());
			}

			/**
			 * Computes the quantity of the current component corresponding to
			 * the specified `concentration`.
			 *
			 * Notice that `concentration` might be different from
			 * concentration(), and that this call does **not** modify the
			 * quantity of the current component.
			 *
			 * @param concentration Concentration of the component
			 * @return Computed quantity of the component
			 */
			virtual double quantity(double concentration) const = 0;

			/**
			 * Sets the activity of the chemical species. The quantity() and
			 * concentration() values are updated accordingly.
			 */
			virtual void setActivity(double activity) = 0;

			virtual ~ChemicalSpecies() {
			}
	};

	
	/**
	 * ChemicalSpecies implementation for a species with fixed quantity,
	 * activity and concentration in the chemical system.
	 *
	 * The quantity, activity and concentration of such species will never
	 * change for any reaction process. This does not mean that a reaction with
	 * a fixed reagent cannot process, but that the fixed activity is always
	 * used in the reaction quotient computation for any extent of the reaction.
	 *
	 * Notice that only components can be fixed within a ChemicalSystem:
	 * FixedChemicalSpecies are supposed to be associated to
	 * FixedChemicalComponent.
	 *
	 * @see FixedChemicalComponent
	 * @see ChemicalSystem::fixComponent()
	 */
	class FixedChemicalSpecies : public ChemicalSpecies {
		private:
			double C;
			double A;
			double Q;

		protected:
			/**
			 * FixedChemicalSpecies constructor.
			 *
			 * Even if the quantity, activity and concentration must be manually
			 * specified in the FixedChemicalSpecies constructor, they are
			 * generally bound according to physical laws. Such rules depends on
			 * the nature of the species and are implemented in Solvent,
			 * FixedAqueousSpecies or FixedMineralSpecies subclasses.
			 *
			 * @param name Name of the fixed species
			 * @param index Index used to retrieve the species in data
			 * structures used internally by the ChemicalSystem. The index can
			 * also be used to uniquely identify each species of a system.
			 * @param C Fixed concentration
			 * @param A Fixed activity
			 * @param Q Fixed quantity
			 */
			FixedChemicalSpecies(
					const std::string& name, std::size_t index,
					double C, double A, double Q
					)
				: ChemicalSpecies(name, index), C(C), A(A), Q(Q) {
				}

		public:
			/**
			 * This method has no effect for a fixed component: the fixed
			 * concentration is not changed.
			 */
			void incrementConcentration(double) override {
			}

			/**
			 * Returns the fixed concentration of this chemical species.
			 */
			double concentration() const override {
				return C;
			}

			/**
			 * Returns the fixed quantity of this chemical species, ignoring the
			 * provided concentration.
			 */
			double quantity(double /*concentration*/) const override {
				return Q;
			}

			/**
			 * Returns the fixed activity of this chemical species, ignoring the
			 * provided concentration.
			 */
			double activity(double /*concentration*/) const override {
				return A;
			}

			/**
			 * This method has no effect for a fixed component: the fixed
			 * activity is not changed.
			 */
			void setActivity(double) override {
			}
	};

	/**
	 * Aqueous species implementation.
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
			 * Notice that the concentration() and activity() of aqueous species
			 * do not depend on the volume.
			 *
			 * If quantities for a different volume are required, it is
			 * currently recommend to use concentration() as output and to
			 * manually convert the result to a quantity with the formula
			 * `n=C*V`.
			 *
			 * Notice that chemmisol currently does not handle reactions with a
			 * variable volume.
			 */
			static const double V;

			/**
			 * Defines a new AqueousSpecies and initializes its concentration to
			 * C. It is the responsibility of the user to ensure unit
			 * consistency. Predefined units can be used for this purpose, for
			 * example:
			 *
			 *     using namespace mineral;
			 *     AqueousSpecies species("Na+", i, 0.1*mol/l);
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

			double concentration() const override {
				return C;
			}

			double quantity(double concentration) const override {
				return concentration*V;
			}

			double activity(double concentration) const override {
				return concentration/(1*mol/l);
			}

			void setActivity(double activity) override {
				C = activity * 1*mol/l;
			}
	};

	/**
	 * Fixed aqueous species implementation.
	 */
	class FixedAqueousSpecies : public FixedChemicalSpecies {
		public:
			/**
			 * Defines a fixed aqueous species. Quantity and activity are fixed
			 * according to the provided concentration.
			 *
			 * @param name Name of the species.
			 * @param index Index used to retrieve the species in data
			 * structures used internally by the ChemicalSystem. The index can
			 * also be used to uniquely identify each species of a system.
			 * @param C Fixed concentration.
			 */
			FixedAqueousSpecies(
					const std::string& name, std::size_t index,
					double C)
				: FixedChemicalSpecies(name, index, C, C/(1*mol/l), C*AqueousSpecies::V) {
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
	 * Mineral species implementation.
	 *
	 * A mineral component is a component adsorbed on a mineral (i.e. solid)
	 * phase to form _surface complexes_.
	 *
	 * In order to define this concept, we need to introduce the notion of _free
	 * sites_. A free site is a metal hydroxide component integrated in a
	 * mineral phase, at the interface with the surrounding solution. Free sites
	 * are typically noted `=SOH`, where S is a metal ion such as Fe or Al.
	 *
	 * Mineral species can then be formed fixing aqueous species on free sites.
	 * Here is some example reactions (charges are not represented except for H+
	 * for the sake of simplicity), that form the `=SOH2` and `=SOPO3`
	 * MineralComponents:
	 *
	 *     =SOH + H+ <-> =SOH2
	 *     =SOH + PO4 + H+ <-> =SOPO3 + H2O
	 *
	 * An AqueousSpecies can also occupy several free sites to form bidentate
	 * surface complexes:
	 *
	 *     2=SOH + PO4 + 2H+ <-> (=SOH)2 PO3 + 2H2O
	 *
	 * Such feature is easily represented by stoichiometric coefficients, and
	 * does not require any particular consideration as long as mineral
	 * components activities are properly defined.
	 *
	 * Indeed, the concentration of a MineralSpecies is rather considered as a
	 * molar fraction (without unit) rather than a mole concentration (mol/l).
	 * To define this fraction, we need to define N, the total quantity of
	 * possible free sites. Notice that this is different from the currently
	 * free sites: it is a constant representing **all** the surface sites,
	 * occupied and unoccupied. The molar fraction of the MineralSpecies used as
	 * concentration() is then defined as `n/N` where n denotes the current
	 * quantity of adsorbed MineralSpecies. The activity of the MineralSpecies
	 * is defined equal to this concentration.
	 *
	 * See [Hiemstra 1996](https://doi.org/10.1006/jcis.1996.0242) for more
	 * detailed and theoretical considerations about surface complexes.
	 */
	class MineralSpecies : public ChemicalSpecies {
		private:
			double fraction; // S/N
			double N;

		public:
			/**
			 * Computes the total sites count that corresponds to the provided
			 * parameters.
			 *
			 * @param solid_concentration Quantity of mineral in suspension in
			 * the solution, usually expressed in g/l.
			 * @param specific_surface_area Surface of the solid in contact with
			 * the solution per unit of mass, usually expressed in m2/g.
			 * @param site_concentration Quantity of sites per unit of surface
			 * in contact with the solution, usually expressed as sites/nm2.
			 */
			static double sites_count(
					double solid_concentration,
					double specific_surface_area,
					double site_concentration);

			/**
			 * Defines a MineralSpecies.
			 *
			 * N is calculated as `V * solid_concentration *
			 * specific_surface_area * site_concentration`. Where V is the
			 * arbitrary and constant volume of the solution defined in
			 * AqueousSpecies::V.
			 *
			 * The molar fraction of the component is initialized from the
			 * specified `fraction` parameter.
			 *
			 * @param name Name of the component.
			 * @param index Unique index of the component.
			 * @param solid_concentration Quantity of mineral in suspension in
			 * the solution, usually expressed in g/l.
			 * @param specific_surface_area Surface of the solid in contact with
			 * the solution per unit of mass, usually expressed in m2/g.
			 * @param site_concentration Quantity of sites per unit of surface
			 * in contact with the solution, usually expressed as sites/nm2.
			 * @param fraction Initial molar fraction of this mineral species.
			 *
			 * It is the responsibility of the user to ensure unit
			 * consistency. Predefined units can be used for this purpose, for
			 * example:
			 *
			 *     using namespace mineral;
			 *     MineralSpecies species(
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
			// MineralSpecies
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
	};

	/**
	 * Fixed mineral species implementation.
	 */
	class FixedMineralSpecies : public FixedChemicalSpecies {
		public:
			/**
			 * Defines a fixed mineral species. Quantity and activity are fixed
			 * according to the provided concentration expressed as a molar
			 * fraction.
			 *
			 * @param name Name of the species.
			 * @param index Index used to retrieve the species in data
			 * structures used internally by the ChemicalSystem. The index can
			 * also be used to uniquely identify each species of a system.
			 * @param solid_concentration Quantity of mineral in suspension in
			 * the solution, usually expressed in g/l.
			 * @param specific_surface_area Surface of the solid in contact with
			 * the solution per unit of mass, usually expressed in m2/g.
			 * @param site_concentration Quantity of sites per unit of surface
			 * in contact with the solution, usually expressed as sites/nm2.
			 * @param F Fixed concentration expressed as a molar fraction.
			 */
			FixedMineralSpecies(
					const std::string& name, std::size_t index,
					double solid_concentration,
					double specific_surface_area,
					double site_concentration,
					double F
					)
				: FixedChemicalSpecies(name, index, F, F,
						MineralSpecies::sites_count(
							solid_concentration, specific_surface_area,
							site_concentration
							)*F) {
				}

			Phase getPhase() const override {
				return MINERAL;
			}
	};

	/**
	 * Describes a chemical component that can be used to define reactions in a
	 * chemical system.
	 *
	 * Any component is automatically associated to a ChemicalSpecies with the
	 * same name.
	 *
	 * A reaction system is defined from two kinds of chemical species:
	 * components and compounds.
	 *
	 * A component is a canonical species that cannot be divided, i.e. that
	 * cannot be expressed as a species produced from other components. On the
	 * other hand, compound species can always be expressed from components. It
	 * is the responsibility of the user to properly specify components and
	 * reactions in a ChemicalSystem so that the reaction system is consistent.
	 * Exceptions such as MissingProducedSpeciesInReaction or
	 * InvalidSpeciesInReaction might be thrown when it is not the case.
	 *
	 * Chemical components are notably associated to a _total quantity_, a very
	 * important concept in the context of equilibrium solving.
	 *
	 * When a component is **fixed**, the quantity of its associated species is
	 * always equal to the total quantity of the component.
	 *
	 * Otherwise, the total quantity of a component is defined as the sum of the
	 * quantity of species associated to the component itself and the quantities
	 * of all species compound from this component, weighted by the
	 * stoichiometric coefficient associated to the component in the reaction
	 * producing exactly one unit of the compound.
	 *
	 * Solvents (such as H2O) are currently considered as fixed components in
	 * Chemmisol.
	 *
	 * @par Example
	 * \parblock
	 *
	 * Let's define a system with the following components: H+ and PO4.
	 * 
	 * The following reactions can be defined to introduce the OH- and H3PO4
	 * produced species as follows:
	 *
	 *     H2O       <-> OH- + H+
	 *     PO4 + 3H+ <-> H3PO4
	 *
	 * The first reaction can be rewritten as follows so that it produces
	 * exactly one unit of OH-:
	 *
	 *     H2O - H+ <-> OH-
	 *
	 * The total quantity N of the PO4 and H+ components are then defined as follows:
	 *
	 *     N(PO4) = n(PO4) + n(H3PO4)
	 *     N(H+) = n(H+) - n(OH-) + 3 * n(H3PO4)
	 *
	 * where n denotes the quantity of each chemical species. Notice that n(PO4)
	 * and n(H+) denote the quantities of PO4 and H+ **species**, that are
	 * likely not equal to total quantities of PO4 and H+ **components**.
	 *
	 * It is the purpose of the ChemicalSystem equilibrium solver to dispatch
	 * the total quantities of components between compound species so that the
	 * above relations actually holds. This constraint is known in physics as
	 * the [law of conservation of
	 * mass](https://en.wikipedia.org/wiki/Conservation_of_mass). See
	 * ChemicalSystem::massConservationLaw() for more details.
	 *
	 * \endparblock
	 */
	class ChemicalComponent {
		private:
			ChemicalSpecies* species;
			std::size_t index;
			double total_quantity;

		public:
			/**
			 * Defines a ChemicalComponent.
			 *
			 * @param species The chemical species instance associated to this
			 * component.
			 * @param index Index used to retrieve the component in data
			 * structures used internally by the ChemicalSystem. The index can
			 * also be used to uniquely identify each component of a system.
			 * @param total_quantity Initial total quantity of the component.
			 */
			ChemicalComponent(
					ChemicalSpecies* species,
					std::size_t index, double total_quantity)
				: species(species), index(index), total_quantity(total_quantity) {
				}

			/**
			 * Returns a reference to the species associated to this component.
			 */
			ChemicalSpecies* getSpecies() {
				return species;
			}

			/**
			 * \copydoc getSpecies()
			 */
			const ChemicalSpecies* getSpecies() const {
				return species;
			}

			/**
			 * Index used to retrieve the component in data structures used
			 * internally by the ChemicalSystem. The index can also be used to
			 * uniquely identify each component of a system.
			 */
			std::size_t getIndex() const {
				return index;
			}

			/**
			 * Returns true if this component is fixed. Returns false by
			 * default, but can be overridden by subclasses.
			 */
			virtual bool isFixed() const {
				return false;
			}

			/**
			 * Returns the total quantity of this component.
			 */
			double getTotalQuantity() const {
				return total_quantity;
			}

			/**
			 * Sets the total quantity of this component.
			 *
			 * When the component is part of a ChemicalSystem,
			 * ChemicalSystem::setTotalQuantity() should be used.
			 *
			 * @warning
			 * The quantity of species produced from this component, including
			 * the species associated to this component, are not updated by this
			 * method. The ChemicalSystem equilibrium solver must be run to
			 * update species quantities according the new equilibrium state and
			 * the mass conservation law.
			 */
			void setTotalQuantity(double quantity) {
				total_quantity = quantity;
			}

			/**
			 * Sets the total quantity of the component to a quantity that
			 * corresponds to the quantity of the associated species at the
			 * provided concentration.
			 *
			 * Equivalent to
			 * `setTotalQuantity(species->quantity(concentration))`.
			 *
			 * When the component is part of a ChemicalSystem,
			 * ChemicalSystem::setTotalConcentration() should be used.
			 *
			 * @see setTotalQuantity()
			 */
			void setTotalConcentration(double concentration) {
				setTotalQuantity(species->quantity(concentration));
			}

			virtual ~ChemicalComponent() {
			}
	};

	/**
	 * ChemicalComponent implementation for a fixed component.
	 *
	 * Notice that calling setTotalQuantity() on a FixedChemicalComponent is
	 * valid.
	 *
	 * - When setTotalQuantity() is called on a regular component, the total
	 *   quantity is dispatched among produced species when the equilibrium
	 *   solver is run.
	 * - When setTotalQuantity() is called on a fixed component, the quantity of
	 *   the associated
	 *   species is fixed to the provided value.
	 *
	 * However, fixing the quantity of the associated species is not handled by
	 * the FixedChemicalComponent itself, but by the
	 * ChemicalSystem::fixComponent() method.
	 *
	 * @see ChemicalSystem::fixComponent()
	 */
	class FixedChemicalComponent : public ChemicalComponent {
		public:
			/**
			 * Defines a FixedChemicalComponent.
			 *
			 * The total quantity of the fixed component is initialized from the
			 * quantity of the associated fixed species.
			 *
			 * @param species The fixed chemical species instance associated to
			 * this component.
			 * @param index Index used to retrieve the component in data
			 * structures used internally by the ChemicalSystem. The index can
			 * also be used to uniquely identify each component of a system.
			 */
			FixedChemicalComponent(
					FixedChemicalSpecies* species,
					std::size_t index)
				: ChemicalComponent(species, index, species->ChemicalSpecies::quantity()) {
				}

			virtual bool isFixed() const override {
				return true;
			}
	};

}
#endif /*CHEMMISOL_SPECIES_H*/

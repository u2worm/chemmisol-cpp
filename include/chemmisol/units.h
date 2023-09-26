#ifndef CHEMMISOL_UNIT_H
#define CHEMMISOL_UNIT_H

/**
 * @file chemmisol/units.h
 *
 * Defines units related features.
 *
 * @par Examples
 * \ref basic_chemical_system/main.cpp
 */

#include <unordered_map>
#include <string>

/**
 * \brief Defines available units.
 *
 * Chemmisol defines a unit system inspired from [GAMA
 * units](https://gama-platform.org/wiki/next/UnitsAndConstants).
 *
 * However, even if very useful when directly using the library, the unit system
 * might be problematic when attempting to embed the library in other systems or
 * languages.
 *
 * In order to solve those issues, it is necessary to enforce the specification
 * of the unit system used by chemmisol.
 *
 * Any quantity specified with [units](#UNITS) is converted internally to a
 * value that corresponds to the same value expressed in the _core unit_
 * of the same category.
 *
 * In consequence, it is safe to use raw values without units in chemmisol, as
 * long as it can be ensure that it's value corresponds to a value expressed in
 * a core unit.
 *
 * Internally, core units corresponds to #UNITS() with a value of 1.
 *
 * For example, the core unit for volumes is the liter `l`. This means that the
 * value `1.5*l` is represented internally by `1.5f`, while `1.2*cm3` is
 * represented as `1.2 * 10-3 = 1.2e-3f`. The interest of this system is that
 * unit conversion is automatically performed by the library. For example, it is
 * perfectly possible to write `double distance = 1*km + 15*m` without bothering
 * about unit conversion (in this example, `distance` is automatically converted
 * to the core unit `m` with a value of `1015.0`). Considering the definition of
 * core units, specifying directly `double distance = 1015.0` is equivalent in
 * chemmisol.
 *
 * A value defined with a chemmisol unit can still be expressed in an other
 * unit, dividing the value by the desired unit. For example:
 * ```cpp
 * double distance = 1*km + 15*m;
 * std::cout << "Distance in km: " << distance/km << std::endl;
 * ```
 *
 * The following table sums up all the core units used by chemmisol:
 *
 * | Use case | Unit |
 * |----------|------|
 * | Distance | m    |
 * | Surface  | m2   |
 * | Volume   | l    |
 * | Mass     | g    |
 * | Quantity | mol  |
 *
 * The core unit of derived units can be obtained with the same expression as
 * the derived unit expressed in core units. For example, the core unit for
 * molar concentrations is `mol/l`.
 *
 * @par Example
 *
 * ```cpp
 * ChemicalSystem chemical_system {
 * 	2.5 * g/l,
 * 	24.2 * m2/g,
 * 	0.8 * entities/nm2,
 * 	"=SOH"
 * };
 * ```
 * is equivalent to
 * ```cpp
 *ChemicalSystem chemical_system {
 * 	2.5,
 * 	24.2,
 * 	0.8 / 6.02214076e23 * 1e18, // Manual conversion from entities/nm2 to mol/m2
 * 	"=SOH"
 * };
 * ```
 *
 * @warning
 * Even if the unit system is inspired from [GAMA
 * units](https://gama-platform.org/wiki/next/UnitsAndConstants), it is not
 * guaranteed that the same core units are used. In consequence, a GAMA raw
 * value cannot be used directly in chemmisol, and must be converted to
 * chemmisol core units using the unit system of GAMA.

 */
#define UNITS(F)\
	F(m, 1, LENGTH)\
	F(cm, 1e-2*m, LENGTH)\
	F(nm, 1e-9*m, LENGTH)\
	F(nm2, nm*nm, SURFACE)\
	F(m2, m*m, SURFACE)\
	F(l, 1, VOLUME)\
	F(m3, 1000*l, VOLUME)\
	F(cm3, 10e-3*l, VOLUME)\
	F(g, 1, MASS)\
	F(gram, 1, MASS)\
	F(kg, 1e3*gram, MASS)\
	F(mg, 1e-3*gram, MASS)\
	F(mol, 1, QUANTITY)\
	F(entities, 1/NA, QUANTITY)\
	F(u, 1.66053906660e-27*kg, MASS)

/**
 * Defines an unit.
 */
#define UNIT_DEF(NAME, VALUE, CATEGORY) extern const double NAME;

/**
 * Defines the value of an unit.
 */
#define UNIT_VALUE(NAME, VALUE, CATEGORY) const double NAME = VALUE;

/**
 * Defines the caterory of an unit.
 */
#define UNIT_NAME_TO_CATEGORY(NAME, VALUE, CATEGORY) {#NAME, Category::CATEGORY},

namespace chemmisol {
	/**
	 * The [Avogadro constant](https://en.wikipedia.org/wiki/Avogadro_constant).
	 */
	extern const double NA;

	/**
	 * Describes the category of each unit.
	 */
	enum Category {
		MASS, LENGTH, VOLUME, SURFACE, QUANTITY
	};

	UNITS(UNIT_DEF)
	static const std::unordered_map<std::string, Category> units_category {
		UNITS(UNIT_NAME_TO_CATEGORY)
	};
}
#endif

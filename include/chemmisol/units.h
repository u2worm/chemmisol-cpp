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
 * Defines available units.
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

#include <unordered_map>
#include <string>

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

#define UNIT_DEF(NAME, VALUE, CATEGORY) static const double NAME = VALUE;
#define UNIT_NAME_TO_CATEGORY(NAME, VALUE, CATEGORY) {#NAME, Category::CATEGORY},
#define UNIT_NAME_TO_VAL(NAME, VALUE, CATEGORY) {#NAME, NAME},

namespace chemmisol {
	static const double NA = 6.02214076e23;

	enum Category {
		MASS, LENGTH, VOLUME, SURFACE, QUANTITY
	};

	UNITS(UNIT_DEF)
	static const std::unordered_map<std::string, Category> units_category {
		UNITS(UNIT_NAME_TO_CATEGORY)
	};
	static const std::unordered_map<std::string, double> units_name_to_val {
		UNITS(UNIT_NAME_TO_VAL)
	};
}
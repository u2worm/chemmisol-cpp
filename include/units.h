#include <unordered_map>
#include <string>

#define UNITS(F)\
	F(m, 1, LENGTH)\
	F(cm, 1e-2*m, LENGTH)\
	F(nm, 1e-9*m, LENGTH)\
	F(nm2, nm*nm, SURFACE)\
	F(m2, m*m, SURFACE)\
	F(m3, m*m*m, VOLUME)\
	F(cm3, cm*cm*cm, VOLUME)\
	F(l, 1e-3*m3, VOLUME)\
	F(gram, 1, MASS)\
	F(kg, 1e3*gram, MASS)\
	F(mg, 1e-3*gram, MASS)\
	F(mol, 1, QUANTITY)\
	F(entities, 1/NA, QUANTITY)\
	F(u, 1.66053906660e-27*kg, MASS)

#define UNIT_DEF(NAME, VALUE, CATEGORY) static const double NAME = VALUE;
#define UNIT_NAME_TO_CATEGORY(NAME, VALUE, CATEGORY) {#NAME, Category::CATEGORY},
#define UNIT_NAME_TO_VAL(NAME, VALUE, CATEGORY) {#NAME, NAME},

namespace mineral {
	static const double NA = 6.02214076e23;

	enum Operator {
		DIV, TIMES
	};

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

	struct UnitToken {
		bool compound;
		double unit1;
		Category unit1_cat;
		double unit2;
		Category unit2_cat;
		Operator op;

		UnitToken(const std::string& unit1);
		UnitToken(const std::string& unit1, const std::string& unit2, Operator op);
	};

	UnitToken parse(const std::string& unit);
	double convert(double value, const std::string& from, const std::string& to);
	double convert(const UnitToken& from, const UnitToken& to);
}

#include "units.h"
#include <regex>
//#include <iostream>
#include <stdexcept>

namespace mineral {
	const std::regex unit_regex ("(\\w*)(\\/|\\.)?(\\w*)((?:-1)?)");

	UnitToken::UnitToken(const std::string& unit1) :
		compound(false),
		unit1(units_name_to_val.at(unit1)),
		unit1_cat(units_category.at(unit1)) {
		}

	UnitToken::UnitToken(
			const std::string& unit1, const std::string& unit2, Operator op) :
		compound(true),
		unit1(units_name_to_val.at(unit1)),
		unit1_cat(units_category.at(unit1)),
		unit2(units_name_to_val.at(unit2)),
		unit2_cat(units_category.at(unit2)),
		op(op) {
		}

	UnitToken parse(const std::string& unit) {
		std::smatch unit_match;
		std::regex_match(unit, unit_match, unit_regex);
		//for(auto& item : unit_match)
			//std::cout << item << ", ";
		//std::cout << std::endl;
		if(unit_match[1].length() > 0) {
			std::string base_unit = unit_match[1];
			if(unit_match[2].length() == 0)
				return {base_unit};
			if(unit_match[3].length() == 0) {
				throw std::invalid_argument("Invalid operator");
			} else {
				std::string second_unit = unit_match[3];
				if(unit_match[2] == ".") {
					if(unit_match[4] == "")
						return {base_unit, second_unit, TIMES};
					else if(unit_match[4] == "-1")
						return {base_unit, second_unit, DIV};
					else
						throw std::invalid_argument("Invalid modifier");
				} else if(unit_match[2] == "/") {
					if(unit_match[4] == "")
						return {base_unit, second_unit, DIV};
					else if(unit_match[4] == "-1")
						return {base_unit, second_unit, TIMES};
					else
						throw std::invalid_argument("Invalid modifier");
				}
			}
		}
		throw std::invalid_argument("Unrecognized unit: " + unit);
	}

	double convert(double value, const std::string& from, const std::string& to) {
		auto _from = parse(from);
		auto _to = parse(to);
		return value*convert(_from, _to);
	}

	double convert(const UnitToken& from, const UnitToken& to) {
		if(from.compound != to.compound)
			throw std::invalid_argument("Cannot convert between units with different dimensions");
		double result;
		if(from.unit1_cat == to.unit1_cat)
			result = from.unit1/to.unit1;
		else
			throw std::invalid_argument("Cannot convert between units with different categories");

		if(from.compound) {
			if(from.unit2_cat == to.unit2_cat) {
				if(from.op == to.op) {
					switch(from.op) {
						case DIV:
							result = result/(from.unit2/to.unit2);
							break;
						case TIMES:
							result = result*(from.unit2/to.unit2);
					}
				} else {
					throw std::invalid_argument("Cannot convert between units of different nature");
				}
			} else {
				throw std::invalid_argument("Cannot convert between units with different categories");
			}
		}
		return result;
	}
}

#include "gmock/gmock.h"
#include "units.h"

using namespace mineral;
using namespace testing;

TEST(Unit, parse) {
	auto _cm3 = parse("cm3");
	ASSERT_FALSE(_cm3.compound);
	ASSERT_EQ(_cm3.unit1, mineral::cm3);
	ASSERT_EQ(_cm3.unit1_cat, mineral::VOLUME);

	auto _cm3_kg_1 = parse("cm3/kg");
	ASSERT_TRUE(_cm3_kg_1.compound);
	ASSERT_EQ(_cm3_kg_1.unit1, mineral::cm3);
	ASSERT_EQ(_cm3_kg_1.unit1_cat, mineral::VOLUME);
	ASSERT_EQ(_cm3_kg_1.unit2, mineral::kg);
	ASSERT_EQ(_cm3_kg_1.unit2_cat, mineral::MASS);
	auto _cm3_kg_2 = parse("cm3.kg-1");
	ASSERT_TRUE(_cm3_kg_2.compound);
	ASSERT_EQ(_cm3_kg_2.unit1, mineral::cm3);
	ASSERT_EQ(_cm3_kg_2.unit1_cat, mineral::VOLUME);
	ASSERT_EQ(_cm3_kg_2.unit2, mineral::kg);
	ASSERT_EQ(_cm3_kg_2.unit2_cat, mineral::MASS);
}

TEST(Unit, convert) {
	double cm3 = convert(10, "m3", "cm3");
	ASSERT_FLOAT_EQ(cm3, 10*1e6);
}

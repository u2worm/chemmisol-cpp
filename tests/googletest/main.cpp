#include <gtest/gtest.h>
#include "chemmisol/logging.h"

int main(int argc, char *argv[])
{
	::testing::InitGoogleTest(&argc, argv);
	START_EASYLOGGINGPP(argc, argv);
	return RUN_ALL_TESTS();
}

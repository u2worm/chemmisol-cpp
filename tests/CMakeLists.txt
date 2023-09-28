add_subdirectory(googletest)

enable_testing()

add_executable(
  chemmisol-tests
  googletest/main.cpp
  chemmisol/tests/math/linear.cpp
  chemmisol/tests/math/newton.cpp
  chemmisol/tests/math/gauss.cpp
  chemmisol/tests/math/regula_falsi.cpp
  chemmisol/tests/math/homotopy.cpp
  chemmisol/tests/math/arithmetics.cpp
  chemmisol/tests/chemical/system/utils.cpp
  chemmisol/tests/chemical/system/basic.cpp
  chemmisol/tests/chemical/system/nacl_test.cpp
  chemmisol/tests/chemical/system/po4_test.cpp
  chemmisol/tests/chemical/system/agcl_test.cpp
  chemmisol/tests/chemical/system/soh_test.cpp
  chemmisol/tests/chemical/reaction.cpp
)
target_include_directories(
	chemmisol-tests
	PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
	)

target_link_libraries(
  chemmisol-tests
  chemmisol
  GTest::gtest
  GTest::gmock
)

include(GoogleTest)
gtest_discover_tests(chemmisol-tests)

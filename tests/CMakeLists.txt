add_subdirectory(googletest)

enable_testing()

add_executable(
  chemmisol-tests
  chemmisol/math/linear.cpp
  chemmisol/math/newton.cpp
  chemmisol/math/gauss.cpp
  chemmisol/math/regula_falsi.cpp
  chemmisol/chemical/system.cpp
  chemmisol/chemical/experiment.cpp
)
target_link_libraries(
  chemmisol-tests
  chemmisol
  GTest::gtest_main
  GTest::gmock_main
)

include(GoogleTest)
gtest_discover_tests(chemmisol-tests)

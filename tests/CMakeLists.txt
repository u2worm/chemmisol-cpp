add_subdirectory(googletest)

enable_testing()

add_executable(
  chemmisol-tests
  chemmisol/linear.cpp
  chemmisol/newton.cpp
  chemmisol/chemical.cpp
  chemmisol/gauss.cpp
  chemmisol/experiment.cpp
  chemmisol/regula_falsi.cpp
)
target_link_libraries(
  chemmisol-tests
  chemmisol
  GTest::gtest_main
  GTest::gmock_main
)

include(GoogleTest)
gtest_discover_tests(chemmisol-tests)

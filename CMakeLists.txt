cmake_minimum_required(VERSION 3.14)

project(newton-chemical)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "\
${CMAKE_CXX_FLAGS_DEBUG} -std=c++14 \
-Wall -Wextra -Wpedantic -Werror -Wno-error=pedantic \
-Woverloaded-virtual -Wnon-virtual-dtor\
")

add_library(newton-chemical-lib
	src/chemical.cpp
	src/gauss.cpp
	src/experiment.cpp
	)
target_include_directories(newton-chemical-lib PUBLIC include)
add_executable(newton-chemical
	src/main.cpp)
target_link_libraries(newton-chemical newton-chemical-lib)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        a3580180d16923d6d5f488e20b3814608a892f17
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  newton_test
  tests/linear.cpp
  tests/newton.cpp
  tests/chemical.cpp
  tests/gauss.cpp
)
target_link_libraries(
  newton_test
  newton-chemical-lib
  GTest::gtest_main
  GTest::gmock_main
)

include(GoogleTest)
gtest_discover_tests(newton_test)


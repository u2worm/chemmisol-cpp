#cmake_minimum_required(VERSION 3.14)

#project(basic_chemical_system_example)

#set(CMAKE_CXX_STANDARD 14)

add_executable(basic_chemical_system_example main.cpp)
target_link_libraries(basic_chemical_system_example chemmisol)
target_include_directories(basic_chemical_system_example PUBLIC ../include)

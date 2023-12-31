cmake_minimum_required(VERSION 3.14)

project(chemmisol VERSION 0.1)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)

set(CHEMMISOL_MEMORY_SANITIZE NO CACHE BOOL "Enables -fsanitize=address.")

set(CMAKE_CXX_FLAGS_DEBUG "\
${CMAKE_CXX_FLAGS_DEBUG} -std=c++14 \
-Wall -Wextra -Wpedantic -Wno-error=pedantic \
-Woverloaded-virtual -Wnon-virtual-dtor \
-Wno-range-loop-construct \
")
if(CHEMMISOL_MEMORY_SANITIZE)
	message(STATUS "Enabling -fsanitize=address")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
endif()

set(CHEMMISOL_LOCAL_EASYLOGGING "" CACHE PATH
	"Path to easylogging++ sources. If not specified, sources will be downloaded at configure time."
)

if(CHEMMISOL_LOCAL_EASYLOGGING)
	set(EASYLOGGING_SRC_PATH ${CHEMMISOL_LOCAL_EASYLOGGING})
else()
	# Download easylogging++ sources
	configure_file(easylogging++/CMakeLists.txt.in
		easylogging-download/CMakeLists.txt)
	execute_process(COMMAND ${CMAKE_COMMAND} .
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/easylogging-download)
	execute_process(COMMAND ${CMAKE_COMMAND} --build .
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/easylogging-download) 
	set(EASYLOGGING_SRC_PATH
		${CMAKE_CURRENT_BINARY_DIR}/_deps/easyloggingcpp-src/src
		)
endif()

find_package(OpenMP)

add_library(chemmisol
	src/chemmisol/units.cpp
	src/chemmisol/math/gauss.cpp
	src/chemmisol/chemical/species.cpp
	src/chemmisol/chemical/reaction.cpp
	src/chemmisol/chemical/solver.cpp
	src/chemmisol/chemical/system.cpp
	src/chemmisol/logging.cpp
	${EASYLOGGING_SRC_PATH}/easylogging++.cc
	)
if(OpenMP_CXX_FOUND)
	target_compile_definitions(chemmisol PUBLIC
		CHEMMISOL_OPENMP
		ELPP_THREAD_SAFE
		)
	target_link_libraries(chemmisol PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(chemmisol PUBLIC
	$<BUILD_INTERFACE:${EASYLOGGING_SRC_PATH}>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
)

target_compile_definitions(chemmisol PUBLIC
	$<$<CONFIG:Release>:ELPP_DISABLE_TRACE_LOGS>)
target_compile_definitions(chemmisol PUBLIC ELPP_STL_LOGGING)
target_compile_definitions(chemmisol PUBLIC ELPP_LOG_STD_ARRAY)

install(TARGETS chemmisol EXPORT ChemmisolTargets
	LIBRARY DESTINATION lib
	INCLUDES DESTINATION include
	)

install(
	DIRECTORY include
	DESTINATION .
)

install(
	FILES ${EASYLOGGING_SRC_PATH}/easylogging++.h
	DESTINATION ./include/chemmisol
	)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/ChemmisolConfigVersion.cmake"
	COMPATIBILITY SameMajorVersion
	)

export(EXPORT ChemmisolTargets
	FILE "${CMAKE_CURRENT_BINARY_DIR}/ChemmisolTargets.cmake"
	NAMESPACE chemmisol::
	)
add_library(chemmisol::chemmisol ALIAS chemmisol)
configure_file(cmake/ChemmisolConfig.cmake
	"${CMAKE_CURRENT_BINARY_DIR}/ChemmisolConfig.cmake"
	COPYONLY
	)

set(ConfigPackageLocation lib/cmake/Chemmisol)
install(EXPORT ChemmisolTargets
	FILE
	ChemmisolTargets.cmake
	NAMESPACE
	chemmisol::
	DESTINATION
	${ConfigPackageLocation}
	)
install(
	FILES
	cmake/ChemmisolConfig.cmake
	"${CMAKE_CURRENT_BINARY_DIR}/ChemmisolConfigVersion.cmake"
	DESTINATION
	${ConfigPackageLocation}
	)

add_subdirectory(tests)
add_subdirectory(examples)

# Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
	file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/docs)
	set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/docs)
	#set(DOXYGEN_MACRO_EXPANSION "YES")
	set(DOXYGEN_HTML_OUTPUT ".")
	set(DOXYGEN_SORT_MEMBER_DOCS "NO")
	set(DOXYGEN_HAVE_DOT "YES")
	#set(DOXYGEN_TEMPLATE_RELATIONS "YES")
	set(DOXYGEN_ALIASES
		[[addComponent=@ref addComponent(const std::string&,Phase,double) \"addComponent()\"]]
		[[fixComponent=@ref fixComponent(const std::string&,Phase,double) \"fixComponent()\"]]
		[[addSolvent=@ref addSolvent(const std::string&) \"addSolvent()\"]]
		[[addReaction=@ref addReaction(const std::string&, double, const std::vector<Reagent>&) \"addReaction()\"]]
		[[massConservationLaw=@ref ChemicalSystem::massConservationLaw(const std::vector<double>&, std::vector<double>&) const \"ChemicalSystem::massConservationLaw()\"]]
		[[distanceToEquilibrium=@ref ChemicalSystem::distanceToEquilibrium(const std::vector<double>&, const Reaction&) const \"ChemicalSystem::distanceToEquilibrium()\"]]
		)
	set(DOXYGEN_EXAMPLE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/examples)
	set(DOXYGEN_EXAMPLE_PATTERNS "*.cpp")
	doxygen_add_docs(doc
		${CMAKE_CURRENT_SOURCE_DIR}/include
		)
endif()

cmake_minimum_required(VERSION 3.14)

project(chemmisol VERSION 0.1)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "\
${CMAKE_CXX_FLAGS_DEBUG} -std=c++14 \
-Wall -Wextra -Wpedantic -Werror -Wno-error=pedantic \
-Woverloaded-virtual -Wnon-virtual-dtor \
-fsanitize=address \
")

add_library(chemmisol
	src/chemmisol/units.cpp
	src/chemmisol/chemical.cpp
	src/chemmisol/gauss.cpp
	src/chemmisol/experiment.cpp
	)

#target_include_directories(chemmisol PUBLIC include)
target_include_directories(chemmisol PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

install(TARGETS chemmisol EXPORT ChemmisolTargets
	LIBRARY DESTINATION lib
	INCLUDES DESTINATION include
	)

install(
	DIRECTORY include
	DESTINATION .
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
	doxygen_add_docs(doc
		${CMAKE_CURRENT_SOURCE_DIR}/include
		)
endif()

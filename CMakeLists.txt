cmake_minimum_required (VERSION 2.8)

project (barotropic-model CXX)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-c++11-extensions")

# ------------------------------------------------------------------------------
# external libraries
# GeoMTK
add_subdirectory ("${PROJECT_SOURCE_DIR}/external/geomtk")
include_directories (
    ${GEOMTK_INCLUDE_DIRS}
)
link_directories (
    ${GEOMTK_LIBRARY_DIRS}
)

# ------------------------------------------------------------------------------
# executable targets
add_executable (main.exe
    "${PROJECT_SOURCE_DIR}/src/commons.h"
    "${PROJECT_SOURCE_DIR}/src/BarotropicTestCase.h"
    "${PROJECT_SOURCE_DIR}/src/RossbyHaurwitzTestCase.h"
    "${PROJECT_SOURCE_DIR}/src/RossbyHaurwitzTestCase.cpp"
    "${PROJECT_SOURCE_DIR}/src/BarotropicModel.h"
	"${PROJECT_SOURCE_DIR}/src/BarotropicModel.cpp"
	"${PROJECT_SOURCE_DIR}/src/main.cpp"
)
target_link_libraries (main.exe geomtk)

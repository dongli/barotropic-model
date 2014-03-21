cmake_minimum_required (VERSION 2.8)

project (barotropic-model CXX)

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
add_executable (run_model
    "${PROJECT_SOURCE_DIR}/src/commons.h"
    "${PROJECT_SOURCE_DIR}/src/BarotropicTestCase.h"
    "${PROJECT_SOURCE_DIR}/src/RossbyHaurwitzTestCase.h"
    "${PROJECT_SOURCE_DIR}/src/RossbyHaurwitzTestCase.cpp"
    "${PROJECT_SOURCE_DIR}/src/BarotropicModel.h"
    "${PROJECT_SOURCE_DIR}/src/BarotropicModel_A_ImplicitMidpoint.h"
	"${PROJECT_SOURCE_DIR}/src/BarotropicModel_A_ImplicitMidpoint.cpp"
	"${PROJECT_SOURCE_DIR}/src/run_model.cpp"
)
target_link_libraries (run_model geomtk)

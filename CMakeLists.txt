cmake_minimum_required (VERSION 2.8)

# ------------------------------------------------------------------------------
# check if there is upper project
if (DEFINED PROJECT_NAME)
    message ("@@ barotropic-model is used as a submodule.")
    set (use_as_submodule TRUE)
endif ()

project (barotropic-model CXX)

if (NOT use_as_submodule)
    # --------------------------------------------------------------------------
    # external libraries
    # GeoMTK
    add_subdirectory ("${PROJECT_SOURCE_DIR}/external/geomtk")
    include_directories (
        ${GEOMTK_INCLUDE_DIRS}
    )
    link_directories (
        ${GEOMTK_LIBRARY_DIRS}
    )
endif ()

# ------------------------------------------------------------------------------
# collect sources and headers
set (sources
    "${PROJECT_SOURCE_DIR}/src/barotropic_model_commons.h"
    "${PROJECT_SOURCE_DIR}/src/GeostrophicRelation.h"
    "${PROJECT_SOURCE_DIR}/src/GeostrophicRelation.cpp"
    "${PROJECT_SOURCE_DIR}/src/BarotropicTestCase.h"
    "${PROJECT_SOURCE_DIR}/src/RossbyHaurwitzTestCase.h"
    "${PROJECT_SOURCE_DIR}/src/RossbyHaurwitzTestCase.cpp"
    "${PROJECT_SOURCE_DIR}/src/ToyTestCase.h"
    "${PROJECT_SOURCE_DIR}/src/ToyTestCase.cpp"
    "${PROJECT_SOURCE_DIR}/src/BarotropicModel.h"
    "${PROJECT_SOURCE_DIR}/src/BarotropicModel_A_ImplicitMidpoint.h"
    "${PROJECT_SOURCE_DIR}/src/BarotropicModel_A_ImplicitMidpoint.cpp"
    "${PROJECT_SOURCE_DIR}/src/BarotropicModel_C_ImplicitMidpoint.h"
    "${PROJECT_SOURCE_DIR}/src/BarotropicModel_C_ImplicitMidpoint.cpp"
)
if (use_as_submodule)
    set (BAROTROPIC_MODEL_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/src" PARENT_SCOPE)
    set (BAROTROPIC_MODEL_LIBRARY_DIRS "${PROJECT_SOURCE_DIR}/src" PARENT_SCOPE)
endif ()

# ------------------------------------------------------------------------------
# library targets
add_library (barotropic-model STATIC ${sources})
add_dependencies (barotropic-model geomtk)

# ------------------------------------------------------------------------------
# executable targets
add_executable (run_model
	"${PROJECT_SOURCE_DIR}/src/run_model.cpp"
)
target_link_libraries (run_model geomtk barotropic-model)

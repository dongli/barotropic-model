cmake_minimum_required (VERSION 2.8)

# Check if this project is used as a submodule.
if (DEFINED PROJECT_NAME)
    message ("@@ BAROTROPIC-MODEL is used as a submodule.")
    set (use_as_submodule TRUE)
else ()
    option (FLAG_OPENMP "Turn OpenMP compiler flag ON or OFF" OFF)
    option (FLAG_SHARED "Turn building shared libraries ON of OFF" OFF)

    if (FLAG_OPENMP)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    else ()
        message ("@@ LASM does not use OpenMP compiler flag.")
    endif ()
    if (FLAG_SHARED)
        set (shared_or_static SHARED)
    else ()
        set (shared_or_static STATIC)
    endif ()

    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif ()

project (barotropic-model CXX)

if (NOT use_as_submodule)
    # Add external libraries.
    # GEOMTK 
    add_subdirectory ("${PROJECT_SOURCE_DIR}/external/geomtk")
    include_directories (${GEOMTK_INCLUDE_DIRS})
endif ()

# Collect sources and headers.
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

# Record the source directories into <PROJECT_NAME>_INCLUDE_DIRS for upper
# project to use.
if (use_as_submodule)
    set (BAROTROPIC_MODEL_INCLUDE_DIRS
        "${PROJECT_SOURCE_DIR}/src"
        PARENT_SCOPE
    )
endif ()

# Add library targets.
add_library (barotropic-model ${shared_or_static} ${sources})
target_link_libraries (barotropic-model geomtk)
add_dependencies (barotropic-model geomtk)

# Add executable targets.
add_executable (run_model
    "${PROJECT_SOURCE_DIR}/src/run_model.cpp"
)
target_link_libraries (run_model
    geomtk
    barotropic-model
)

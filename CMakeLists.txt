# CMake 3.18 required to compile Cuda code with c++17. If only an earlier version of
# Cmake is available we need to set all the standards to c++14 instead (if using Cuda).
cmake_minimum_required(VERSION 3.18)
set(CPP_VERSION "17")

# set project name and info
project(amdem
        VERSION 1.0
        DESCRIPTION "A DEM code to model powder deposition for additive manufacturing"
        LANGUAGES CXX)

project(amdem-lib
        VERSION 1.0
        DESCRIPTION "A library of DEM helper functions/calsses"
        LANGUAGES CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/TPLs" ${CMAKE_MODULE_PATH})

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR " You cannot build in a source directory (or any directory with a CMakeLists.txt file). "
        "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles." )
endif()

# Set default to Debug for now... can change this once everything is working
set(default_build_type "Debug") 
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    # set the available options for the cmake-gui 
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 "Debug" "Release" "MinSeizeRel" "RelWithDebInfo")
endif()

# see full compile and linking commands for compilation debugging
set(CMAKE_VERBOSE_MAKEFILE ON)

#------------------------------------------------------------------------------
# Setup BLT
#------------------------------------------------------------------------------
if (DEFINED BLT_SOURCE_DIR)
    # Support having a shared BLT outside of the repository if given a BLT_SOURCE_DIR
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR "Given BLT_SOURCE_DIR does not contain SetupBLT.cmake")
    endif()
else()
    # Use internal BLT if no BLT_SOURCE_DIR is given
    set(BLT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/TPLs/blt" CACHE PATH "")
    if (NOT EXISTS ${BLT_SOURCE_DIR}/SetupBLT.cmake)
        message(FATAL_ERROR
            "The BLT git submodule is not present. "
            "Either run the following two commands in your git repository: \n"
            "    git submodule init\n"
            "    git submodule update\n"
            "Or add -DBLT_SOURCE_DIR=/path/to/blt to your CMake command." )
    endif()
endif()

if (NOT BLT_CXX_STD)
    set(BLT_CXX_STD "c++${CPP_VERSION}" CACHE STRING "")
endif()

include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

#------------------------------------------------------------------------------
# Linking to Kokkos, et al.
# N.B. We set the KOKKOS_DIR variable in the appropriate host-config file as
# different Kokkos compilations are required depending on if we link to CUDA or
# other build configurations 
#------------------------------------------------------------------------------
if(NOT EXISTS ${KOKKOS_DIR}/KokkosConfig.cmake)
    message(FATAL_ERROR "Missing required file \"KOKKOS_DIR/KokkosConfig.cmake\" associated with an installed Kokkos, "
        "ensure directory is correct in host-config file")
endif()
find_package(Kokkos REQUIRED
             NO_DEFAULT_PATH
             PATHS ${KOKKOS_DIR})

if (Kokkos_ENABLE_CUDA)
    # fatal error if not enabled
    kokkos_check(OPTIONS CUDA_LAMBDA)
    set(amdem_dependencies ${amdem_dependencies} openmp cuda)
else()
    set(amdem_dependencies ${amdem_dependencies} openmp)
endif()

set(amdem_dependencies ${amdem_dependencies} Kokkos::kokkos)

#==========================================================================
# if(NOT EXISTS ${KOKKOS_KERNELS_DIR}/KokkosKernelsConfig.cmake)
#     message(FATAL_ERROR "Missing required file \"${KOKKOS_KERNELS_DIR}/KokkosKernelsConfig.cmake\" associated with an installed KokkosKernels, ensure directory is correct in host-config file")
# endif()
# find_package(KokkosKernels REQUIRED
#              NO_DEFAULT_PATH
#              PATHS ${KOKKOS_KERNELS_DIR})
# 
# set(amdem_dependencies ${amdem_dependencies} Kokkos::kokkoskernels)


#------------------------------------------------------------------------------
# Link to fmt library
#------------------------------------------------------------------------------
SET(FMT_ROOT "${CMAKE_CURRENT_LIST_DIR}/TPLs/fmt/build/fmt_install" CACHE PATH "")
SET(FMT_DIR "${FMT_ROOT}/lib/cmake/fmt" CACHE PATH "")
if(NOT EXISTS ${FMT_DIR}/fmt-config.cmake)
    message(FATAL_ERROR "Missing required file \"FMT_DIR/fmt-config.cmake\" associated with an installed fmt, "
        "ensure directory is correct in CMakeLists.txt or that the TPLs/fmt git submodule is installed.")
endif()

find_package(fmt REQUIRED
             NO_DEFAULT_PATH
             PATHS ${FMT_DIR})

set(amdem_dependencies ${amdem_dependencies} fmt::fmt)

#------------------------------------------------------------------------------
# lean-vtk library for writing VTK files
#------------------------------------------------------------------------------
if(NOT DEFINED LEAN-VTK_DIR)
    SET(LEAN-VTK_DIR "${CMAKE_CURRENT_LIST_DIR}/TPLs/lean-vtk/build/lean-vtk_install" CACHE PATH "")
    SET(LEAN-VTK_LIB_DIR "${LEAN-VTK_DIR}/lib/LeanVTK-1.0" CACHE PATH "")
endif()
if(NOT EXISTS ${LEAN-VTK_DIR}/lib/LeanVTK-1.0/LeanVTK-config.cmake)
    message(FATAL_ERROR "Missing required 'LEAN-VTK_DIR' variable pointing to an installed lean-vtk")
endif()
if(NOT DEFINED LEAN-VTK_INCLUDE_DIR)
    SET(LEAN-VTK_INCLUDE_DIR "${LEAN-VTK_DIR}/include/LeanVTK-1.0" CACHE PATH "")
endif()

find_library(LEAN-VTK_LIBRARY NAMES LeanVTK
             PATHS ${LEAN-VTK_LIB_DIR}
             NO_DEFAULT_PATH
             NO_CMAKE_ENVIRONMENT_PATH
             NO_CMAKE_PATH
             NO_SYSTEM_ENVIRONMENT_PATH
             NO_CMAKE_SYSTEM_PATH)

blt_import_library(NAME       LeanVTK
                   INCLUDES   ${LEAN-VTK_INCLUDE_DIR}
                   LIBRARIES  ${LEAN-VTK_LIBRARY}
                   EXPORTABLE ON)

set(amdem_dependencies ${amdem_dependencies} LeanVTK)

#------------------------------------------------------------------------------
# Include the header only CLI11 library
#------------------------------------------------------------------------------
set(amdem_include_dirs ${amdem_include_dirs} "${CMAKE_CURRENT_LIST_DIR}/TPLs/CLI11")


##------------------------------------------------------------------------------
## Create actual library and exe
##------------------------------------------------------------------------------

# Add directories and compile commands
add_subdirectory(src)

# add library and exe (which just has main then links to the lib)
blt_add_library(NAME       amdem-lib
                SOURCES    ${amdem-lib_source_files}
                HEADERS    ${amdem-lib_header_files}
                INCLUDES   ${amdem_include_dirs}
                DEFINES    ${amdem_defines}
                DEPENDS_ON ${amdem_dependencies})

blt_add_executable(NAME       amdem
                   SOURCES    ${amdem-exe_source_files}
                   INCLUDES   ${amdem_include_dirs}
                   DEFINES    ${amdem_defines}
                   DEPENDS_ON ${amdem_dependencies} amdem-lib)


# the next line ensures -std=c++xx instead of -std=g++xx
set_target_properties(amdem PROPERTIES CXX_EXTENSIONS OFF)
set_target_properties(amdem-lib PROPERTIES CXX_EXTENSIONS OFF)


##------------------------------------------------------------------------------
## Set up unit tests
##------------------------------------------------------------------------------
if(ENABLE_TESTS) # on by default with BLT

    add_subdirectory(tests)
    set(testing_dependencies amdem-lib gtest)
    set(test_include_dirs ${amdem_include_dirs} tests src)

    # 1. test_particles
    blt_add_executable(NAME         test_particles
                       SOURCES      ${test_particles_source_files}
                       INCLUDES     ${test_include_dirs}
                       DEPENDS_ON   ${testing_dependencies})
    # experiment to see compiler requirements...
    # blt_print_target_properties(TARGET test_particles)

    # the command option can be used to pass args to the test exe
    blt_add_test(NAME test_particles
                 COMMAND test_particles)

    # 2. test_bins
    blt_add_executable(NAME         test_bins
                       SOURCES      ${test_bins_source_files}
                       INCLUDES     ${test_include_dirs}
                       DEPENDS_ON   ${testing_dependencies})

    # the command option can be used to pass args to the test exe
    blt_add_test(NAME test_bins
                 COMMAND test_bins)


endif()

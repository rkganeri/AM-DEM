#------------------------------------------------------------------------------
# host-config file for Mac OS 10.13.6 (High Sierra) built with the Clang11
# compiler 
#------------------------------------------------------------------------------
#
# This file provides CMake with paths / details for:
#  C/C++:   Clang 11.1.0
#
#------------------------------------------------------------------------------
#---------------------------------------
# Definitions
#---------------------------------------
set(amdem_defines ${amdem_defines} AMDEM_HOST_ONLY CACHE PATH "")

#---------------------------------------
# Compilers
#---------------------------------------
set(CMAKE_C_COMPILER "/opt/local/bin/clang-mp-11" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/opt/local/bin/clang++-mp-11" CACHE PATH "")

#---------------------------------------
# MPI
#---------------------------------------
set(ENABLE_MPI OFF CACHE BOOL "")

#------------------------------------------------------------------------------
# OpenMP support
#------------------------------------------------------------------------------
set(ENABLE_OPENMP ON CACHE BOOL "")

#------------------------------------------------------------------------------
# Cuda
#------------------------------------------------------------------------------
set(ENABLE_CUDA OFF CACHE BOOL "")
# Cuda requires setting the right flags and linking with nvcc. Unfortunately,
# my 2011 MBP laptop does not have a compute compatible GPU.

#------------------------------------------------------------------------------
# Third Party Library Dependencies
# We only set lib paths here that have separate installs depending on build type
#------------------------------------------------------------------------------
SET(KOKKOS_ROOT "${CMAKE_CURRENT_LIST_DIR}/../TPLs/kokkos/build_clang11_omp_debug/kokkos_install" CACHE PATH "")
SET(KOKKOS_DIR "${KOKKOS_ROOT}/lib/cmake/Kokkos" CACHE PATH "")

#SET(KOKKOS_KERNELS_ROOT "${CMAKE_CURRENT_LIST_DIR}/../TPLs/kokkos-kernels/build_clang11_omp_debug/kokkos-kernels_install" CACHE PATH "")
#SET(KOKKOS_KERNELS_DIR "${KOKKOS_KERNELS_ROOT}/lib/cmake/KokkosKernels" CACHE PATH "")


#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

#set(CLANGFORMAT_EXECUTABLE "/opt/local/bin/clang-format-mp-11" CACHE PATH "")

#set(CLANGTIDY_EXECUTABLE "/opt/local/bin/clang-tidy-mp-11" CACHE PATH "")

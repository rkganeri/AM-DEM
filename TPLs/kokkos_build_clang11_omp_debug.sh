# command from https://github.com/kokkos/kokkos/blob/master/BUILD.md

cmake \
  -DCMAKE_CXX_COMPILER=/opt/local/bin/clang++-mp-11 \
  -DKokkos_CXX_STANDARD="17" \
  -DCMAKE_INSTALL_PREFIX="kokkos_install" \
  -DKokkos_ENABLE_SERIAL=On \
  -DKokkos_ENABLE_OPENMP=On \
  -DKokkos_ENABLE_DEBUG=ON \
  -DKokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON \
  ..
#  -DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=On \
#  -DKokkos_ENABLE_TESTS=On \
#  -DKokkos_ENABLE_EXAMPLES=On 
#  -DCMAKE_INSTALL_PREFIX=$path_to_install \
#  -DKokkos_ENABLE_HWLOC=On \
#  -DKokkos_HWLOC_DIR=$path_to_hwloc
#  -DCMAKE_CXX_COMPILER=$(pwd)/../bin/nvcc_wrapper \
#  -DCMAKE_CXX_COMPILER=g++ \

make install

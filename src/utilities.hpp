#ifndef AMDEM_UTILITIES_HPP
#define AMDEM_UTILITIES_HPP

#include <cstdlib>  
#include <cmath> 

#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

namespace amdem {

namespace utilities {

// declare an inline host or device function to calculate the norm of the difference
// between 2 vectors
KOKKOS_INLINE_FUNCTION double diffNorm(double* a, double* b, int len=3) {
    double val = 0.0;
    for (int i=0; i<len; i++) {
        val += pow(a[i]-b[i], 2.0);
    }
    return sqrt(val);
}
    



            xd = h_coords(i,0) - h_coords(j,0);
            yd = h_coords(i,1) - h_coords(j,1);
            zd = h_coords(i,2) - h_coords(j,2);
            dist = sqrt(pow(xd,2.) + pow(yd,2.) + pow(zd,2.));
    


// A Functor for generating uint64_t random numbers templated on the
// GeneratorPool type, see for details
// https://github.com/kokkos/kokkos/blob/master/example/tutorial/Algorithms/01_random_numbers/random_numbers.cpp
template <class GeneratorPool>
struct GenerateRandom {
  // Output View for the random numbers
  Kokkos::View<uint64_t*> vals_;

  // The GeneratorPool
  GeneratorPool rand_pool_;

  int samples_;

  // Initialize all members
  GenerateRandom(Kokkos::View<uint64_t*> vals, GeneratorPool rand_pool,
                  int samples)
      : vals_(vals), rand_pool_(rand_pool), samples_(samples) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    // Get a random number state from the pool for the active thread
    typename GeneratorPool::generator_type rand_gen = rand_pool_.get_state();

    // Draw samples numbers from the pool as urand64 between 0 and
    // rand_pool.MAX_URAND64 Note there are function calls to get other type of
    // scalars, and also to specify Ranges or get a normal distributed float.
    for (int k = 0; k < samples; k++)
      vals_(i * samples_ + k) = rand_gen.urand64();

    // Give the state back, which will allow another thread to acquire it
    rand_pool_.free_state(rand_gen);
  }
};

} // namespace utilities

} // namespace amdem

#endif // AMDEM_UTILITIES_HPP




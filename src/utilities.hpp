#ifndef AMDEM_UTILITIES_HPP
#define AMDEM_UTILITIES_HPP

#include <cstdlib>  
#include <cmath> 
#include <type_traits>
#include <typeinfo>
#ifndef _MSC_VER
#   include <cxxabi.h>
#endif
#include <memory>
#include <string>
#include <cstdlib>

#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

namespace amdem {

namespace utilities {

// this guy is a handy template for printing out some of those weird type when you can't use auto
// use would be "std::string type = type_name<decltype(var)>();" to see the typename of var
// taken from https://stackoverflow.com/questions/81870/is-it-possible-to-print-a-variables-type-in-standard-c
template <class T> std::string type_name()
{
    typedef typename std::remove_reference<T>::type TR;
    std::unique_ptr<char, void(*)(void*)> own
           (
#ifndef _MSC_VER
                abi::__cxa_demangle(typeid(TR).name(), nullptr,
                                           nullptr, nullptr),
#else
                nullptr,
#endif
                std::free
           );
    std::string r = own != nullptr ? own.get() : typeid(TR).name();
    if (std::is_const<TR>::value)
        r += " const";
    if (std::is_volatile<TR>::value)
        r += " volatile";
    if (std::is_lvalue_reference<T>::value)
        r += "&";
    else if (std::is_rvalue_reference<T>::value)
        r += "&&";
    return r;
}


// declare an inline host or device function to calculate the 2-norm of a vector 
// 1st version receives a c-arrays
KOKKOS_INLINE_FUNCTION double norm2(double* a, int len=3) {
    double val = 0.0;
    for (int i=0; i<len; i++) {
        val += pow(a[i], 2.0);
    }
    return sqrt(val);
}
// declare a 2nd version in case we pass in components directly (useful when dealing with views)
KOKKOS_INLINE_FUNCTION double norm2(double a0, double a1, double a2) {
    return sqrt(pow(a0,2.0) + pow(a1,2.0) + pow(a2,2.0));
}


// inline function to calculate the dot product of 2 vectors (takes in c-arrays)
KOKKOS_INLINE_FUNCTION double dotProduct(double* a, double* b, int len=3) {
    double val = 0.0;
    for (int i=0; i<len; i++) {
        val += a[i]*b[i];
    }
    return val;
}


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
    for (int k = 0; k < samples_; k++)
      vals_(i * samples_ + k) = rand_gen.urand64();

    // Give the state back, which will allow another thread to acquire it
    rand_pool_.free_state(rand_gen);
  }
};

} // namespace utilities

} // namespace amdem

#endif // AMDEM_UTILITIES_HPP




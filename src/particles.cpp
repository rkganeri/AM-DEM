#include <string>   
#include <ctime>

#include "Kokkos_Core.hpp"
#include "Kokkos_Random.hpp"

#include "particles.hpp"
#include "utilities.hpp"

namespace amdem {

Particles::Particles(const int num_particles)
    : num_particles_(num_particles),
      radius_("radius",num_particles),
      mass_("mass",num_particles),
      volume_("volume",num_particles),
      vn_("vn",num_particles),  // (num_particles,3)
      vnp1_("vnp1",num_particles),  // (num_particles,3)
      coords_("coordsn",num_particles),  // (num_particles,3)
      coordsnp1_("coordsnp1",num_particles),  // (num_particles,3)
{   }  

// our static method to actually instantiate this singleton-esque object
Particles& Particles::GetInstance(const int num_particles) {
    static Particles* particles = new Particles(num_particles);
    return *particles;
}

void Particles::initLocation(GlobalSettings global_settings) {

    // unpack local vars
    const double xmin = 0.0;
    const double xmax = global_settings.length_;
    const double ymin = 0.0;
    const double ymax = global_settings.width_;
    const double zmin = 0.0;
    const double zmax = global_settings.height_;

    const double mean_rad = global_settings.mean_rad_;
    const double stdev_rad = global_settings.stdev_rad_;
    const double min_rad = global_settings.min_rad_;
    const double max_rad = global_settings.max_rad_;

    // first task is to generate random numbers. We do this in parallel using
    // built-in Kokkos helper classes, for details see 
    // https://github.com/kokkos/kokkos/blob/master/example/tutorial/Algorithms/01_random_numbers/random_numbers.cpp
    static_cast<uint64_t> rand_seed = time(NULL)
    Kokkos::Random_XorShift64_Pool<> rand_pool64(rand_seed);

    Kokkos::View<uint64_t*> rand_nums("rand_nums",num_particles_*3);

    // the Kokkos::parallel_for loop construct will default to OpenMP threading if compiled CPU-only
    // or create a Cuda/HIP device function call (through the GenerateRandom functor) if compiled with Cuda or HIP
    Kokkos::parallel_for("generate_random_nums",num_particles_,
                         utilities::GenerateRandom<Kokkos::Random_XorShift64_Pool<> >(
                         rand_nums, rand_pool64, 3));
    
    
    


    // the Kokkos::parallel_for loop construct will default to OpenMP threading if compiled CPU-only
    // or create a Cuda/HIP device function call (through lambda functions) if compiled with Cuda or HIP
    // the lambda defaults to capture by value, which is what we want since Kokkos Views are treated
    // as pointers
    Kokkos::parallel_for("init_DEM_coords", num_particles_, KOKKOS_LAMBDA(int i) {

        // first task is to generate a random number


    
    

        const double length_;
        const double width_;
        const double height_;

        // particle size distribution
        const double mean_rad_;
        const double stdev_rad_;
        const double min_rad_;
        const double max_rad_;
    

}






} // namespace amdem

